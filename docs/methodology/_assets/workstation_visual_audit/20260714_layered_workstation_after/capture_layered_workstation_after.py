#!/usr/bin/env python3
"""Full-scope Chromium audit for the layered reconstruction workstation.

The runner is deliberately read-only with respect to the eight source HTML
files. It captures screenshot evidence before DOM interpretation, tests the
offline interaction surface, and exits non-zero when any hard gate fails.

Exit codes:
    0: every hard gate passed
    1: the audit completed but one or more gates failed
    2: runner/configuration failure
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import struct
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple
from urllib.parse import unquote, urlparse

try:
    from playwright.sync_api import Browser, Locator, Page, sync_playwright
except ImportError as exc:
    Browser = Locator = Page = Any  # type: ignore[misc,assignment]
    sync_playwright = None  # type: ignore[assignment]
    PLAYWRIGHT_IMPORT_ERROR: Optional[str] = str(exc)
else:
    PLAYWRIGHT_IMPORT_ERROR = None


SCRIPT_PATH = Path(__file__).resolve()
DEFAULT_OUTPUT_DIR = SCRIPT_PATH.parent
DEFAULT_INPUT_DIR = SCRIPT_PATH.parents[2] / "layered_workstation"
DEFAULT_METRICS = DEFAULT_OUTPUT_DIR / "metrics.json"
DEFAULT_BEFORE_METRICS = SCRIPT_PATH.parent.parent / "20260714_layered_workstation_before" / "metrics.json"
SCREENSHOT_DIR_NAME = "screenshots"

INDEX_FILE = "index.html"
SAMPLE_FILES = [
    "HCC1395.html",
    "COLO829.html",
    "H1437.html",
    "H2009.html",
    "HCC1395_DORADO.html",
    "HCC1937.html",
    "HCC1954.html",
]
ALL_FILES = [INDEX_FILE, *SAMPLE_FILES]
VIEWPORTS = {
    "desktop": {"width": 1440, "height": 1000},
    "mobile": {"width": 390, "height": 844},
    "narrow": {"width": 320, "height": 720},
}
INTERACTION_VIEWPORTS = {"desktop", "mobile"}
FORBIDDEN_VISIBLE = [
    "posterior",
    "softmax",
    "單 clone",
    "中間 clone",
    "CCF pigeonhole",
    "L3 PENDING",
    "no-germline",
]
EXPECTED_AGGREGATE = {
    "W_tree": 51815,
    "W_primary": 50215,
    "complete": 42240,
    "incomplete": 7975,
    "topology_exact": 11582,
    "topology_shape": 10737,
    "topology_multiple": 19921,
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path) -> Dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path),
        "bytes": stat.st_size,
        "mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
        "sha256": sha256(path),
    }


def png_metadata(path: Path) -> Dict[str, Any]:
    result: Dict[str, Any] = {"path": str(path), "exists": path.is_file()}
    if not path.is_file():
        return result
    result["bytes"] = path.stat().st_size
    try:
        with path.open("rb") as handle:
            header = handle.read(24)
        if len(header) >= 24 and header[:8] == b"\x89PNG\r\n\x1a\n":
            result["width"], result["height"] = struct.unpack(">II", header[16:24])
    except OSError as exc:
        result["metadata_error"] = str(exc)
    return result


def slug(value: str) -> str:
    value = re.sub(r"\.html$", "", value, flags=re.IGNORECASE)
    return re.sub(r"[^0-9A-Za-z_-]+", "_", value).strip("_").lower()


def build_screenshot_plan() -> Dict[Tuple[str, str, str], Dict[str, Any]]:
    plan: Dict[Tuple[str, str, str], Dict[str, Any]] = {}
    step = 1
    for viewport_name in VIEWPORTS:
        key = (INDEX_FILE, viewport_name, "full_page")
        plan[key] = {
            "step": step,
            "description": f"cohort index · {viewport_name} · initial full page",
            "filename": f"{step:02d}_{viewport_name}_index_full_page.png",
        }
        step += 1
    for sample_file in SAMPLE_FILES:
        for viewport_name in VIEWPORTS:
            kinds = [("full_page", "initial full page")]
            if viewport_name in INTERACTION_VIEWPORTS:
                kinds.extend(
                    [
                        ("region_view", "selected incomplete 8-site region"),
                        ("candidate_network", "selected candidate network"),
                    ]
                )
                if sample_file == "HCC1395_DORADO.html" and viewport_name == "desktop":
                    kinds.append(("no_primary_region", "no-primary auxiliary/control evidence caption"))
            for kind, description in kinds:
                key = (sample_file, viewport_name, kind)
                plan[key] = {
                    "step": step,
                    "description": f"{sample_file[:-5]} · {viewport_name} · {description}",
                    "filename": f"{step:02d}_{viewport_name}_{slug(sample_file)}_{kind}.png",
                }
                step += 1
    if step != 54:
        raise AssertionError(f"Expected 53 screenshot steps, got {step - 1}")
    return plan


SCREENSHOT_PLAN = build_screenshot_plan()


def add_check(
    checks: List[Dict[str, Any]],
    name: str,
    passed: bool,
    expected: Any,
    actual: Any,
    *,
    page_name: Optional[str] = None,
    viewport: Optional[str] = None,
    severity: str = "error",
) -> None:
    item: Dict[str, Any] = {
        "name": name,
        "status": "pass" if passed else "fail",
        "severity": severity,
        "expected": expected,
        "actual": actual,
    }
    if page_name:
        item["page"] = page_name
    if viewport:
        item["viewport"] = viewport
    checks.append(item)


def capture_page(page: Page, path: Path, *, full_page: bool) -> Dict[str, Any]:
    try:
        page.screenshot(
            path=str(path),
            full_page=full_page,
            animations="disabled",
            caret="hide",
            timeout=120_000,
        )
        result = png_metadata(path)
        result["status"] = "pass" if result.get("bytes", 0) > 0 else "fail"
        return result
    except Exception as exc:
        return {"path": str(path), "exists": False, "status": "fail", "error": str(exc)}


def capture_locator(locator: Locator, path: Path) -> Dict[str, Any]:
    try:
        locator.scroll_into_view_if_needed(timeout=30_000)
        locator.screenshot(
            path=str(path),
            animations="disabled",
            caret="hide",
            timeout=120_000,
        )
        result = png_metadata(path)
        result["status"] = "pass" if result.get("bytes", 0) > 0 else "fail"
        return result
    except Exception as exc:
        return {"path": str(path), "exists": False, "status": "fail", "error": str(exc)}


def register_screenshot(
    screenshots: List[Dict[str, Any]],
    key: Tuple[str, str, str],
    metadata: Dict[str, Any],
) -> Dict[str, Any]:
    planned = SCREENSHOT_PLAN[key]
    record = {
        "step": planned["step"],
        "description": planned["description"],
        "page": key[0],
        "viewport": key[1],
        "kind": key[2],
        **metadata,
    }
    screenshots.append(record)
    return record


def attach_runtime_handlers(page: Page, runtime: Dict[str, Any]) -> None:
    page.on(
        "console",
        lambda message: runtime["console_errors"].append(
            {"text": message.text, "location": message.location}
        )
        if message.type == "error"
        else None,
    )
    page.on("pageerror", lambda error: runtime["page_errors"].append(str(error)))
    page.on(
        "requestfailed",
        lambda request: runtime["request_failures"].append(
            {"url": request.url, "failure": request.failure}
        ),
    )
    page.on(
        "response",
        lambda response: runtime["response_errors"].append(
            {"url": response.url, "status": response.status}
        )
        if response.status >= 400
        else None,
    )


def wait_ready(page: Page, page_kind: str, timeout_ms: int) -> None:
    if page_kind == "index":
        page.wait_for_selector("#cohort-table tbody tr", state="visible", timeout=timeout_ms)
        page.wait_for_function(
            "document.querySelectorAll('#cohort-table tbody tr').length === 7",
            timeout=timeout_ms,
        )
    else:
        page.wait_for_function(
            "document.querySelectorAll('.chrom-button').length === 22 && "
            "document.querySelectorAll('#result-list .result-row').length > 0",
            timeout=timeout_ms,
        )
    page.wait_for_timeout(180)
    page.evaluate("scrollTo(0, 0)")


def visible_forbidden_hits(page: Page) -> Dict[str, Any]:
    text = page.locator("body").inner_text(timeout=30_000)
    lowered = text.lower()
    counts = {term: lowered.count(term.lower()) for term in FORBIDDEN_VISIBLE}
    return {
        "counts": counts,
        "total": sum(counts.values()),
        "claim_boundary_not_ccf_visible": "不是 ccf" in lowered,
    }


def collect_layout(page: Page) -> Dict[str, Any]:
    return page.evaluate(
        """() => {
            const all = selector => Array.from(document.querySelectorAll(selector));
            const visible = el => {
                const closed = el.closest('details:not([open])');
                if (closed && el !== closed && !el.matches('summary') && !el.closest('summary')) return false;
                const style = getComputedStyle(el);
                const rect = el.getBoundingClientRect();
                return style.display !== 'none' && style.visibility !== 'hidden' &&
                    Number(style.opacity || 1) !== 0 && rect.width > 0 && rect.height > 0;
            };
            const root = document.documentElement;
            const body = document.body;
            const scrollRegions = all('.table-scroll,.scroll-region,.network-scroll').filter(visible).map(el => {
                const style = getComputedStyle(el);
                return {
                    class: typeof el.className === 'string' ? el.className : '',
                    role: el.getAttribute('role'),
                    tabindex: el.getAttribute('tabindex'),
                    aria_label: el.getAttribute('aria-label'),
                    overflow_x: style.overflowX,
                    client_width: el.clientWidth,
                    scroll_width: el.scrollWidth,
                    scrollable: el.scrollWidth > el.clientWidth + 1
                };
            });
            const networks = all('#detail .network-scroll svg').filter(visible).map(svg => {
                const rect = svg.getBoundingClientRect();
                const intrinsic = Number(svg.getAttribute('width')) || svg.viewBox.baseVal.width || rect.width;
                const textNodes = Array.from(svg.querySelectorAll('text'));
                return {
                    intrinsic_width: intrinsic,
                    rendered_width: rect.width,
                    scale_ratio: intrinsic ? rect.width / intrinsic : 1,
                    text_nodes: textNodes.length,
                    font_sizes: textNodes.map(node => parseFloat(getComputedStyle(node).fontSize || '0')),
                    role: svg.getAttribute('role'),
                    aria_label: svg.getAttribute('aria-label'),
                    has_title: !!svg.querySelector('title')
                };
            });
            const rawLinks = all('a[href*=".json"]');
            const visibleRaw = rawLinks.filter(visible);
            const headings = all('h1,h2,h3,h4,h5,h6').filter(visible).map(el => ({
                level: Number(el.tagName.slice(1)),
                text: (el.innerText || '').replace(/\s+/g, ' ').trim()
            }));
            return {
                document: {
                    title: document.title,
                    lang: root.lang,
                    ready_state: document.readyState,
                    viewport: {width: innerWidth, height: innerHeight},
                    page: {width: root.scrollWidth, height: root.scrollHeight},
                    body_overflow_x_px: Math.max(0, Math.max(root.scrollWidth, body.scrollWidth) - innerWidth),
                    body_overflow_x_css: getComputedStyle(body).overflowX
                },
                chrom_buttons: all('.chrom-button').length,
                result_rows: all('#result-list .result-row').length,
                raw_json_links: rawLinks.length,
                raw_json_visible: visibleRaw.length,
                raw_json_outside_evidence_drawer: rawLinks.filter(link => !link.closest('.evidence-drawer')).map(link => link.href),
                evidence_drawers: all('.evidence-drawer').map(el => ({open: el.open, links: el.querySelectorAll('a[href*=".json"]').length})),
                scroll_regions: scrollRegions,
                networks,
                headings,
                broken_images: all('img').filter(img => !img.complete || !img.naturalWidth).map(img => img.src)
            };
        }"""
    )


def heading_issues(headings: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    issues: List[Dict[str, Any]] = []
    if sum(item["level"] == 1 for item in headings) != 1:
        issues.append({"type": "h1_count", "actual": sum(item["level"] == 1 for item in headings)})
    previous: Optional[int] = None
    for item in headings:
        level = int(item["level"])
        if previous is not None and level > previous + 1:
            issues.append({"type": "heading_skip", "from": previous, "to": level, "text": item["text"]})
        previous = level
    return issues


def validate_local_targets(page: Page, source_dir: Path) -> Dict[str, Any]:
    records = page.evaluate(
        """() => ({
            links: Array.from(document.querySelectorAll('a[href]')).map(el => ({
                text: (el.innerText || el.getAttribute('aria-label') || '').replace(/\s+/g, ' ').trim(),
                href: el.getAttribute('href') || '',
                resolved: el.href
            })),
            datasetOptions: Array.from(document.querySelectorAll('#dataset-switch option')).map(option => ({
                text: option.textContent.trim(), value: option.value
            }))
        })"""
    )
    results: List[Dict[str, Any]] = []
    for link in records["links"]:
        href = link["href"]
        item = dict(link)
        if href.startswith("#"):
            exists = page.evaluate(
                "href => { try { return href.length > 1 && !!document.querySelector(href); } catch (_) { return false; } }",
                href,
            )
            item.update({"kind": "anchor", "exists": bool(exists)})
        else:
            parsed = urlparse(link["resolved"])
            if parsed.scheme == "file":
                target = Path(unquote(parsed.path))
                item.update({"kind": "local_file", "target": str(target), "exists": target.is_file()})
            else:
                item.update({"kind": parsed.scheme or "relative", "exists": True})
        results.append(item)
    for option in records["datasetOptions"]:
        target = (source_dir / option["value"]).resolve()
        results.append(
            {
                "kind": "dataset_option",
                "text": option["text"],
                "target": str(target),
                "exists": target.is_file(),
            }
        )
    failed = [item for item in results if not item["exists"]]
    return {"tested": len(results), "failed": len(failed), "failures": failed, "results": results}


def choose_representative_region(page: Page) -> Dict[str, Any]:
    candidate = page.evaluate(
        """() => {
            const data = JSON.parse(document.getElementById('workstation-data').textContent);
            const chunks = {};
            const detail = row => {
                if (!chunks[row.chrom]) chunks[row.chrom] = JSON.parse(document.getElementById('chunk-' + row.chrom).textContent);
                return chunks[row.chrom][row.chunk_index];
            };
            const scoreRows = data.region_index.map((row, order) => {
                if (row.n_sSNV !== 8 || row.topology_class !== 'incomplete') return null;
                const region = detail(row);
                const primary = (region.lineages || []).filter(line => line.is_primary_lineage === true);
                const treeLines = primary.filter(line => (line.trees || []).some(tree => (tree.edges || []).length));
                const stored = treeLines.reduce((sum, line) => sum + Number(line.n_trees_stored || 0), 0);
                const shapes = treeLines.reduce((sum, line) => sum + Number(line.n_distinct_shapes_stored || line.n_distinct_shapes_exact || 0), 0);
                const maxEdges = treeLines.reduce((best, line) => Math.max(best, ...(line.trees || []).map(tree => (tree.edges || []).length), 0), 0);
                const score = (treeLines.length ? 1000000 : 0) + shapes * 10000 + Math.min(stored, 1000) * 10 + maxEdges;
                return {row, order, score, primary_tree_lanes: treeLines.length, stored_trees: stored, stored_shapes: shapes, max_edges: maxEdges};
            }).filter(Boolean).sort((a, b) => b.score - a.score || a.order - b.order);
            const picked = scoreRows[0];
            if (!picked) return null;
            return {
                region: picked.row.region,
                chrom: picked.row.chrom,
                n_sSNV: picked.row.n_sSNV,
                topology_class: picked.row.topology_class,
                evidence_mode: picked.row.evidence_mode,
                has_recurrence: picked.row.has_recurrence,
                primary_tree_lanes: picked.primary_tree_lanes,
                stored_trees: picked.stored_trees,
                stored_shapes: picked.stored_shapes,
                max_edges: picked.max_edges,
                preferred_incomplete_8_site: true
            };
        }"""
    )
    if not candidate:
        raise RuntimeError("No incomplete 8-site representative region was found")
    page.locator("#fq").fill(candidate["region"])
    row = page.locator(f'.result-row[data-region="{candidate["region"]}"]')
    row.wait_for(state="visible", timeout=30_000)
    started = time.monotonic()
    row.click(timeout=30_000)
    page.wait_for_function(
        "region => document.querySelector('#detail-title')?.textContent.trim() === region",
        arg=candidate["region"],
        timeout=30_000,
    )
    page.wait_for_selector("#detail .network-scroll svg", state="visible", timeout=30_000)
    page.wait_for_timeout(120)
    candidate["selection_milliseconds"] = round((time.monotonic() - started) * 1000, 1)
    candidate["rendered_networks"] = page.locator("#detail .network-scroll svg").count()
    candidate["focus_after_selection"] = page.evaluate(
        """() => {
            const active = document.activeElement;
            const detail = document.getElementById('detail');
            return {
                tag: active?.tagName.toLowerCase() || null,
                id: active?.id || null,
                class: typeof active?.className === 'string' ? active.className : null,
                detail_or_descendant: !!detail && (active === detail || detail.contains(active))
            };
        }"""
    )
    return candidate


def mark_candidate_network(page: Page) -> Dict[str, Any]:
    state = page.evaluate(
        """() => {
            const scrollers = Array.from(document.querySelectorAll('#detail .lane:not(.auxiliary) .candidate-view .network-scroll'));
            scrollers.forEach(scroller => scroller.removeAttribute('data-audit-network'));
            const ranked = scrollers.map((scroller, index) => {
                const svg = scroller.querySelector('svg');
                const score = Number(svg?.getAttribute('width') || 0) + Number(svg?.getAttribute('height') || 0) +
                    (svg?.querySelectorAll('line').length || 0) * 50;
                return {index, score, scroller, svg};
            }).sort((a, b) => b.score - a.score);
            if (!ranked.length) return null;
            const picked = ranked[0];
            picked.scroller.setAttribute('data-audit-network', 'true');
            return {
                exact_network_index: picked.index,
                intrinsic_width: Number(picked.svg?.getAttribute('width') || 0),
                intrinsic_height: Number(picked.svg?.getAttribute('height') || 0),
                edge_lines: picked.svg?.querySelectorAll('line').length || 0,
                text_nodes: picked.svg?.querySelectorAll('text').length || 0,
                scroller_client_width: picked.scroller?.clientWidth || 0,
                scroller_scroll_width: picked.scroller?.scrollWidth || 0
            };
        }"""
    )
    if not state:
        raise RuntimeError("No candidate network card could be marked")
    return state


def test_skip_link(page: Page, expected_href: str) -> Dict[str, Any]:
    page.evaluate(
        """() => {
            document.documentElement.dataset.auditScrollBehavior = document.documentElement.style.scrollBehavior || '';
            document.documentElement.style.scrollBehavior = 'auto';
            document.body.setAttribute('tabindex', '-1');
            document.body.focus();
            scrollTo(0, 0);
        }"""
    )
    page.keyboard.press("Tab")
    page.wait_for_timeout(80)
    focused = page.evaluate(
        """() => {
            const el = document.activeElement;
            if (!el) return null;
            const r = el.getBoundingClientRect();
            const s = getComputedStyle(el);
            return {
                tag: el.tagName.toLowerCase(),
                class: typeof el.className === 'string' ? el.className : '',
                href: el.getAttribute('href'),
                text: (el.innerText || el.getAttribute('aria-label') || '').trim(),
                in_view: r.top >= 0 && r.bottom <= innerHeight,
                focus_indicator: (s.outlineStyle !== 'none' && parseFloat(s.outlineWidth || 0) > 0) ||
                    (s.boxShadow && s.boxShadow !== 'none')
            };
        }"""
    )
    target: Dict[str, Any] = {}
    if focused and str(focused.get("href", "")).startswith("#"):
        page.keyboard.press("Enter")
        page.wait_for_timeout(100)
        target = page.evaluate(
            """href => {
                const el = document.querySelector(href);
                const r = el?.getBoundingClientRect();
                return {
                    exists: !!el,
                    in_view: !!r && r.top < innerHeight && r.bottom > 0,
                    hash: location.hash
                };
            }""",
            focused["href"],
        )
    passed = bool(
        focused
        and focused.get("href") == expected_href
        and focused.get("in_view")
        and focused.get("focus_indicator")
        and target.get("exists")
        and target.get("in_view")
    )
    page.evaluate(
        """() => {
            document.body.removeAttribute('tabindex');
            document.activeElement?.blur();
            history.replaceState(null, '', location.pathname + location.search);
            scrollTo(0, 0);
            document.documentElement.style.scrollBehavior = document.documentElement.dataset.auditScrollBehavior || '';
            delete document.documentElement.dataset.auditScrollBehavior;
        }"""
    )
    return {"status": "pass" if passed else "fail", "focused": focused, "target": target}


def test_keyboard_focus(page: Page, limit: int = 400) -> Dict[str, Any]:
    available = page.evaluate(
        """() => Array.from(document.querySelectorAll(
            'a[href],button,input,select,textarea,summary,[tabindex]:not([tabindex="-1"])'
        )).filter(el => {
            const closed = el.closest('details:not([open])');
            if (closed) {
                const ownSummary = closed.querySelector(':scope > summary');
                if (!ownSummary || (el !== ownSummary && !ownSummary.contains(el))) return false;
                if (closed.parentElement?.closest('details:not([open])')) return false;
            }
            const r = el.getBoundingClientRect();
            const s = getComputedStyle(el);
            return !el.disabled && el.tabIndex >= 0 && s.display !== 'none' && s.visibility !== 'hidden' && r.width > 0 && r.height > 0;
        }).length"""
    )
    page.evaluate(
        """() => {
            document.documentElement.style.scrollBehavior = 'auto';
            document.body.setAttribute('tabindex', '-1');
            document.body.focus();
            scrollTo(0, 0);
        }"""
    )
    sequence: List[Dict[str, Any]] = []
    for step in range(min(int(available) + 1, limit)):
        page.keyboard.press("Tab")
        page.wait_for_timeout(18)
        item = page.evaluate(
            """step => {
                const el = document.activeElement;
                if (!el) return {step, missing: true};
                const r = el.getBoundingClientRect();
                const s = getComputedStyle(el);
                return {
                    step,
                    tag: el.tagName.toLowerCase(),
                    id: el.id || null,
                    class: typeof el.className === 'string' ? el.className.slice(0, 100) : '',
                    role: el.getAttribute('role'),
                    name: (el.innerText || el.getAttribute('aria-label') || el.value || '').replace(/\s+/g, ' ').trim().slice(0, 180),
                    focus_indicator: (s.outlineStyle !== 'none' && parseFloat(s.outlineWidth || 0) > 0) ||
                        (s.boxShadow && s.boxShadow !== 'none'),
                    visible: r.width > 0 && r.height > 0 && s.display !== 'none' && s.visibility !== 'hidden'
                };
            }""",
            step + 1,
        )
        sequence.append(item)
        if step > 0 and item.get("tag") == "body":
            break
    page.evaluate(
        """() => {
            document.body.removeAttribute('tabindex');
            document.activeElement?.blur();
            document.documentElement.style.scrollBehavior = '';
            scrollTo(0, 0);
        }"""
    )
    real = [item for item in sequence if item.get("tag") != "body" and not item.get("missing")]
    missing = [item for item in real if item.get("visible") and not item.get("focus_indicator")]
    return {
        "available": available,
        "visited": len(real),
        "complete": len(real) >= available,
        "missing_indicator_count": len(missing),
        "missing": missing,
        "sequence": sequence,
    }


def test_details(page: Page) -> Dict[str, Any]:
    results: List[Dict[str, Any]] = []
    details = page.locator("details")
    for index in range(details.count()):
        detail = details.nth(index)
        original_states = page.evaluate(
            """index => {
                const all = Array.from(document.querySelectorAll('details'));
                const el = all[index];
                if (!el) return null;
                const ancestors = [];
                let parent = el.parentElement?.closest('details');
                while (parent) { ancestors.push({el: parent, open: parent.open}); parent = parent.parentElement?.closest('details'); }
                ancestors.reverse().forEach(item => item.el.open = true);
                return {open: el.open, ancestors: ancestors.map(item => item.open)};
            }""",
            index,
        )
        if original_states is None:
            continue
        try:
            summary = detail.locator("summary").first
            label = summary.inner_text(timeout=5_000).strip()[:200]
            before = detail.evaluate("el => el.open")
            summary.click(timeout=10_000)
            after = detail.evaluate("el => el.open")
            summary.click(timeout=10_000)
            restored = detail.evaluate("el => el.open") == before
            passed = after != before and restored
            results.append(
                {
                    "index": index,
                    "summary": label,
                    "toggled": after != before,
                    "restored": restored,
                    "status": "pass" if passed else "fail",
                }
            )
        except Exception as exc:
            results.append({"index": index, "status": "fail", "error": str(exc)})
        finally:
            page.evaluate(
                """payload => {
                    const all = Array.from(document.querySelectorAll('details'));
                    const el = all[payload.index];
                    if (!el) return;
                    el.open = payload.state.open;
                    const ancestors = [];
                    let parent = el.parentElement?.closest('details');
                    while (parent) { ancestors.push(parent); parent = parent.parentElement?.closest('details'); }
                    ancestors.reverse().forEach((node, i) => node.open = payload.state.ancestors[i]);
                }""",
                {"index": index, "state": original_states},
            )
    failed = [item for item in results if item["status"] != "pass"]
    return {"tested": len(results), "failed": len(failed), "results": results}


def reset_sample_filters(page: Page) -> None:
    page.evaluate(
        """() => {
            document.documentElement.style.scrollBehavior = 'auto';
            ['fchr','ftopo','fevidence','fsignal'].forEach(id => {
                const el = document.getElementById(id);
                if (el) { el.value = ''; el.dispatchEvent(new Event('change', {bubbles:true})); }
            });
            const query = document.getElementById('fq');
            if (query) { query.value = ''; query.dispatchEvent(new Event('input', {bubbles:true})); }
        }"""
    )
    page.wait_for_timeout(60)


def result_count(page: Page) -> int:
    text = page.locator("#result-count").inner_text().strip()
    match = re.search(r"[\d,]+", text)
    return int(match.group(0).replace(",", "")) if match else -1


def test_chromosome_buttons(page: Page) -> Dict[str, Any]:
    reset_sample_filters(page)
    expected = page.evaluate(
        """() => {
            const data = JSON.parse(document.getElementById('workstation-data').textContent);
            const counts = {};
            for (let i=1; i<=22; i++) counts['chr'+i] = data.region_index.filter(row => row.chrom === 'chr'+i).length;
            return counts;
        }"""
    )
    results: List[Dict[str, Any]] = []
    for number in range(1, 23):
        chrom = f"chr{number}"
        button = page.locator(f'.chrom-button[data-chrom="{chrom}"]')
        started = time.monotonic()
        button.click(timeout=20_000)
        page.wait_for_function(
            "chrom => document.getElementById('fchr').value === chrom",
            arg=chrom,
            timeout=20_000,
        )
        elapsed_ms = round((time.monotonic() - started) * 1000, 1)
        actual_count = result_count(page)
        rows = page.locator("#result-list .result-row").evaluate_all(
            "els => els.map(el => el.dataset.region)"
        )
        aria_state = button.evaluate(
            """el => ({
                pressed: el.getAttribute('aria-pressed'),
                current: el.getAttribute('aria-current'),
                selected: el.getAttribute('aria-selected'),
                active: el.classList.contains('active')
            })"""
        )
        semantic_active = (
            aria_state.get("pressed") == "true"
            or aria_state.get("current") in {"true", "page"}
            or aria_state.get("selected") == "true"
        )
        passed = (
            actual_count == int(expected[chrom])
            and 0 <= len(rows) <= 81
            and all(region.startswith(chrom + ":") for region in rows)
            and aria_state.get("active")
            and semantic_active
        )
        results.append(
            {
                "chrom": chrom,
                "expected_count": expected[chrom],
                "actual_count": actual_count,
                "rendered_rows": len(rows),
                "elapsed_ms": elapsed_ms,
                "aria_state": aria_state,
                "status": "pass" if passed else "fail",
            }
        )
    page.locator("#all-genome").click(timeout=20_000)
    page.wait_for_function("document.getElementById('fchr').value === ''", timeout=20_000)
    failed = [item for item in results if item["status"] != "pass"]
    return {
        "tested": len(results),
        "failed": len(failed),
        "first_interaction_ms": results[0]["elapsed_ms"] if results else None,
        "results": results,
    }


def expected_facet_count(page: Page, control_id: str, value: str) -> int:
    return int(
        page.evaluate(
            """payload => {
                const rows = JSON.parse(document.getElementById('workstation-data').textContent).region_index;
                if (payload.id === 'ftopo') return rows.filter(row => row.topology_class === payload.value).length;
                if (payload.id === 'fevidence') return rows.filter(row => row.evidence_mode === payload.value).length;
                if (payload.id === 'fsignal' && payload.value === 'recurrence') return rows.filter(row => row.has_recurrence).length;
                if (payload.id === 'fsignal' && payload.value === 'hidden') return rows.filter(row => row.hidden_positive).length;
                if (payload.id === 'fsignal' && payload.value === 'multi_hp') return rows.filter(row => row.hp_multiplicity === 2).length;
                return -1;
            }""",
            {"id": control_id, "value": value},
        )
    )


def test_facets(page: Page, sample_name: str) -> Dict[str, Any]:
    reset_sample_filters(page)
    results: List[Dict[str, Any]] = []
    for control_id in ("ftopo", "fevidence", "fsignal"):
        values = page.locator(f"#{control_id} option").evaluate_all(
            "options => options.map(option => option.value).filter(Boolean)"
        )
        for value in values:
            expected = expected_facet_count(page, control_id, value)
            page.locator(f"#{control_id}").select_option(value)
            page.wait_for_timeout(45)
            actual = result_count(page)
            regions = page.locator("#result-list .result-row").evaluate_all(
                "els => els.map(el => el.dataset.region)"
            )
            matches = page.evaluate(
                """payload => {
                    const rows = JSON.parse(document.getElementById('workstation-data').textContent).region_index;
                    const map = new Map(rows.map(row => [row.region, row]));
                    return payload.regions.every(region => {
                        const row = map.get(region);
                        if (!row) return false;
                        if (payload.id === 'ftopo') return row.topology_class === payload.value;
                        if (payload.id === 'fevidence') return row.evidence_mode === payload.value;
                        if (payload.value === 'recurrence') return row.has_recurrence;
                        if (payload.value === 'hidden') return row.hidden_positive;
                        if (payload.value === 'multi_hp') return row.hp_multiplicity === 2;
                        return false;
                    });
                }""",
                {"id": control_id, "value": value, "regions": regions},
            )
            passed = actual == expected and matches and len(regions) <= 81
            results.append(
                {
                    "control": control_id,
                    "value": value,
                    "expected_count": expected,
                    "actual_count": actual,
                    "rendered_rows": len(regions),
                    "rows_match": matches,
                    "status": "pass" if passed else "fail",
                }
            )
            page.locator(f"#{control_id}").select_option("")
            page.wait_for_timeout(35)
    combined: Optional[Dict[str, Any]] = None
    if sample_name == "H2009":
        page.locator("#ftopo").select_option("incomplete")
        page.locator("#fsignal").select_option("recurrence")
        page.wait_for_timeout(60)
        expected = int(
            page.evaluate(
                """() => JSON.parse(document.getElementById('workstation-data').textContent).region_index
                    .filter(row => row.topology_class === 'incomplete' && row.has_recurrence).length"""
            )
        )
        actual = result_count(page)
        regions = page.locator("#result-list .result-row").count()
        combined = {
            "filters": ["incomplete", "recurrence"],
            "expected_count": expected,
            "actual_count": actual,
            "rendered_rows": regions,
            "status": "pass" if expected > 0 and actual == expected and regions == min(80, expected) else "fail",
        }
        results.append({"control": "H2009 recurrence+incomplete", **combined})
    reset_sample_filters(page)
    failed = [item for item in results if item["status"] != "pass"]
    return {"tested": len(results), "failed": len(failed), "combined": combined, "results": results}


def find_tab_rich_region(page: Page) -> Optional[Dict[str, Any]]:
    return page.evaluate(
        """() => {
            const data = JSON.parse(document.getElementById('workstation-data').textContent);
            const chunks = {};
            const signature = edges => {
                const children = {}; const childNodes = new Set();
                (edges || []).forEach(edge => {(children[edge[0]] = children[edge[0]] || []).push(edge[1]); childNodes.add(edge[1]);});
                const roots = Object.keys(children).filter(node => !childNodes.has(node));
                const walk = node => '(' + ((children[node] || []).map(walk).sort().join('')) + ')';
                return roots.length ? roots.map(walk).sort().join('|') : '()';
            };
            for (const row of data.region_index) {
                if (row.topology_class === 'incomplete' || row.topology_class === 'no_primary_lineage') continue;
                if (!chunks[row.chrom]) chunks[row.chrom] = JSON.parse(document.getElementById('chunk-' + row.chrom).textContent);
                const region = chunks[row.chrom][row.chunk_index];
                for (const line of (region.lineages || [])) {
                    if (!line.is_primary_lineage || (line.trees || []).length < 3) continue;
                    const groups = {};
                    (line.trees || []).forEach((tree, index) => {(groups[signature(tree.edges || [])] = groups[signature(tree.edges || [])] || []).push(index);});
                    const sizes = Object.values(groups).map(group => group.length);
                    if (sizes.length > 1 && sizes.some(size => size > 1)) {
                        return {region: row.region, chrom: row.chrom, shapes: sizes.length, max_exact_in_shape: Math.max(...sizes)};
                    }
                }
            }
            return null;
        }"""
    )


def select_region_by_search(page: Page, region: str) -> None:
    reset_sample_filters(page)
    page.locator("#fq").fill(region)
    row = page.locator(f'.result-row[data-region="{region}"]')
    row.wait_for(state="visible", timeout=30_000)
    row.click(timeout=30_000)
    page.wait_for_function(
        "region => document.querySelector('#detail-title')?.textContent.trim() === region",
        arg=region,
        timeout=30_000,
    )
    page.wait_for_timeout(80)


def test_shape_exact_tabs(page: Page) -> Dict[str, Any]:
    target = find_tab_rich_region(page)
    if not target:
        return {"status": "fail", "reason": "No region with multiple shape and exact tabs"}
    select_region_by_search(page, target["region"])
    picked_lane: Optional[Locator] = None
    shape_texts: List[str] = []
    target_shape = None
    lanes = page.locator("#detail .lane:not(.auxiliary)")
    for lane_index in range(lanes.count()):
        lane = lanes.nth(lane_index)
        lane_tabs = lane.locator(".shape-tab")
        texts = lane_tabs.all_inner_texts()
        if len(texts) < 2:
            continue
        for index, text in enumerate(texts):
            match = re.search(r"·\s*([\d,]+)\s+(?:stored\s+)?exact", text)
            if match and int(match.group(1).replace(",", "")) > 1:
                picked_lane = lane
                shape_texts = texts
                target_shape = index
                break
        if picked_lane is not None:
            break
    if picked_lane is None or target_shape is None:
        return {"status": "fail", "target": target, "shape_tabs": shape_texts}
    shape_tabs = picked_lane.locator(".shape-tab")
    shape_tabs.nth(target_shape).click(timeout=20_000)
    page.wait_for_timeout(80)
    shape_selected = shape_tabs.nth(target_shape).evaluate(
        "el => el.getAttribute('aria-selected') === 'true' || el.getAttribute('aria-pressed') === 'true'"
    )
    tree_tabs = picked_lane.locator(".tree-tab")
    tree_count = tree_tabs.count()
    tree_selected = False
    before_label = picked_lane.locator(".candidate-view .network-head b").inner_text()
    after_label = before_label
    if tree_count > 1:
        tree_tabs.nth(1).click(timeout=20_000)
        page.wait_for_timeout(80)
        tree_selected = tree_tabs.nth(1).evaluate(
            "el => el.getAttribute('aria-selected') === 'true' || el.getAttribute('aria-pressed') === 'true'"
        )
        after_label = picked_lane.locator(".candidate-view .network-head b").inner_text()
    passed = shape_selected and tree_count > 1 and tree_selected and before_label != after_label
    return {
        "status": "pass" if passed else "fail",
        "target": target,
        "shape_tabs": shape_tabs.count(),
        "target_shape": target_shape,
        "shape_selected": shape_selected,
        "tree_tabs": tree_count,
        "tree_selected": tree_selected,
        "candidate_before": before_label,
        "candidate_after": after_label,
    }


def find_shape_scope_cases(page: Page) -> Dict[str, Any]:
    return page.evaluate(
        """() => {
            const data = JSON.parse(document.getElementById('workstation-data').textContent);
            const chunks = {};
            const out = {all_exact_shapes_stored: null, unshown_exact_shapes: null};
            for (const row of data.region_index) {
                if (out.all_exact_shapes_stored && out.unshown_exact_shapes) break;
                if (!chunks[row.chrom]) chunks[row.chrom] = JSON.parse(document.getElementById('chunk-' + row.chrom).textContent);
                const region = chunks[row.chrom][row.chunk_index];
                for (const line of (region.lineages || [])) {
                    if (!line.is_primary_lineage || line.analysis_candidate_set_complete !== true || line.display_trees_complete === true) continue;
                    const stored = Number(line.n_distinct_shapes_stored);
                    const exact = Number(line.n_distinct_shapes_exact);
                    if (!Number.isFinite(stored) || !Number.isFinite(exact)) continue;
                    const record = {region: row.region, fam_label: line.fam_label, stored_shapes: stored, exact_shapes: exact};
                    if (stored === exact && !out.all_exact_shapes_stored) out.all_exact_shapes_stored = record;
                    if (stored < exact && !out.unshown_exact_shapes) out.unshown_exact_shapes = record;
                }
            }
            return out;
        }"""
    )


def shape_scope_text_for_case(page: Page, case: Dict[str, Any]) -> str:
    select_region_by_search(page, case["region"])
    lanes = page.locator("#detail .lane:not(.auxiliary)")
    for index in range(lanes.count()):
        lane = lanes.nth(index)
        label = lane.locator(".lane-title b").inner_text().strip()
        if label == case["fam_label"]:
            return lane.locator(".shape-scope").inner_text().strip()
    raise RuntimeError(f"Cannot find lane {case['fam_label']} for {case['region']}")


def test_shape_scope_split(page: Page, sample_name: str) -> Dict[str, Any]:
    if sample_name != "H2009":
        return {"status": "not_applicable"}
    cases = find_shape_scope_cases(page)
    equal_case = cases.get("all_exact_shapes_stored")
    subset_case = cases.get("unshown_exact_shapes")
    if not equal_case or not subset_case:
        return {"status": "fail", "reason": "Both shape-scope cases are required", "cases": cases}
    equal_text = shape_scope_text_for_case(page, equal_case)
    subset_text = shape_scope_text_for_case(page, subset_case)
    equal_pass = "all" in equal_text and "exact topology groups" in equal_text and "unshown exact shapes" not in equal_text
    subset_pass = "unshown exact shapes remain" in subset_text
    return {
        "status": "pass" if equal_pass and subset_pass else "fail",
        "all_exact_shapes_stored": {**equal_case, "text": equal_text, "passed": equal_pass},
        "unshown_exact_shapes": {**subset_case, "text": subset_text, "passed": subset_pass},
    }


def test_copy_button(page: Page) -> Dict[str, Any]:
    button = page.locator("#copy-link")
    before = button.inner_text().strip()
    button.click(timeout=20_000)
    page.wait_for_timeout(100)
    after = button.inner_text().strip()
    passed = after in {"已複製檢視連結", "請複製網址列"}
    return {"status": "pass" if passed else "fail", "before": before, "after": after}


def load_more_state(page: Page) -> Dict[str, Any]:
    return page.locator("#load-more").evaluate(
        """el => {
            const style = getComputedStyle(el);
            const rect = el.getBoundingClientRect();
            return {
                hidden: el.hidden,
                display: style.display,
                text: el.textContent.trim(),
                visible: style.display !== 'none' && style.visibility !== 'hidden' && rect.width > 0 && rect.height > 0
            };
        }"""
    )


def test_load_more_contract(page: Page, representative_region: str) -> Dict[str, Any]:
    reset_sample_filters(page)
    initial = load_more_state(page)
    page.locator("#fq").fill(representative_region)
    page.wait_for_timeout(60)
    one_count = result_count(page)
    one = load_more_state(page)
    page.locator("#fq").fill("chr22:999999999-1000000000")
    page.wait_for_timeout(60)
    zero_count = result_count(page)
    zero = load_more_state(page)
    reset_sample_filters(page)
    nonnegative = all(not re.search(r"\(-\d", item["text"]) for item in (initial, one, zero))
    passed = (
        initial["hidden"] is False
        and initial["visible"]
        and one_count == 1
        and one["hidden"] is True
        and one["display"] == "none"
        and not one["visible"]
        and zero_count == 0
        and zero["hidden"] is True
        and zero["display"] == "none"
        and not zero["visible"]
        and nonnegative
    )
    return {
        "status": "pass" if passed else "fail",
        "initial": initial,
        "filtered_one": {"count": one_count, **one},
        "filtered_zero": {"count": zero_count, **zero},
        "nonnegative_text": nonnegative,
    }


def test_prev_next_back(page: Page) -> Dict[str, Any]:
    reset_sample_filters(page)
    rows = page.locator("#result-list .result-row")
    if rows.count() < 3:
        return {"status": "fail", "reason": "Fewer than three initial rows"}
    rows.nth(1).click(timeout=20_000)
    page.wait_for_selector("#next-region:not([disabled])", timeout=20_000)
    first = page.locator("#detail-title").inner_text().strip()
    page.locator("#next-region").click(timeout=20_000)
    page.wait_for_function(
        "before => document.querySelector('#detail-title')?.textContent.trim() !== before",
        arg=first,
        timeout=20_000,
    )
    second = page.locator("#detail-title").inner_text().strip()
    page.locator("#prev-region").click(timeout=20_000)
    page.wait_for_function(
        "wanted => document.querySelector('#detail-title')?.textContent.trim() === wanted",
        arg=first,
        timeout=20_000,
    )
    page.evaluate("document.documentElement.style.scrollBehavior='auto'")
    page.locator("#back-results").click(timeout=20_000)
    page.wait_for_timeout(100)
    back_in_view = page.locator("#results-panel").evaluate(
        "el => { const r=el.getBoundingClientRect(); return r.top < innerHeight && r.bottom > 0; }"
    )
    passed = first != second and page.locator("#detail-title").inner_text().strip() == first and back_in_view
    return {
        "status": "pass" if passed else "fail",
        "first": first,
        "next": second,
        "prev_restored": page.locator("#detail-title").inner_text().strip(),
        "back_results_in_view": back_in_view,
    }


def test_browser_history_rehydrate(page: Page, source: Path, timeout_ms: int) -> Dict[str, Any]:
    page.goto("about:blank", wait_until="load", timeout=timeout_ms)
    page.goto(source.as_uri(), wait_until="load", timeout=timeout_ms)
    wait_ready(page, "sample", timeout_ms)
    reset_sample_filters(page)
    page.locator("#ftopo").select_option("incomplete")
    page.wait_for_timeout(60)
    rows = page.locator("#result-list .result-row")
    if rows.count() < 1:
        return {"status": "fail", "reason": "No incomplete rows"}
    first = rows.nth(0).get_attribute("data-region")
    rows.nth(0).click(timeout=20_000)
    page.wait_for_function(
        "region => document.querySelector('#detail-title')?.textContent.trim() === region",
        arg=first,
        timeout=timeout_ms,
    )
    # Change the filter without creating a separate history entry; selecting
    # the second region must snapshot this distinct filter state.
    page.evaluate(
        """() => {
            document.getElementById('ftopo').value = 'exact_and_topology_unique';
            applyFilters();
        }"""
    )
    page.wait_for_timeout(60)
    second_row = page.locator("#result-list .result-row").first
    if second_row.count() < 1:
        return {"status": "fail", "reason": "No exact-and-topology-unique rows"}
    second = second_row.get_attribute("data-region")
    second_row.click(timeout=20_000)
    page.wait_for_function(
        "region => document.querySelector('#detail-title')?.textContent.trim() === region",
        arg=second,
        timeout=timeout_ms,
    )
    page.evaluate("history.back()")
    page.wait_for_function(
        "region => document.querySelector('#detail-title')?.textContent.trim() === region && "
        "new URLSearchParams(location.hash.slice(1)).get('region') === region",
        arg=first,
        timeout=timeout_ms,
    )
    back_state = page.evaluate(
        """() => ({
            detail: document.querySelector('#detail-title')?.textContent.trim(),
            hash: new URLSearchParams(location.hash.slice(1)).get('region'),
            chrom: document.getElementById('fchr')?.value,
            topology: document.getElementById('ftopo')?.value,
            detail_focused: document.activeElement === document.getElementById('detail') ||
                document.getElementById('detail')?.contains(document.activeElement)
        })"""
    )
    page.evaluate("history.forward()")
    page.wait_for_function(
        "region => document.querySelector('#detail-title')?.textContent.trim() === region && "
        "new URLSearchParams(location.hash.slice(1)).get('region') === region",
        arg=second,
        timeout=timeout_ms,
    )
    forward_state = page.evaluate(
        """() => ({
            detail: document.querySelector('#detail-title')?.textContent.trim(),
            hash: new URLSearchParams(location.hash.slice(1)).get('region'),
            chrom: document.getElementById('fchr')?.value,
            topology: document.getElementById('ftopo')?.value,
            detail_focused: document.activeElement === document.getElementById('detail') ||
                document.getElementById('detail')?.contains(document.activeElement)
        })"""
    )
    passed = (
        back_state["detail"] == first
        and back_state["hash"] == first
        and back_state["chrom"] == str(first).split(":", 1)[0]
        and back_state["topology"] == "incomplete"
        and back_state["detail_focused"]
        and forward_state["detail"] == second
        and forward_state["hash"] == second
        and forward_state["chrom"] == str(second).split(":", 1)[0]
        and forward_state["topology"] == "exact_and_topology_unique"
        and forward_state["detail_focused"]
    )
    return {
        "status": "pass" if passed else "fail",
        "first": first,
        "second": second,
        "back": back_state,
        "forward": forward_state,
    }


def test_deep_link(page: Page, source: Path, region: str, timeout_ms: int) -> Dict[str, Any]:
    from urllib.parse import quote

    url = source.as_uri() + "#region=" + quote(region, safe="")
    # A same-document fragment navigation does not rerun init(); force a fresh
    # document so this verifies the real cold deep-link contract.
    page.goto("about:blank", wait_until="load", timeout=timeout_ms)
    page.goto(url, wait_until="load", timeout=timeout_ms)
    wait_ready(page, "sample", timeout_ms)
    page.wait_for_function(
        "region => document.querySelector('#detail-title')?.textContent.trim() === region",
        arg=region,
        timeout=timeout_ms,
    )
    state = page.evaluate(
        """region => {
            const chrom = region.split(':')[0];
            return {
                detail: document.querySelector('#detail-title')?.textContent.trim(),
                filter_chrom: document.getElementById('fchr')?.value,
                active_chrom: document.querySelector('.chrom-button.active')?.dataset.chrom,
                hash_region: new URLSearchParams(location.hash.slice(1)).get('region'),
                current_rows: document.querySelectorAll('.result-row[aria-current="true"]').length
            };
        }""",
        region,
    )
    chrom = region.split(":", 1)[0]
    passed = (
        state.get("detail") == region
        and state.get("filter_chrom") == chrom
        and state.get("active_chrom") == chrom
        and state.get("hash_region") == region
    )
    return {"status": "pass" if passed else "fail", "url": url, "state": state}


def test_no_primary_caption(page: Page, source: Path, timeout_ms: int) -> Dict[str, Any]:
    if source.name != "HCC1395_DORADO.html":
        return {"status": "not_applicable"}
    region = "chr1:190064024-190196077"
    deep = test_deep_link(page, source, region, timeout_ms)
    caption = page.locator("#detail .site-table caption").inner_text().strip()
    detail_text = page.locator("#detail").inner_text()
    passed = (
        deep["status"] == "pass"
        and "auxiliary／control units 加總" in caption
        and "primary HP lineages 加總" not in caption
        and "不進 C／Topo" in caption
        and any(marker in detail_text for marker in ("無 mutation-bearing HP1／HP2", "沒有 mutation-bearing HP1／HP2"))
    )
    return {
        "status": "pass" if passed else "fail",
        "region": region,
        "caption": caption,
        "deep_link": deep,
    }


def test_sample_meta(page: Page) -> Dict[str, Any]:
    return page.evaluate(
        """() => {
            const data = JSON.parse(document.getElementById('workstation-data').textContent);
            const actual = {
                current_summary: document.querySelector('meta[name="intersubmod-current-summary-sha256"]')?.content || null,
                region_view: document.querySelector('meta[name="intersubmod-region-view-sha256"]')?.content || null,
                canonical_sample: document.querySelector('meta[name="intersubmod-canonical-sample"]')?.content || null
            };
            const expected = {
                current_summary: data.source.machine_summary_sha256,
                region_view: data.source.region_view_sha256,
                canonical_sample: data.sample
            };
            return {
                actual,
                expected,
                exact_matches: Object.keys(expected).filter(key => actual[key] === expected[key]).length,
                status: Object.keys(expected).every(key => actual[key] === expected[key]) ? 'pass' : 'fail'
            };
        }"""
    )


def evaluate_scroller_contract(layout: Dict[str, Any]) -> Dict[str, Any]:
    regions = layout["scroll_regions"]
    failures = []
    for item in regions:
        if (
            item.get("role") != "region"
            or item.get("tabindex") != "0"
            or not str(item.get("aria_label") or "").strip()
            or item.get("overflow_x") not in {"auto", "scroll"}
        ):
            failures.append(item)
    networks = layout["networks"]
    network_failures = []
    for item in networks:
        fonts = item.get("font_sizes", [])
        min_font = min(fonts) if fonts else 999
        if (
            abs(float(item.get("scale_ratio", 0)) - 1.0) > 0.01
            or min_font < 10
            or item.get("role") != "img"
            or not item.get("aria_label")
            or not item.get("has_title")
        ):
            network_failures.append({**item, "min_font": min_font})
    network_region_names = [
        str(item.get("aria_label") or "")
        for item in regions
        if "network-scroll" in str(item.get("class") or "")
    ]
    network_svg_names = [str(item.get("aria_label") or "") for item in networks]
    duplicate_network_region_names = sorted(
        {name for name in network_region_names if network_region_names.count(name) > 1}
    )
    duplicate_network_svg_names = sorted(
        {name for name in network_svg_names if network_svg_names.count(name) > 1}
    )
    return {
        "regions_tested": len(regions),
        "region_failures": failures,
        "networks_tested": len(networks),
        "network_failures": network_failures,
        "duplicate_network_region_names": duplicate_network_region_names,
        "duplicate_network_svg_names": duplicate_network_svg_names,
        "status": "pass"
        if not failures
        and not network_failures
        and not duplicate_network_region_names
        and not duplicate_network_svg_names
        else "fail",
    }


def test_local_scroller_keyboard(page: Page) -> Dict[str, Any]:
    scrollers = page.locator("#detail .scroll-region, #detail .network-scroll")
    results: List[Dict[str, Any]] = []
    for index in range(scrollers.count()):
        scroller = scrollers.nth(index)
        state = scroller.evaluate(
            """el => {
                const box=el.getBoundingClientRect();const style=getComputedStyle(el);
                return {
                    class: typeof el.className === 'string' ? el.className : '',
                    label: el.getAttribute('aria-label'),
                    visible: !el.closest('details:not([open])') && box.width>0 && box.height>0 &&
                        style.display!=='none' && style.visibility!=='hidden',
                    scrollable: el.scrollWidth > el.clientWidth + 1,
                    before: el.scrollLeft,
                    client_width: el.clientWidth,
                    scroll_width: el.scrollWidth
                };
            }"""
        )
        if not state["visible"]:
            continue
        if not state["scrollable"]:
            continue
        scroller.focus()
        page.keyboard.press("ArrowRight")
        page.wait_for_timeout(70)
        after = scroller.evaluate("el => el.scrollLeft")
        scroller.evaluate("el => el.scrollLeft=0")
        results.append(
            {
                **state,
                "after": after,
                "status": "pass" if after > state["before"] else "fail",
            }
        )
    failed = [item for item in results if item["status"] != "pass"]
    return {
        "tested": len(results),
        "failed": len(failed),
        "status": "pass" if not failed else "fail",
        "not_applicable_no_overflow": not results,
        "results": results,
    }


def test_candidate_group_semantics(page: Page) -> Dict[str, Any]:
    groups = page.evaluate(
        """() => Array.from(document.querySelectorAll('#detail .shape-tabs,#detail .tree-tabs')).map(group => {
            const buttons = Array.from(group.querySelectorAll('button'));
            const role = group.getAttribute('role');
            const label = group.getAttribute('aria-label');
            const buttonGroup = role === 'group' && buttons.every(button =>
                button.getAttribute('role') !== 'tab' && ['true','false'].includes(button.getAttribute('aria-pressed'))
            );
            const tabPattern = role === 'tablist' && buttons.every(button => {
                const controlled = button.getAttribute('aria-controls');
                const panel = controlled ? document.getElementById(controlled) : null;
                return button.getAttribute('role') === 'tab' &&
                    ['true','false'].includes(button.getAttribute('aria-selected')) &&
                    !!panel && panel.getAttribute('role') === 'tabpanel';
            });
            return {
                class: group.className,
                role,
                label,
                buttons: buttons.length,
                button_group_valid: buttonGroup,
                tab_pattern_valid: tabPattern,
                valid: !!label && buttons.length > 0 && (buttonGroup || tabPattern)
            };
        })"""
    )
    failures = [item for item in groups if not item["valid"]]
    return {
        "tested": len(groups),
        "failed": len(failures),
        "status": "pass" if groups and not failures else "fail",
        "groups": groups,
    }


def test_acquisition_label_layout(page: Page) -> Dict[str, Any]:
    networks = page.evaluate(
        """() => Array.from(document.querySelectorAll('#detail .candidate-view .network-scroll svg')).map((svg, svgIndex) => {
            const texts = Array.from(svg.querySelectorAll('text')).map((node, index) => {
                const value = Array.from(node.childNodes)
                    .filter(child => child.nodeType === Node.TEXT_NODE)
                    .map(child => child.nodeValue || '')
                    .join('').trim();
                const box = node.getBoundingClientRect();
                return {
                    index,
                    value,
                    acquisition: /^\+\s*S\d+(?:\s*,\s*S\d+)*$/.test(value),
                    starts_plus: value.startsWith('+'),
                    long_suffix: /repeated|site|acquisition/i.test(value),
                    box: {left:box.left, right:box.right, top:box.top, bottom:box.bottom, width:box.width, height:box.height}
                };
            });
            const acquisition = texts.filter(item => item.starts_plus);
            const overlaps = [];
            acquisition.forEach(label => texts.forEach(other => {
                if (label.index === other.index) return;
                const x = Math.min(label.box.right, other.box.right) - Math.max(label.box.left, other.box.left);
                const y = Math.min(label.box.bottom, other.box.bottom) - Math.max(label.box.top, other.box.top);
                if (x > 1 && y > 1) overlaps.push({label:label.value, other:other.value, overlap_x:x, overlap_y:y});
            }));
            return {
                svg_index: svgIndex,
                labels: acquisition.map(item => item.value),
                invalid_short_labels: acquisition.filter(item => !item.acquisition).map(item => item.value),
                long_suffix_labels: acquisition.filter(item => item.long_suffix).map(item => item.value),
                overlaps
            };
        })"""
    )
    labels = [label for network in networks for label in network["labels"]]
    invalid = [label for network in networks for label in network["invalid_short_labels"]]
    suffixes = [label for network in networks for label in network["long_suffix_labels"]]
    overlaps = [item for network in networks for item in network["overlaps"]]
    passed = bool(networks and labels and not invalid and not suffixes and not overlaps)
    return {
        "status": "pass" if passed else "fail",
        "networks": len(networks),
        "labels": labels,
        "invalid_short_labels": invalid,
        "long_suffix_labels": suffixes,
        "overlaps": overlaps,
        "details": networks,
    }


def test_index(page: Page) -> Dict[str, Any]:
    text = page.locator("main").inner_text()
    digits = re.sub(r"[^0-9]", "", text)
    normalized = re.sub(r"[\s,]", "", text)
    extracted = page.evaluate(
        """() => ({
            metrics: Array.from(document.querySelectorAll('.metric-card')).map(card => ({
                code: card.querySelector('.metric-code')?.textContent.trim() || '',
                number: card.querySelector('.metric-number')?.textContent.trim() || ''
            })),
            aggregateTopologyTitles: Array.from(document.querySelector('.topology-main .topology-bar')?.querySelectorAll('.bar-segment') || [])
                .map(segment => segment.getAttribute('title'))
        })"""
    )
    aggregate = {
        "W_tree": bool(re.search(r"W_tree[^\d]{0,80}51,815", text, re.I | re.S)),
        "W_primary": "W_primary50215" in normalized,
        "complete": bool(re.search(r"Complete[^\d]{0,80}42,240", text, re.I | re.S)),
        "incomplete": bool(re.search(r"Incomplete[^\d]{0,80}7,975", text, re.I | re.S)),
        "topology_exact": "11582" in digits,
        "topology_shape": "10737" in digits,
        "topology_multiple": "19921" in digits,
    }
    scope = "chr1–22" in text and ("全基因" in text or "Whole-genome" in text)
    launchers = page.locator("nav.genome-launchers a.genome-link").count()
    table = page.locator("#cohort-table table")
    table_semantics = {
        "captions": table.locator("caption").count(),
        "headers": table.locator("th").count(),
        "scoped_headers": table.locator("th[scope]").count(),
    }
    table_semantics["status"] = (
        "pass"
        if table_semantics["captions"] == 1
        and table_semantics["headers"] == table_semantics["scoped_headers"]
        else "fail"
    )
    scroller = page.locator("#cohort-table .table-scroll")
    horizontal = scroller.evaluate(
        "el => ({scrollable:el.scrollWidth>el.clientWidth+1,before:el.scrollLeft,client:el.clientWidth,width:el.scrollWidth})"
    )
    scroller.focus()
    page.keyboard.press("ArrowRight")
    page.wait_for_timeout(80)
    horizontal["after"] = scroller.evaluate("el => el.scrollLeft")
    horizontal["status"] = (
        "pass"
        if (not horizontal["scrollable"] or horizontal["after"] > horizontal["before"])
        else "fail"
    )
    scroller.evaluate("el => el.scrollLeft=0")
    return {
        "aggregate": aggregate,
        "extracted_values": extracted,
        "aggregate_status": "pass" if all(aggregate.values()) else "fail",
        "scope_visible": scope,
        "launchers": launchers,
        "table_semantics": table_semantics,
        "horizontal_keyboard_scroll": horizontal,
    }


def test_sample_controls(
    page: Page,
    source: Path,
    representative: Dict[str, Any],
    timeout_ms: int,
) -> Dict[str, Any]:
    sample_name = source.stem
    def stage(name: str, operation: Any) -> Any:
        try:
            return operation()
        except Exception as exc:
            raise RuntimeError(f"sample control stage {name}: {exc}") from exc

    chrom = stage("chromosome_buttons", lambda: test_chromosome_buttons(page))
    facets = stage("facets", lambda: test_facets(page, sample_name))
    tabs = stage("shape_exact_tabs", lambda: test_shape_exact_tabs(page))
    shape_scope = stage("shape_scope_split", lambda: test_shape_scope_split(page, sample_name))
    load_more = stage("load_more", lambda: test_load_more_contract(page, representative["region"]))
    copy = stage("copy", lambda: test_copy_button(page))
    navigation = stage("prev_next_back", lambda: test_prev_next_back(page))
    browser_history = stage("browser_history", lambda: test_browser_history_rehydrate(page, source, timeout_ms))
    deep_link = stage("deep_link", lambda: test_deep_link(page, source, representative["region"], timeout_ms))
    no_primary = stage("no_primary_caption", lambda: test_no_primary_caption(page, source, timeout_ms))
    stage("restore_representative", lambda: select_region_by_search(page, representative["region"]))
    details = stage("details", lambda: test_details(page))
    tests = {
        "chromosome_buttons": chrom,
        "facets": facets,
        "shape_exact_tabs": tabs,
        "shape_scope_split": shape_scope,
        "load_more": load_more,
        "copy": copy,
        "prev_next_back": navigation,
        "browser_history": browser_history,
        "deep_link": deep_link,
        "no_primary_caption": no_primary,
        "details": details,
    }
    failed = []
    if chrom["failed"]:
        failed.append("chromosome_buttons")
    if facets["failed"]:
        failed.append("facets")
    for name in ("shape_exact_tabs", "load_more", "copy", "prev_next_back", "browser_history", "deep_link"):
        if tests[name].get("status") != "pass":
            failed.append(name)
    if shape_scope.get("status") not in {"pass", "not_applicable"}:
        failed.append("shape_scope_split")
    if no_primary.get("status") not in {"pass", "not_applicable"}:
        failed.append("no_primary_caption")
    if details["failed"]:
        failed.append("details")
    return {"status": "pass" if not failed else "fail", "failed": failed, "tests": tests}


def add_common_checks(
    checks: List[Dict[str, Any]],
    result: Dict[str, Any],
    *,
    page_name: str,
    viewport_name: str,
) -> None:
    for screenshot in result["screenshots"]:
        add_check(
            checks,
            f"screenshot_step_{screenshot['step']:02d}",
            screenshot.get("status") == "pass" and screenshot.get("bytes", 0) > 0,
            "non-empty PNG",
            screenshot,
            page_name=page_name,
            viewport=viewport_name,
        )
    for error_kind, errors in result["runtime"].items():
        add_check(
            checks,
            error_kind,
            not errors,
            {"count": 0},
            {"count": len(errors), "items": errors},
            page_name=page_name,
            viewport=viewport_name,
        )
    initial = result["initial"]
    document = initial["layout"]["document"]
    add_check(
        checks,
        "body_horizontal_overflow",
        document["body_overflow_x_px"] == 0,
        0,
        document["body_overflow_x_px"],
        page_name=page_name,
        viewport=viewport_name,
    )
    add_check(
        checks,
        "body_overflow_css",
        document["body_overflow_x_css"] in {"hidden", "clip"}
        or document["body_overflow_x_px"] == 0,
        "hidden/clip, or visible only when measured horizontal overflow is 0",
        {
            "overflow_x_css": document["body_overflow_x_css"],
            "overflow_x_px": document["body_overflow_x_px"],
        },
        page_name=page_name,
        viewport=viewport_name,
    )
    add_check(
        checks,
        "heading_hierarchy",
        not initial["heading_issues"],
        "one H1 and no skipped visible levels",
        initial["heading_issues"],
        page_name=page_name,
        viewport=viewport_name,
    )
    add_check(
        checks,
        "local_href_targets",
        initial["local_targets"]["failed"] == 0,
        {"failed": 0},
        {"tested": initial["local_targets"]["tested"], "failures": initial["local_targets"]["failures"]},
        page_name=page_name,
        viewport=viewport_name,
    )
    raw = initial["layout"]
    expected_raw_links = 9 if page_name == INDEX_FILE else 4
    add_check(
        checks,
        "raw_json_links_collapsed",
        raw["raw_json_links"] == expected_raw_links
        and raw["raw_json_visible"] == 0
        and not raw["raw_json_outside_evidence_drawer"]
        and all(not drawer["open"] for drawer in raw["evidence_drawers"]),
        f"{expected_raw_links} raw JSON links, all inside closed .evidence-drawer; initially visible 0",
        {
            "links": raw["raw_json_links"],
            "visible": raw["raw_json_visible"],
            "outside": raw["raw_json_outside_evidence_drawer"],
            "drawers": raw["evidence_drawers"],
        },
        page_name=page_name,
        viewport=viewport_name,
    )
    add_check(
        checks,
        "visible_forbidden_text_initial",
        initial["forbidden"]["total"] == 0,
        {"total": 0},
        initial["forbidden"],
        page_name=page_name,
        viewport=viewport_name,
    )
    add_check(
        checks,
        "broken_images",
        not raw["broken_images"],
        {"count": 0},
        raw["broken_images"],
        page_name=page_name,
        viewport=viewport_name,
    )
    add_check(
        checks,
        "offline_mode",
        result["offline"],
        True,
        result["offline"],
        page_name=page_name,
        viewport=viewport_name,
    )


def audit_page_viewport(
    browser: Browser,
    source: Path,
    source_dir: Path,
    output_dir: Path,
    viewport_name: str,
    timeout_ms: int,
    screenshots: List[Dict[str, Any]],
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    page_name = source.name
    page_kind = "index" if page_name == INDEX_FILE else "sample"
    viewport = VIEWPORTS[viewport_name]
    checks: List[Dict[str, Any]] = []
    result: Dict[str, Any] = {
        "page": page_name,
        "page_kind": page_kind,
        "viewport": viewport_name,
        "configured_viewport": viewport,
        "offline": True,
        "capture_order": [
            "01_load_local_file_with_context_offline",
            "02_capture_initial_full_page_before_dom_interpretation",
            "03_collect_initial_dom_and_link_metrics",
            "04_select_representative_region_for_sample",
            "05_capture_region_and_candidate_network",
            "06_test_interactions_and_accessibility",
        ],
        "runtime": {"console_errors": [], "page_errors": [], "request_failures": [], "response_errors": []},
        "screenshots": [],
    }
    context = browser.new_context(
        viewport=viewport,
        color_scheme="light",
        reduced_motion="reduce",
        locale="zh-TW",
        device_scale_factor=1,
    )
    context.set_offline(True)
    page = context.new_page()
    page.set_default_timeout(timeout_ms)
    attach_runtime_handlers(page, result["runtime"])
    started = time.monotonic()
    try:
        nav_started = time.monotonic()
        page.goto(source.as_uri(), wait_until="load", timeout=timeout_ms)
        wait_ready(page, page_kind, timeout_ms)
        result["load"] = {
            "status": "pass",
            "milliseconds": round((time.monotonic() - nav_started) * 1000, 1),
            "ready_state": page.evaluate("document.readyState"),
            "url": page.url,
        }

        full_key = (page_name, viewport_name, "full_page")
        full_path = output_dir / SCREENSHOT_DIR_NAME / SCREENSHOT_PLAN[full_key]["filename"]
        full_record = register_screenshot(
            screenshots,
            full_key,
            capture_page(page, full_path, full_page=True),
        )
        result["screenshots"].append(full_record)

        initial_layout = collect_layout(page)
        result["initial"] = {
            "layout": initial_layout,
            "heading_issues": heading_issues(initial_layout["headings"]),
            "local_targets": validate_local_targets(page, source_dir),
            "forbidden": visible_forbidden_hits(page),
        }
        result["skip_link"] = (
            test_skip_link(page, "#main-content")
            if page_kind == "index"
            else test_skip_link(page, "#genome-overview")
        )

        if page_kind == "index":
            result["index_tests"] = test_index(page)
            result["initial_scroller_contract"] = evaluate_scroller_contract(initial_layout)
            result["keyboard_focus"] = test_keyboard_focus(page)
        else:
            result["meta"] = test_sample_meta(page)
            result["scope_visible"] = page.evaluate(
                """() => {
                    const text = document.body.innerText;
                    return text.includes('chr1–22') && (text.includes('Whole-genome scope') || text.includes('全基因視圖'));
                }"""
            )
            result["cohort_count_visible"] = page.evaluate(
                """() => {
                    const text = document.querySelector('.hero')?.innerText || '';
                    return /7\s+datasets?/i.test(text) && /6\s+biological\s+samples?/i.test(text);
                }"""
            )
            if viewport_name in INTERACTION_VIEWPORTS:
                representative = choose_representative_region(page)
                result["representative"] = representative
                page.locator("#detail").scroll_into_view_if_needed(timeout=30_000)
                page.evaluate("scrollBy(0,-8)")
                page.wait_for_timeout(100)
                region_key = (page_name, viewport_name, "region_view")
                region_path = output_dir / SCREENSHOT_DIR_NAME / SCREENSHOT_PLAN[region_key]["filename"]
                region_record = register_screenshot(
                    screenshots,
                    region_key,
                    capture_page(page, region_path, full_page=False),
                )
                result["screenshots"].append(region_record)

                result["network_capture"] = mark_candidate_network(page)
                network_key = (page_name, viewport_name, "candidate_network")
                network_path = output_dir / SCREENSHOT_DIR_NAME / SCREENSHOT_PLAN[network_key]["filename"]
                network_locator = page.locator('[data-audit-network="true"]')
                network_locator.scroll_into_view_if_needed(timeout=30_000)
                page.evaluate("scrollBy(0,-110)")
                page.wait_for_timeout(80)
                network_record = register_screenshot(
                    screenshots,
                    network_key,
                    capture_locator(network_locator, network_path),
                )
                result["screenshots"].append(network_record)

                selected_layout = collect_layout(page)
                result["selected"] = {
                    "layout": selected_layout,
                    "scroller_contract": evaluate_scroller_contract(selected_layout),
                    "keyboard_scroll": test_local_scroller_keyboard(page),
                    "candidate_group_semantics": test_candidate_group_semantics(page),
                    "acquisition_label_layout": test_acquisition_label_layout(page),
                    "forbidden": visible_forbidden_hits(page),
                }

                if page_name == "HCC1395_DORADO.html" and viewport_name == "desktop":
                    result["no_primary_precheck"] = test_no_primary_caption(page, source, timeout_ms)
                    page.locator("#detail").scroll_into_view_if_needed(timeout=30_000)
                    page.evaluate("scrollBy(0,-8)")
                    no_primary_key = (page_name, viewport_name, "no_primary_region")
                    no_primary_path = output_dir / SCREENSHOT_DIR_NAME / SCREENSHOT_PLAN[no_primary_key]["filename"]
                    no_primary_record = register_screenshot(
                        screenshots,
                        no_primary_key,
                        capture_page(page, no_primary_path, full_page=False),
                    )
                    result["screenshots"].append(no_primary_record)
                    result["no_primary_forbidden"] = visible_forbidden_hits(page)
                    select_region_by_search(page, representative["region"])

                if viewport_name == "desktop":
                    result["controls"] = test_sample_controls(page, source, representative, timeout_ms)
                    select_region_by_search(page, representative["region"])
                else:
                    result["mobile_controls"] = {
                        "copy": test_copy_button(page),
                        "details": test_details(page),
                    }
                    select_region_by_search(page, representative["region"])
                result["keyboard_focus"] = test_keyboard_focus(page)

        add_check(
            checks,
            "document_load",
            True,
            "local file loaded and rendered while browser context was offline",
            result["load"],
            page_name=page_name,
            viewport=viewport_name,
        )
        add_common_checks(checks, result, page_name=page_name, viewport_name=viewport_name)
        add_check(
            checks,
            "skip_link_keyboard",
            result["skip_link"]["status"] == "pass",
            "first Tab exposes skip link and Enter reaches its target",
            result["skip_link"],
            page_name=page_name,
            viewport=viewport_name,
        )

        if page_kind == "index":
            index_tests = result["index_tests"]
            add_check(checks, "index_exact_aggregate", index_tests["aggregate_status"] == "pass", EXPECTED_AGGREGATE, index_tests["aggregate"], page_name=page_name, viewport=viewport_name)
            add_check(checks, "index_whole_genome_scope", index_tests["scope_visible"], True, index_tests["scope_visible"], page_name=page_name, viewport=viewport_name)
            add_check(checks, "index_dataset_launchers", index_tests["launchers"] == 7, 7, index_tests["launchers"], page_name=page_name, viewport=viewport_name)
            add_check(checks, "index_table_semantics", index_tests["table_semantics"]["status"] == "pass", "caption=1 and all th have scope", index_tests["table_semantics"], page_name=page_name, viewport=viewport_name)
            add_check(checks, "index_table_keyboard_scroll", index_tests["horizontal_keyboard_scroll"]["status"] == "pass", "ArrowRight scrolls when table is wider than its region", index_tests["horizontal_keyboard_scroll"], page_name=page_name, viewport=viewport_name)
            contract = result["initial_scroller_contract"]
            add_check(checks, "index_scroll_region_contract", contract["status"] == "pass", "labeled role=region tabindex=0 local horizontal containment", contract, page_name=page_name, viewport=viewport_name)
        else:
            initial = result["initial"]["layout"]
            add_check(checks, "sample_chromosome_buttons", initial["chrom_buttons"] == 22, 22, initial["chrom_buttons"], page_name=page_name, viewport=viewport_name)
            add_check(checks, "sample_initial_dom_rows", initial["result_rows"] <= 80, {"maximum": 80}, initial["result_rows"], page_name=page_name, viewport=viewport_name)
            add_check(checks, "sample_three_exact_meta", result["meta"]["status"] == "pass" and result["meta"]["exact_matches"] == 3, 3, result["meta"], page_name=page_name, viewport=viewport_name)
            add_check(checks, "sample_whole_genome_scope", result["scope_visible"], True, result["scope_visible"], page_name=page_name, viewport=viewport_name)
            add_check(checks, "sample_cohort_count_visible", result["cohort_count_visible"], "7 datasets / 6 biological samples visible", result["cohort_count_visible"], page_name=page_name, viewport=viewport_name)
            if viewport_name in INTERACTION_VIEWPORTS:
                representative = result["representative"]
                add_check(checks, "representative_incomplete_8_site", representative["preferred_incomplete_8_site"] and representative["n_sSNV"] == 8 and representative["topology_class"] == "incomplete" and representative["rendered_networks"] > 0, "incomplete 8-site region with rendered candidate network", representative, page_name=page_name, viewport=viewport_name)
                add_check(checks, "initial_interaction_time_measured", isinstance(representative["selection_milliseconds"], (int, float)) and representative["selection_milliseconds"] >= 0, "numeric milliseconds", representative["selection_milliseconds"], page_name=page_name, viewport=viewport_name)
                if viewport_name == "mobile":
                    add_check(checks, "mobile_selection_moves_focus_to_detail", representative["focus_after_selection"]["detail_or_descendant"], "active focus is #detail or a descendant after selecting a mobile result", representative["focus_after_selection"], page_name=page_name, viewport=viewport_name)
                add_check(checks, "selected_scroller_and_network_contract", result["selected"]["scroller_contract"]["status"] == "pass", "wide tables/networks locally contained; role/name/tabindex; network text >=10px and scale=1", result["selected"]["scroller_contract"], page_name=page_name, viewport=viewport_name)
                add_check(checks, "selected_local_keyboard_scroll", result["selected"]["keyboard_scroll"]["status"] == "pass", "each visible wide table/network scrolls with ArrowRight", result["selected"]["keyboard_scroll"], page_name=page_name, viewport=viewport_name)
                add_check(checks, "candidate_selector_semantics", result["selected"]["candidate_group_semantics"]["status"] == "pass", "shape/exact selectors use a complete button-group or tab pattern", result["selected"]["candidate_group_semantics"], page_name=page_name, viewport=viewport_name)
                add_check(checks, "candidate_acquisition_label_layout", result["selected"]["acquisition_label_layout"]["status"] == "pass", "visual acquisition labels are short +Si tokens with no long suffix or text-box collision", result["selected"]["acquisition_label_layout"], page_name=page_name, viewport=viewport_name)
                if viewport_name == "mobile":
                    scrolled_classes = [item["class"] for item in result["selected"]["keyboard_scroll"]["results"]]
                    mobile_wide_pass = any("network-scroll" in value for value in scrolled_classes) and any("scroll-region" in value for value in scrolled_classes)
                    add_check(checks, "mobile_wide_table_and_network_scrollable", mobile_wide_pass, "at least one site table and one candidate network are locally scrollable", scrolled_classes, page_name=page_name, viewport=viewport_name)
                add_check(checks, "visible_forbidden_text_selected", result["selected"]["forbidden"]["total"] == 0, {"total": 0}, result["selected"]["forbidden"], page_name=page_name, viewport=viewport_name)
                if viewport_name == "desktop":
                    add_check(checks, "sample_full_control_suite", result["controls"]["status"] == "pass", {"failed": []}, result["controls"], page_name=page_name, viewport=viewport_name)
                    chrom_ms = result["controls"]["tests"]["chromosome_buttons"]["first_interaction_ms"]
                    add_check(checks, "chrom_filter_interaction_time_measured", isinstance(chrom_ms, (int, float)) and chrom_ms >= 0, "numeric milliseconds", chrom_ms, page_name=page_name, viewport=viewport_name)
                else:
                    mobile = result["mobile_controls"]
                    mobile_pass = mobile["copy"]["status"] == "pass" and mobile["details"]["failed"] == 0
                    add_check(checks, "mobile_copy_and_details", mobile_pass, "copy works and all details toggle", mobile, page_name=page_name, viewport=viewport_name)
                if page_name == "HCC1395_DORADO.html" and viewport_name == "desktop":
                    add_check(checks, "no_primary_auxiliary_caption", result["no_primary_precheck"]["status"] == "pass", "auxiliary/control caption and no primary wording", result["no_primary_precheck"], page_name=page_name, viewport=viewport_name)
                    add_check(checks, "no_primary_forbidden_text", result["no_primary_forbidden"]["total"] == 0, {"total": 0}, result["no_primary_forbidden"], page_name=page_name, viewport=viewport_name)

        if "keyboard_focus" in result:
            focus = result["keyboard_focus"]
            add_check(checks, "keyboard_focus_complete", focus["complete"] and focus["missing_indicator_count"] == 0, "all visible focusables visited with visible indicator", {"available": focus["available"], "visited": focus["visited"], "missing": focus["missing_indicator_count"], "targets": focus["missing"]}, page_name=page_name, viewport=viewport_name)
    except Exception as exc:
        result["fatal_error"] = str(exc)
        add_check(checks, "page_viewport_audit", False, "completed without fatal exception", str(exc), page_name=page_name, viewport=viewport_name)
    finally:
        result["duration_seconds"] = round(time.monotonic() - started, 3)
        context.close()
    return result, checks


def load_before_comparison(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        return {"available": False, "path": str(path)}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return {"available": False, "path": str(path), "error": str(exc)}
    summary = payload.get("summary", {})
    failed = summary.get("failed_checks", [])
    categories: Dict[str, int] = {}
    for item in failed:
        name = item.get("name", "unknown")
        categories[name] = categories.get(name, 0) + 1
    return {
        "available": True,
        "path": str(path),
        "audit_kind": payload.get("audit_kind"),
        "summary": summary,
        "failure_categories": categories,
        "source_hashes_unchanged": payload.get("input", {}).get("unchanged"),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the complete 8-page layered workstation Chromium after audit.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Exit codes: 0 all hard gates pass, 1 findings remain, 2 runner/config error.",
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--metrics", type=Path, default=DEFAULT_METRICS)
    parser.add_argument("--before-metrics", type=Path, default=DEFAULT_BEFORE_METRICS)
    parser.add_argument("--timeout-ms", type=int, default=90_000)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    source_dir = args.input_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    metrics_path = args.metrics.expanduser().resolve()
    before_metrics = args.before_metrics.expanduser().resolve()

    if PLAYWRIGHT_IMPORT_ERROR:
        print(json.dumps({"status": "error", "error": PLAYWRIGHT_IMPORT_ERROR}, ensure_ascii=False))
        return 2
    if args.timeout_ms <= 0:
        print(json.dumps({"status": "error", "error": "--timeout-ms must be positive"}, ensure_ascii=False))
        return 2
    if output_dir == source_dir or source_dir in output_dir.parents:
        print(json.dumps({"status": "error", "error": "Output cannot be inside the source HTML directory"}, ensure_ascii=False))
        return 2
    if metrics_path.parent != output_dir:
        print(json.dumps({"status": "error", "error": f"Metrics must be directly under {output_dir}"}, ensure_ascii=False))
        return 2

    sources = [source_dir / filename for filename in ALL_FILES]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        print(json.dumps({"status": "error", "missing": missing}, ensure_ascii=False, indent=2))
        return 2

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / SCREENSHOT_DIR_NAME).mkdir(parents=True, exist_ok=True)
    files_before = {path.name: file_record(path) for path in sources}
    report: Dict[str, Any] = {
        "schema_version": "2.0",
        "audit_kind": "layered_workstation_complete_after_visual_audit",
        "task_type": "B_comprehensive_validation",
        "partial": False,
        "generated_at": utc_now(),
        "scope": {
            "pages_expected": 8,
            "pages": ALL_FILES,
            "sample_pages_expected": 7,
            "viewports": VIEWPORTS,
            "page_viewports_expected": 24,
            "screenshots_expected": 53,
            "desktop_mobile_region_screenshots_expected": 14,
            "desktop_mobile_network_screenshots_expected": 14,
            "narrow_320_smokes_expected": 8,
            "partial": False,
        },
        "input": {"directory": str(source_dir), "files_before": files_before},
        "output": {
            "directory": str(output_dir),
            "metrics": str(metrics_path),
            "screenshots": str(output_dir / SCREENSHOT_DIR_NAME),
        },
        "before_comparison": load_before_comparison(before_metrics),
        "browser": None,
        "page_viewports": [],
        "screenshots": [],
        "checks": [],
        "fatal_errors": [],
    }
    playwright = None
    browser = None
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
            "offline_context": True,
        }
        for source in sources:
            for viewport_name in VIEWPORTS:
                page_result, page_checks = audit_page_viewport(
                    browser,
                    source,
                    source_dir,
                    output_dir,
                    viewport_name,
                    args.timeout_ms,
                    report["screenshots"],
                )
                report["page_viewports"].append(page_result)
                report["checks"].extend(page_checks)
    except Exception as exc:
        report["fatal_errors"].append(str(exc))
    finally:
        if browser is not None:
            try:
                browser.close()
            except Exception as exc:
                report["fatal_errors"].append(f"browser.close: {exc}")
        if playwright is not None:
            try:
                playwright.stop()
            except Exception as exc:
                report["fatal_errors"].append(f"playwright.stop: {exc}")

    files_after = {path.name: file_record(path) for path in sources}
    hash_changes = []
    for name in ALL_FILES:
        if files_before[name]["sha256"] != files_after[name]["sha256"]:
            hash_changes.append(
                {"file": name, "before": files_before[name]["sha256"], "after": files_after[name]["sha256"]}
            )
    report["input"]["files_after"] = files_after
    report["input"]["hash_changes"] = hash_changes
    report["input"]["unchanged"] = not hash_changes
    add_check(report["checks"], "source_html_hashes_unchanged", not hash_changes, {"changes": 0}, hash_changes)

    report["screenshots"].sort(key=lambda item: item["step"])
    screenshot_steps = [item["step"] for item in report["screenshots"]]
    add_check(
        report["checks"],
        "complete_screenshot_plan",
        screenshot_steps == list(range(1, 54))
        and all(item.get("status") == "pass" for item in report["screenshots"]),
        {"steps": list(range(1, 54)), "failed": 0},
        {
            "steps": screenshot_steps,
            "failed": [item["step"] for item in report["screenshots"] if item.get("status") != "pass"],
        },
    )
    audited_pairs = {
        (item.get("page"), item.get("viewport"))
        for item in report["page_viewports"]
        if not item.get("fatal_error")
    }
    expected_pairs = {(name, viewport) for name in ALL_FILES for viewport in VIEWPORTS}
    add_check(
        report["checks"],
        "complete_page_viewport_scope",
        audited_pairs == expected_pairs,
        sorted(expected_pairs),
        sorted(audited_pairs),
    )
    narrow_records = [
        item for item in report["page_viewports"] if item.get("viewport") == "narrow" and not item.get("fatal_error")
    ]
    add_check(report["checks"], "complete_320px_smoke_scope", len(narrow_records) == 8, 8, len(narrow_records))
    sample_desktop_controls = [
        item for item in report["page_viewports"]
        if item.get("page_kind") == "sample" and item.get("viewport") == "desktop" and item.get("controls")
    ]
    add_check(
        report["checks"],
        "seven_sample_full_interaction_suites",
        len(sample_desktop_controls) == 7 and all(item["controls"]["status"] == "pass" for item in sample_desktop_controls),
        {"samples": 7, "failed": 0},
        {
            "samples": len(sample_desktop_controls),
            "failed": [item["page"] for item in sample_desktop_controls if item["controls"]["status"] != "pass"],
        },
    )

    failed = [item for item in report["checks"] if item["status"] != "pass"]
    if report["fatal_errors"]:
        status = "error"
        exit_code = 2
    elif failed:
        status = "fail"
        exit_code = 1
    else:
        status = "pass"
        exit_code = 0
    before_summary = report["before_comparison"].get("summary", {})
    report["comparison"] = {
        "before": {
            "checks_total": before_summary.get("checks_total"),
            "checks_failed": before_summary.get("checks_failed"),
            "screenshots": before_summary.get("screenshots"),
            "page_viewports": before_summary.get("page_viewports_audited"),
            "failure_categories": report["before_comparison"].get("failure_categories", {}),
        },
        "after": {
            "checks_total": len(report["checks"]),
            "checks_failed": len(failed),
            "screenshots": len(report["screenshots"]),
            "page_viewports": len(report["page_viewports"]),
        },
        "failed_check_delta": (
            len(failed) - int(before_summary.get("checks_failed", 0))
            if before_summary.get("checks_failed") is not None
            else None
        ),
    }
    report["summary"] = {
        "status": status,
        "exit_code": exit_code,
        "pages_audited": len({item.get("page") for item in report["page_viewports"]}),
        "page_viewports_audited": len(report["page_viewports"]),
        "screenshots": len(report["screenshots"]),
        "checks_total": len(report["checks"]),
        "checks_passed": len(report["checks"]) - len(failed),
        "checks_failed": len(failed),
        "failed_checks": [
            {
                "name": item["name"],
                "page": item.get("page"),
                "viewport": item.get("viewport"),
                "severity": item.get("severity"),
                "actual": item["actual"],
            }
            for item in failed
        ],
    }
    report["finished_at"] = utc_now()
    metrics_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": status,
                "exit_code": exit_code,
                "input": str(source_dir / "*.html"),
                "metrics": str(metrics_path),
                "screenshots": len(report["screenshots"]),
                "pages": report["summary"]["pages_audited"],
                "page_viewports": report["summary"]["page_viewports_audited"],
                "checks_total": report["summary"]["checks_total"],
                "checks_failed": report["summary"]["checks_failed"],
                "source_hashes_unchanged": report["input"]["unchanged"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
