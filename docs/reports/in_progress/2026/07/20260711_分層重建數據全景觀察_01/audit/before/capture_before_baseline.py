#!/usr/bin/env python3
"""Capture reproducible before/after QA for the standalone panorama report.

The script opens the source with a ``file://`` URL in headless Chromium.  For
each viewport it captures the initial viewport and full page before collecting
DOM metrics or exercising safe, reversible interactions.  It never writes to
the source HTML; all outputs stay next to this script.

Exit codes:
    0: capture completed and all baseline checks passed
    1: capture completed but one or more report checks failed
    2: capture runner/configuration failed
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

try:
    from playwright.sync_api import Browser, Locator, Page, sync_playwright
except ImportError as exc:
    Browser = Locator = Page = Any  # type: ignore[misc,assignment]
    sync_playwright = None  # type: ignore[assignment]
    PLAYWRIGHT_IMPORT_ERROR: Optional[str] = str(exc)
else:
    PLAYWRIGHT_IMPORT_ERROR = None


SCRIPT_PATH = Path(__file__).resolve()
OUTPUT_DIR = SCRIPT_PATH.parent
AUDIT_ROOT = SCRIPT_PATH.parent.parent
REPORT_DIR = SCRIPT_PATH.parents[2]
DEFAULT_INPUT = REPORT_DIR / "20260711_分層重建全景數據觀察_01.standalone.html"
DEFAULT_METRICS = OUTPUT_DIR / "metrics.json"
VIEWPORTS = {
    "desktop": {"width": 1440, "height": 1000},
    "mobile": {"width": 390, "height": 844},
}
SECTION_LIMIT = 14
FOCUS_LIMIT = 100
DETAILS_LIMIT = 40
NAV_LIMIT = 12


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def slugify(text: str, fallback: str) -> str:
    text = re.sub(r"\s+", "_", text.strip())
    text = re.sub(r"[^0-9A-Za-z\u4e00-\u9fff_.-]+", "", text)
    return text[:64] or fallback


def png_metadata(path: Path) -> Dict[str, Any]:
    result: Dict[str, Any] = {"path": str(path), "exists": path.is_file()}
    if not path.is_file():
        return result
    result["bytes"] = path.stat().st_size
    try:
        with path.open("rb") as handle:
            header = handle.read(24)
        if header[:8] == b"\x89PNG\r\n\x1a\n" and len(header) >= 24:
            result["width"], result["height"] = struct.unpack(">II", header[16:24])
    except OSError as exc:
        result["metadata_error"] = str(exc)
    return result


def add_check(
    checks: List[Dict[str, Any]],
    name: str,
    passed: bool,
    expected: Any,
    actual: Any,
    *,
    viewport: Optional[str] = None,
    severity: str = "error",
) -> None:
    check = {
        "name": name,
        "status": "pass" if passed else "fail",
        "severity": severity,
        "expected": expected,
        "actual": actual,
    }
    if viewport:
        check["viewport"] = viewport
    checks.append(check)


def capture_png(page: Page, path: Path, *, full_page: bool) -> Dict[str, Any]:
    try:
        page.screenshot(
            path=str(path),
            full_page=full_page,
            animations="disabled",
            caret="hide",
        )
        result = png_metadata(path)
        result["status"] = "pass"
        result["full_page"] = full_page
        return result
    except Exception as exc:
        return {
            "path": str(path),
            "exists": False,
            "status": "fail",
            "full_page": full_page,
            "error": str(exc),
        }


def heading_issues(headings: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    issues: List[Dict[str, Any]] = []
    h1_count = sum(item["level"] == 1 for item in headings)
    if h1_count != 1:
        issues.append({"type": "h1_count", "expected": 1, "actual": h1_count})
    previous: Optional[int] = None
    for item in headings:
        level = int(item["level"])
        if previous is not None and level > previous + 1:
            issues.append(
                {
                    "type": "heading_level_skip",
                    "from": previous,
                    "to": level,
                    "text": item["text"],
                    "index": item["index"],
                }
            )
        previous = level
    return issues


def collect_dom_metrics(page: Page) -> Dict[str, Any]:
    return page.evaluate(
        """() => {
            const visible = el => {
                const closed = el.closest('details:not([open])');
                if (closed && el !== closed && !el.matches('summary') && !el.closest('summary')) return false;
                const r = el.getBoundingClientRect();
                const s = getComputedStyle(el);
                return s.display !== 'none' && s.visibility !== 'hidden' &&
                    Number(s.opacity || 1) !== 0 && r.width > 0 && r.height > 0;
            };
            const box = el => {
                const r = el.getBoundingClientRect();
                return {
                    x: Math.round(r.x * 10) / 10,
                    y: Math.round((r.y + scrollY) * 10) / 10,
                    width: Math.round(r.width * 10) / 10,
                    height: Math.round(r.height * 10) / 10
                };
            };
            const short = (value, n=180) => String(value || '').replace(/\s+/g, ' ').trim().slice(0, n);
            const all = selector => Array.from(document.querySelectorAll(selector));
            const root = document.documentElement;
            const body = document.body;
            const viewportWidth = innerWidth;
            const rootOverflow = Math.max(root.scrollWidth, body ? body.scrollWidth : 0) - viewportWidth;
            const overflowOffenders = all('body *').filter(el => {
                if (!visible(el)) return false;
                const r = el.getBoundingClientRect();
                return r.left < -1 || r.right > viewportWidth + 1;
            }).slice(0, 24).map(el => ({
                tag: el.tagName.toLowerCase(),
                id: el.id || null,
                class: typeof el.className === 'string' ? short(el.className, 120) : null,
                box: box(el),
                    local_scroller: !!el.closest('[style*="overflow"],.table-wrap,.table-scroll,.chart-wrap,.svg-explainer,[class*="scroll"]')
            }));
            const headings = all('h1,h2,h3,h4,h5,h6').filter(visible).map((el, index) => ({
                index,
                level: Number(el.tagName.slice(1)),
                tag: el.tagName.toLowerCase(),
                id: el.id || null,
                text: short(el.innerText, 500),
                box: box(el)
            }));
            const tables = all('table').filter(visible).map((el, index) => {
                const rect = el.getBoundingClientRect();
                const headers = Array.from(el.querySelectorAll('thead th, tr:first-child th'))
                    .map(th => short(th.innerText, 120));
                const rows = Array.from(el.rows || []);
                return {
                    index,
                    id: el.id || null,
                    class: typeof el.className === 'string' ? short(el.className, 120) : null,
                    caption: short(el.caption ? el.caption.innerText : ''),
                    headers,
                    rows: rows.length,
                    body_rows: el.tBodies ? Array.from(el.tBodies).reduce((n, tb) => n + tb.rows.length, 0) : null,
                    columns_max: rows.reduce((n, row) => Math.max(n, row.cells.length), 0),
                    box: box(el),
                    exceeds_viewport_width: rect.width > innerWidth + 1,
                    clipped_by_parent: !!el.parentElement && el.scrollWidth > el.parentElement.clientWidth + 1,
                    parent_overflow_x: el.parentElement ? getComputedStyle(el.parentElement).overflowX : null
                };
            });
            const rawCharts = all('svg,canvas,img');
            const seenCharts = new Set();
            const charts = rawCharts.filter(el => {
                if (!visible(el) || seenCharts.has(el)) return false;
                seenCharts.add(el);
                return true;
            }).map((el, index) => {
                const rect = el.getBoundingClientRect();
                return {
                    index,
                    tag: el.tagName.toLowerCase(),
                    id: el.id || null,
                    class: typeof el.className === 'string' ? short(el.className, 120) : null,
                    label: short(el.getAttribute('aria-label') || el.getAttribute('alt') ||
                        (el.getAttribute('aria-labelledby') || '').split(/\s+/).map(id =>
                            document.getElementById(id)?.innerText || '').join(' ') ||
                        el.querySelector?.('figcaption')?.innerText || el.querySelector?.('title')?.textContent),
                    box: box(el),
                    too_small_for_data_display: rect.width < 240 || rect.height < 120,
                    intrinsic_width: el.naturalWidth || el.width?.baseVal?.value || el.width || null,
                    intrinsic_height: el.naturalHeight || el.height?.baseVal?.value || el.height || null
                };
            });
            const chartHosts = all('[data-recharts-chart]').map((host, index) => {
                const figure = host.closest('figure');
                const labelOwner = figure || host;
                const labelledBy = (labelOwner.getAttribute('aria-labelledby') || '').split(/\s+/).filter(Boolean);
                const svg = host.querySelector('[data-recharts-live] svg');
                return {
                    index,
                    key: host.getAttribute('data-recharts-chart'),
                    title: short(labelledBy.map(id => document.getElementById(id)?.innerText || '').join(' '), 500),
                    live_svg: !!svg,
                    box: svg ? box(svg) : box(host),
                    horizontally_scrollable: !!host.closest('.chart-wrap,[class*="scroll"]')
                };
            });
            const internalLinks = all('a[href^="#"]').filter(visible).map((el, index) => {
                const href = el.getAttribute('href') || '';
                let target = null;
                if (href.length > 1) {
                    try { target = document.querySelector(href); } catch (_) { target = null; }
                }
                return {
                    index,
                    text: short(el.innerText || el.getAttribute('aria-label'), 240),
                    href,
                    target_exists: !!target,
                    box: box(el)
                };
            });
            const details = all('details').filter(visible).map((el, index) => ({
                index,
                id: el.id || null,
                open: el.open,
                summary: short(el.querySelector('summary')?.innerText, 300),
                box: box(el)
            }));
            const focusSelector = 'a[href],button,input,select,textarea,summary,[tabindex]:not([tabindex="-1"])';
            const focusables = all(focusSelector).filter(visible).map((el, index) => ({
                index,
                tag: el.tagName.toLowerCase(),
                id: el.id || null,
                role: el.getAttribute('role'),
                tabindex: el.getAttribute('tabindex'),
                name: short(el.innerText || el.getAttribute('aria-label') || el.getAttribute('title') || el.value, 240),
                box: box(el)
            }));
            const brokenImages = all('img').filter(img => !img.complete || img.naturalWidth === 0).map(img => ({
                src: img.getAttribute('src'),
                alt: img.getAttribute('alt')
            }));
            const skipLinks = all('a[href^="#"]').filter(el => /跳到|跳至|skip/i.test(
                el.innerText || el.getAttribute('aria-label') || '')).map(el => ({
                    text: short(el.innerText || el.getAttribute('aria-label')),
                    href: el.getAttribute('href')
                }));
            const companionDataLinks = all('.data-drawer a[href$=".json"],a.chart-source[href$=".json"]')
                .map(el => ({
                    href: el.getAttribute('href'),
                    text: short(el.innerText || el.getAttribute('aria-label'), 160),
                    in_collapsed_drawer: !!el.closest('.data-drawer:not([open])'),
                    chart_source: el.classList.contains('chart-source')
                }));
            const visibleMainText = document.querySelector('main')?.innerText || '';
            return {
                document: {
                    title: document.title,
                    lang: root.lang || null,
                    ready_state: document.readyState,
                    viewport: {width: innerWidth, height: innerHeight, device_pixel_ratio: devicePixelRatio},
                    page: {width: root.scrollWidth, height: root.scrollHeight},
                    body: body ? box(body) : null,
                    body_overflow_x_px: Math.max(0, Math.round(rootOverflow * 10) / 10),
                    overflow_offenders: overflowOffenders
                },
                counts: {
                    elements: all('*').length,
                    h1: all('h1').filter(visible).length,
                    h2: all('h2').filter(visible).length,
                    h3: all('h3').filter(visible).length,
                    headings: headings.length,
                    main: all('main').filter(visible).length,
                    nav: all('nav,[role="navigation"]').filter(visible).length,
                    sections: all('section,[role="region"]').filter(visible).length,
                    tables: tables.length,
                    table_rows: tables.reduce((n, table) => n + table.rows, 0),
                    charts: charts.length,
                    chart_hosts: all('[data-recharts-chart]').length,
                    live_chart_svg: all('[data-recharts-live] svg').length,
                    role_img: all('[role="img"]').filter(visible).length,
                    svg: all('svg').filter(visible).length,
                    canvas: all('canvas').filter(visible).length,
                    images: all('img').filter(visible).length,
                    links: all('a[href]').filter(visible).length,
                    internal_links: internalLinks.length,
                    buttons: all('button').filter(visible).length,
                    inputs: all('input,select,textarea').filter(visible).length,
                    details: details.length,
                    details_open: details.filter(item => item.open).length,
                    focusables: focusables.length,
                    positive_tabindex: focusables.filter(item => Number(item.tabindex) > 0).length
                },
                headings,
                tables,
                charts,
                chart_hosts: chartHosts,
                navigation: {
                    internal_links: internalLinks,
                    missing_internal_targets: internalLinks.filter(item => item.href.length > 1 && !item.target_exists),
                    skip_links: skipLinks
                },
                source_disclosure: {
                    inline_source_tooltips: all('.source-tooltip.inline-source').length,
                    sourced_value_annotations: all('.sourced-value[data-source-key]').length,
                    companion_data_links: companionDataLinks,
                    unique_companion_hrefs: Array.from(new Set(companionDataLinks.map(item => item.href))),
                    visible_json_filename_mentions: (visibleMainText.match(/\S+\.json\b/g) || [])
                },
                details,
                focusables,
                broken_images: brokenImages
            };
        }"""
    )


def discover_sections(page: Page) -> List[Dict[str, Any]]:
    return page.evaluate(
        """limit => {
            const visible = el => {
                const closed = el.closest('details:not([open])');
                if (closed && el !== closed && !el.matches('summary') && !el.closest('summary')) return false;
                const r = el.getBoundingClientRect();
                const s = getComputedStyle(el);
                return s.display !== 'none' && s.visibility !== 'hidden' && r.width > 0 && r.height >= 100;
            };
            const short = value => String(value || '').replace(/\s+/g, ' ').trim();
            const hasData = el => !!el.querySelector('table,svg,canvas,figure,img,.kpi,.metric,[class*="chart" i],[class*="plot" i]');
            const candidates = [];
            const seen = new Set();
            const add = (el, heading, source) => {
                if (!el || seen.has(el) || !visible(el)) return;
                const rect = el.getBoundingClientRect();
                const text = short(el.innerText);
                if (!hasData(el) && text.length < 350) return;
                if (rect.height > document.documentElement.scrollHeight * 0.85) return;
                seen.add(el);
                candidates.push({el, heading: short(heading), source, rect});
            };
            document.querySelectorAll('main section,article section,body > section,[role="region"],.report-section,.data-section').forEach(el => {
                const h = el.querySelector('h2,h3,h4');
                add(el, h?.innerText || el.getAttribute('aria-label') || el.id, 'semantic-section');
            });
            document.querySelectorAll('h2,h3').forEach(h => {
                let el = h.closest('section,[role="region"],article,.panel,.card,.chapter');
                if (!el) {
                    el = h.parentElement;
                    while (el && el !== document.body && !hasData(el) && el.innerText.length < 350) {
                        el = el.parentElement;
                    }
                }
                add(el, h.innerText, 'heading-container');
            });
            candidates.sort((a, b) => a.rect.top - b.rect.top || a.rect.height - b.rect.height);
            return candidates.slice(0, limit).map((item, index) => {
                item.el.setAttribute('data-audit-capture-section', String(index));
                const rect = item.el.getBoundingClientRect();
                return {
                    index,
                    heading: item.heading || `section-${index + 1}`,
                    source: item.source,
                    tag: item.el.tagName.toLowerCase(),
                    id: item.el.id || null,
                    class: typeof item.el.className === 'string' ? item.el.className.slice(0, 180) : null,
                    box: {
                        x: Math.round(rect.x * 10) / 10,
                        y: Math.round((rect.y + scrollY) * 10) / 10,
                        width: Math.round(rect.width * 10) / 10,
                        height: Math.round(rect.height * 10) / 10
                    }
                };
            });
        }""",
        SECTION_LIMIT,
    )


def capture_sections(
    page: Page,
    viewport_name: str,
    *,
    reuse_screenshots: bool,
) -> List[Dict[str, Any]]:
    sections = discover_sections(page)
    results: List[Dict[str, Any]] = []
    for section in sections:
        name = slugify(section["heading"], f"section_{section['index'] + 1:02d}")
        path = OUTPUT_DIR / f"{viewport_name}_section_{section['index'] + 1:02d}_{name}.png"
        locator = page.locator(
            f'[data-audit-capture-section="{section["index"]}"]'
        ).first
        capture: Dict[str, Any] = {"section": section, "path": str(path)}
        if reuse_screenshots:
            capture.update(png_metadata(path))
            capture["status"] = "pass" if capture.get("exists") and capture.get("bytes", 0) > 0 else "fail"
            capture["reused"] = True
            results.append(capture)
            continue
        try:
            locator.scroll_into_view_if_needed(timeout=30_000)
            locator.screenshot(
                path=str(path),
                animations="disabled",
                caret="hide",
                timeout=60_000,
            )
            capture.update(png_metadata(path))
            capture["status"] = "pass"
        except Exception as exc:
            capture.update({"status": "fail", "error": str(exc), "exists": False})
        results.append(capture)
    return results


def test_details(page: Page) -> Dict[str, Any]:
    results: List[Dict[str, Any]] = []
    details = page.locator("details")
    count = min(details.count(), DETAILS_LIMIT)
    for index in range(count):
        detail = details.nth(index)
        summary = detail.locator("summary").first
        before = detail.evaluate("el => el.open")
        text = ""
        try:
            detail.evaluate(
                """el => {
                    let parent = el.parentElement;
                    while (parent) {
                        if (parent.tagName === 'DETAILS' && !parent.open) {
                            parent.dataset.qaOpenedAncestor = 'true';
                            parent.open = true;
                        }
                        parent = parent.parentElement;
                    }
                }"""
            )
            text = summary.inner_text(timeout=5_000).strip()[:300]
            summary.scroll_into_view_if_needed(timeout=10_000)
            summary.click(timeout=10_000)
            page.wait_for_timeout(80)
            after = detail.evaluate("el => el.open")
            toggled = after != before
            summary.click(timeout=10_000)
            page.wait_for_timeout(50)
            restored = detail.evaluate("el => el.open") == before
            results.append(
                {
                    "index": index,
                    "summary": text,
                    "before_open": before,
                    "after_open": after,
                    "toggled": toggled,
                    "restored": restored,
                    "status": "pass" if toggled and restored else "fail",
                }
            )
        except Exception as exc:
            results.append(
                {
                    "index": index,
                    "summary": text,
                    "before_open": before,
                    "status": "fail",
                    "error": str(exc),
                }
            )
        finally:
            detail.evaluate(
                """el => {
                    let parent = el.parentElement;
                    while (parent) {
                        if (parent.tagName === 'DETAILS' && parent.dataset.qaOpenedAncestor === 'true') {
                            parent.open = false;
                            delete parent.dataset.qaOpenedAncestor;
                        }
                        parent = parent.parentElement;
                    }
                }"""
            )
    return {
        "tested": count,
        "available": details.count(),
        "passed": sum(item["status"] == "pass" for item in results),
        "failed": sum(item["status"] != "pass" for item in results),
        "results": results,
    }


def test_internal_navigation(page: Page) -> Dict[str, Any]:
    candidates = page.evaluate(
        """limit => Array.from(document.querySelectorAll('a[href^="#"]'))
            .map((el, index) => {
                const href = el.getAttribute('href') || '';
                let target = null;
                try { target = href.length > 1 ? document.querySelector(href) : null; } catch (_) {}
                const r = el.getBoundingClientRect();
                const s = getComputedStyle(el);
                const visible = s.display !== 'none' && s.visibility !== 'hidden' && r.width > 0 && r.height > 0;
                return {
                    dom_index: index,
                    href,
                    text: (el.innerText || el.getAttribute('aria-label') || '').replace(/\s+/g, ' ').trim().slice(0, 240),
                    visible,
                    is_skip_link: el.classList.contains('skip-link') || /跳到|skip/i.test(
                        el.innerText || el.getAttribute('aria-label') || ''),
                    target_exists: !!target
                };
            })
            .filter(item => item.visible && !item.is_skip_link && item.target_exists && item.href.length > 1)
            .slice(0, limit)""",
        NAV_LIMIT,
    )
    results: List[Dict[str, Any]] = []
    for candidate in candidates:
        href = candidate["href"]
        try:
            link = page.locator(f'a[href="{href}"]').first
            link.scroll_into_view_if_needed(timeout=10_000)
            link.click(timeout=10_000)
            page.wait_for_timeout(180)
            state = page.evaluate(
                """href => {
                    let target = null;
                    try { target = document.querySelector(href); } catch (_) {}
                    const r = target ? target.getBoundingClientRect() : null;
                    return {
                        hash: location.hash,
                        target_exists: !!target,
                        target_in_view: !!r && r.top < innerHeight && r.bottom > 0,
                        target_top: r ? Math.round(r.top * 10) / 10 : null
                    };
                }""",
                href,
            )
            passed = state["target_exists"] and state["target_in_view"]
            results.append({**candidate, "state": state, "status": "pass" if passed else "fail"})
        except Exception as exc:
            results.append({**candidate, "status": "fail", "error": str(exc)})
    page.evaluate("history.replaceState(null, '', location.pathname + location.search); scrollTo(0, 0)")
    return {
        "tested": len(results),
        "passed": sum(item["status"] == "pass" for item in results),
        "failed": sum(item["status"] != "pass" for item in results),
        "results": results,
    }


def test_keyboard_focus(page: Page) -> Dict[str, Any]:
    available = page.evaluate(
        """() => Array.from(document.querySelectorAll(
            'a[href],button,input,select,textarea,summary,[tabindex]:not([tabindex="-1"])'
        )).filter(el => {
            const closed = el.closest('details:not([open])');
            if (closed && el !== closed && !el.matches('summary') && !el.closest('summary')) return false;
            const r = el.getBoundingClientRect();
            const s = getComputedStyle(el);
            return !el.disabled && el.tabIndex >= 0 && s.display !== 'none' &&
                s.visibility !== 'hidden' && r.width > 0 && r.height > 0;
        }).length"""
    )
    page.evaluate(
        """() => {
            document.body.setAttribute('tabindex', '-1');
            document.body.focus();
            scrollTo(0, 0);
        }"""
    )
    results: List[Dict[str, Any]] = []
    for step in range(min(FOCUS_LIMIT, available)):
        page.keyboard.press("Tab")
        page.wait_for_timeout(220 if step == 0 else 55)
        state = page.evaluate(
            """step => {
                const el = document.activeElement;
                if (!el) return {step, missing: true};
                const r = el.getBoundingClientRect();
                const s = getComputedStyle(el);
                const name = (el.innerText || el.getAttribute('aria-label') ||
                    el.getAttribute('title') || el.value || '').replace(/\s+/g, ' ').trim().slice(0, 240);
                const outlineVisible = s.outlineStyle !== 'none' && parseFloat(s.outlineWidth || 0) > 0;
                const shadowVisible = s.boxShadow && s.boxShadow !== 'none';
                return {
                    step,
                    tag: el.tagName.toLowerCase(),
                    id: el.id || null,
                    role: el.getAttribute('role'),
                    name,
                    href: el.getAttribute('href'),
                    tabindex: el.getAttribute('tabindex'),
                    visible: r.width > 0 && r.height > 0 && s.display !== 'none' && s.visibility !== 'hidden',
                    in_view: r.top < innerHeight && r.bottom > 0,
                    focus_indicator: !!(outlineVisible || shadowVisible),
                    outline: `${s.outlineStyle} ${s.outlineWidth} ${s.outlineColor}`,
                    box_shadow: s.boxShadow,
                    box: {x: Math.round(r.x), y: Math.round(r.y), width: Math.round(r.width), height: Math.round(r.height)}
                };
            }""",
            step + 1,
        )
        results.append(state)
    page.evaluate("document.body.removeAttribute('tabindex'); document.activeElement?.blur(); scrollTo(0, 0)")
    identities = {
        (item.get("tag"), item.get("id"), item.get("name"), item.get("href"))
        for item in results
        if not item.get("missing")
    }
    return {
        "available": available,
        "steps": len(results),
        "unique_targets": len(identities),
        "visible_targets": sum(bool(item.get("visible")) for item in results),
        "in_view_targets": sum(bool(item.get("in_view")) for item in results),
        "with_focus_indicator": sum(bool(item.get("focus_indicator")) for item in results),
        "sequence": results,
    }


def viewport_audit(
    browser: Browser,
    source: Path,
    viewport_name: str,
    timeout_ms: int,
    reuse_screenshots: bool,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    viewport = VIEWPORTS[viewport_name]
    checks: List[Dict[str, Any]] = []
    result: Dict[str, Any] = {
        "name": viewport_name,
        "configured_viewport": viewport,
        "capture_order": [
            "01_navigate_and_wait_networkidle",
            "02_reuse_initial_viewport" if reuse_screenshots else "02_capture_initial_viewport",
            "03_reuse_full_page" if reuse_screenshots else "03_capture_full_page",
            "04_collect_dom_metrics",
            "05_reuse_main_data_sections" if reuse_screenshots else "05_capture_main_data_sections",
            "06_test_details_navigation_keyboard",
        ],
        "runtime": {"console_errors": [], "page_errors": [], "request_failures": []},
        "screenshots": {},
    }
    context = browser.new_context(
        viewport=viewport,
        color_scheme="light",
        reduced_motion="reduce",
        locale="zh-TW",
        device_scale_factor=1,
    )
    page = context.new_page()
    page.set_default_timeout(timeout_ms)

    page.on(
        "console",
        lambda message: result["runtime"]["console_errors"].append(
            {"text": message.text, "location": message.location}
        )
        if message.type == "error"
        else None,
    )
    page.on("pageerror", lambda error: result["runtime"]["page_errors"].append(str(error)))
    page.on(
        "requestfailed",
        lambda request: result["runtime"]["request_failures"].append(
            {"url": request.url, "failure": request.failure}
        ),
    )

    started = time.monotonic()
    navigation_started = time.monotonic()
    try:
        page.goto(source.as_uri(), wait_until="load", timeout=timeout_ms)
        page.wait_for_load_state("networkidle", timeout=timeout_ms)
        page.wait_for_timeout(1_500)
        try:
            page.wait_for_function(
                "document.querySelectorAll('[data-recharts-live] svg').length >= 14",
                timeout=min(timeout_ms, 15_000),
            )
        except Exception:
            # Runtime checks below retain the observable chart count; capture still proceeds.
            pass
        page.evaluate("scrollTo(0, 0)")
        result["load"] = {
            "status": "pass",
            "milliseconds": round((time.monotonic() - navigation_started) * 1000, 1),
            "ready_state": page.evaluate("document.readyState"),
            "url": page.url,
        }
    except Exception as exc:
        result["load"] = {
            "status": "fail",
            "milliseconds": round((time.monotonic() - navigation_started) * 1000, 1),
            "error": str(exc),
            "url": page.url,
        }
        add_check(checks, "document_load", False, "load + networkidle", result["load"], viewport=viewport_name)
        result["duration_seconds"] = round(time.monotonic() - started, 3)
        context.close()
        return result, checks

    # The first visual evidence is captured before any DOM interpretation or interaction.
    viewport_png = OUTPUT_DIR / f"{viewport_name}_viewport.png"
    full_png = OUTPUT_DIR / f"{viewport_name}_full_page.png"
    if reuse_screenshots:
        for name, path, full_page in (
            ("viewport", viewport_png, False),
            ("full_page", full_png, True),
        ):
            result["screenshots"][name] = png_metadata(path)
            result["screenshots"][name]["status"] = (
                "pass"
                if result["screenshots"][name].get("exists")
                and result["screenshots"][name].get("bytes", 0) > 0
                else "fail"
            )
            result["screenshots"][name]["full_page"] = full_page
            result["screenshots"][name]["reused"] = True
    else:
        result["screenshots"]["viewport"] = capture_png(page, viewport_png, full_page=False)
        result["screenshots"]["full_page"] = capture_png(page, full_png, full_page=True)

    result["dom"] = collect_dom_metrics(page)
    result["dom"]["heading_issues"] = heading_issues(result["dom"]["headings"])
    result["screenshots"]["sections"] = capture_sections(
        page,
        viewport_name,
        reuse_screenshots=reuse_screenshots,
    )
    result["interactions"] = {
        "details": test_details(page),
        "internal_navigation": test_internal_navigation(page),
        "keyboard_focus": test_keyboard_focus(page),
    }

    for kind in ("viewport", "full_page"):
        shot = result["screenshots"][kind]
        add_check(
            checks,
            f"{kind}_screenshot",
            shot["status"] == "pass" and shot.get("bytes", 0) > 0,
            "non-empty PNG",
            shot,
            viewport=viewport_name,
        )
    section_shots = result["screenshots"]["sections"]
    add_check(
        checks,
        "main_section_screenshots",
        bool(section_shots) and all(item["status"] == "pass" for item in section_shots),
        "at least one main data section and no capture failures",
        {
            "captured": len(section_shots),
            "failed": sum(item["status"] != "pass" for item in section_shots),
        },
        viewport=viewport_name,
    )
    add_check(checks, "document_load", True, "load + networkidle", result["load"], viewport=viewport_name)
    for error_type in ("console_errors", "page_errors", "request_failures"):
        errors = result["runtime"][error_type]
        add_check(checks, error_type, not errors, {"count": 0}, {"count": len(errors), "items": errors}, viewport=viewport_name)
    body_overflow = result["dom"]["document"]["body_overflow_x_px"]
    add_check(checks, "body_horizontal_overflow", body_overflow <= 1, {"maximum_px": 1}, body_overflow, viewport=viewport_name)
    issues = result["dom"]["heading_issues"]
    add_check(checks, "heading_hierarchy", not issues, "one H1 and no skipped levels", issues, viewport=viewport_name)
    missing_targets = result["dom"]["navigation"]["missing_internal_targets"]
    add_check(checks, "internal_navigation_targets", not missing_targets, {"missing": 0}, missing_targets, viewport=viewport_name)
    broken_images = result["dom"]["broken_images"]
    add_check(checks, "broken_images", not broken_images, {"count": 0}, broken_images, viewport=viewport_name)
    source_disclosure = result["dom"]["source_disclosure"]
    add_check(
        checks,
        "inline_json_source_clutter",
        source_disclosure["inline_source_tooltips"] == 0 and not source_disclosure["visible_json_filename_mentions"],
        {"inline_source_tooltips": 0, "visible_json_filename_mentions": []},
        {
            "inline_source_tooltips": source_disclosure["inline_source_tooltips"],
            "visible_json_filename_mentions": source_disclosure["visible_json_filename_mentions"],
        },
        viewport=viewport_name,
    )
    missing_companions = []
    for href in source_disclosure["unique_companion_hrefs"]:
        target = (source.parent / href).resolve()
        if not target.is_file():
            missing_companions.append({"href": href, "resolved": str(target)})
    add_check(
        checks,
        "hidden_companion_json_links",
        len(source_disclosure["unique_companion_hrefs"]) >= 4 and not missing_companions,
        {"minimum_unique_links": 4, "missing": 0},
        {"unique_links": source_disclosure["unique_companion_hrefs"], "missing": missing_companions},
        viewport=viewport_name,
    )
    details_result = result["interactions"]["details"]
    add_check(
        checks,
        "details_toggle_restore",
        details_result["failed"] == 0,
        {"failed": 0},
        {"tested": details_result["tested"], "failed": details_result["failed"]},
        viewport=viewport_name,
    )
    nav_result = result["interactions"]["internal_navigation"]
    add_check(
        checks,
        "internal_navigation_interaction",
        nav_result["tested"] > 0 and nav_result["failed"] == 0,
        "at least one working internal navigation link",
        {"tested": nav_result["tested"], "failed": nav_result["failed"]},
        viewport=viewport_name,
    )
    focus = result["interactions"]["keyboard_focus"]
    add_check(
        checks,
        "keyboard_focus_progression",
        focus["unique_targets"] >= min(5, result["dom"]["counts"]["focusables"]),
        {"minimum_unique": min(5, result["dom"]["counts"]["focusables"])},
        {
            "unique_targets": focus["unique_targets"],
            "steps": focus["steps"],
            "with_focus_indicator": focus["with_focus_indicator"],
        },
        viewport=viewport_name,
    )
    skip_links = result["dom"]["navigation"]["skip_links"]
    if skip_links:
        first_focus = focus["sequence"][0] if focus["sequence"] else None
        skip_focus_passed = bool(
            first_focus
            and first_focus.get("href") == skip_links[0].get("href")
            and first_focus.get("in_view")
            and first_focus.get("focus_indicator")
        )
        add_check(
            checks,
            "skip_link_keyboard_visibility",
            skip_focus_passed,
            "first Tab focuses an in-view skip link with a visible indicator",
            first_focus,
            viewport=viewport_name,
        )
    missing_indicators = [
        item
        for item in focus["sequence"]
        if item.get("visible") and item.get("tag") != "body" and not item.get("focus_indicator")
    ]
    add_check(
        checks,
        "keyboard_focus_indicators",
        not missing_indicators,
        {"missing": 0},
        {"missing": len(missing_indicators), "targets": missing_indicators},
        viewport=viewport_name,
    )
    result["duration_seconds"] = round(time.monotonic() - started, 3)
    context.close()
    return result, checks


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Capture desktop/mobile before-baseline screenshots and UI metrics for "
            "the standalone layered reconstruction panorama report."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Exit codes: 0 all checks pass, 1 findings remain, 2 runner/configuration error.",
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Standalone HTML input.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Output directory.")
    parser.add_argument("--metrics", type=Path, default=DEFAULT_METRICS, help="Metrics JSON path.")
    parser.add_argument("--timeout-ms", type=int, default=60_000, help="Playwright timeout in milliseconds.")
    parser.add_argument(
        "--reuse-screenshots",
        action="store_true",
        help="Refresh metrics/interactions while reusing existing PNG evidence.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    global OUTPUT_DIR
    args = build_parser().parse_args(argv)
    source = args.input.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    metrics_path = args.metrics.expanduser().resolve()
    started_at = utc_now()

    if output_dir.parent != AUDIT_ROOT:
        print(json.dumps({"status": "error", "error": f"Output must remain directly under {AUDIT_ROOT}"}, ensure_ascii=False))
        return 2
    if metrics_path.parent != output_dir:
        print(json.dumps({"status": "error", "error": f"Metrics must remain under {output_dir}"}, ensure_ascii=False))
        return 2
    if not source.is_file():
        print(json.dumps({"status": "error", "error": f"Input not found: {source}"}, ensure_ascii=False))
        return 2
    if args.timeout_ms <= 0:
        print(json.dumps({"status": "error", "error": "--timeout-ms must be positive"}, ensure_ascii=False))
        return 2
    if PLAYWRIGHT_IMPORT_ERROR:
        print(json.dumps({"status": "error", "error": PLAYWRIGHT_IMPORT_ERROR}, ensure_ascii=False))
        return 2

    output_dir.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR = output_dir
    source_hash_before = sha256(source)
    report: Dict[str, Any] = {
        "schema_version": "1.0",
        "audit_kind": f"{output_dir.name}_visual_qa",
        "generated_at": started_at,
        "source": {
            "path": str(source),
            "url": source.as_uri(),
            "bytes": source.stat().st_size,
            "mtime": datetime.fromtimestamp(source.stat().st_mtime, timezone.utc).isoformat(),
            "sha256_before": source_hash_before,
            "sha256_after": None,
            "unchanged": None,
        },
        "output": {"directory": str(output_dir), "metrics": str(metrics_path)},
        "browser": None,
        "viewports": [],
        "checks": [],
        "fatal_errors": [],
    }

    playwright = None
    browser = None
    try:
        playwright = sync_playwright().start()  # type: ignore[union-attr]
        browser = playwright.chromium.launch(headless=True, args=["--allow-file-access-from-files"])
        report["browser"] = {"engine": "chromium", "version": browser.version, "headless": True}
        for viewport_name in VIEWPORTS:
            viewport_result, viewport_checks = viewport_audit(
                browser,
                source,
                viewport_name,
                args.timeout_ms,
                args.reuse_screenshots,
            )
            report["viewports"].append(viewport_result)
            report["checks"].extend(viewport_checks)
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

    source_hash_after = sha256(source)
    report["source"]["sha256_after"] = source_hash_after
    report["source"]["unchanged"] = source_hash_before == source_hash_after
    add_check(
        report["checks"],
        "source_html_unchanged",
        report["source"]["unchanged"],
        source_hash_before,
        source_hash_after,
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
    report["summary"] = {
        "status": status,
        "checks_total": len(report["checks"]),
        "checks_passed": len(report["checks"]) - len(failed),
        "checks_failed": len(failed),
        "failed_checks": [
            {"name": item["name"], "viewport": item.get("viewport"), "actual": item["actual"]}
            for item in failed
        ],
        "exit_code": exit_code,
    }
    report["finished_at"] = utc_now()
    metrics_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": status,
                "exit_code": exit_code,
                "metrics": str(metrics_path),
                "checks_total": report["summary"]["checks_total"],
                "checks_failed": report["summary"]["checks_failed"],
                "screenshots": sum(
                    2 + len(viewport.get("screenshots", {}).get("sections", []))
                    for viewport in report["viewports"]
                ),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
