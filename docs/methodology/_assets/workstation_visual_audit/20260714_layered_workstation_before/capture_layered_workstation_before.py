#!/usr/bin/env python3
"""Capture and measure the complete layered workstation before baseline.

The runner opens the static HTML files with local headless Chromium.  It writes
only under the audit output directory and verifies that all source HTML hashes
remain unchanged.  A full-page screenshot is always captured before DOM
interpretation; sample interactions then select a real tree-bearing region and
capture the region browser state.

Exit codes:
    0: capture completed and all checks passed
    1: capture completed and one or more audit findings remain
    2: runner or configuration failure
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
SPECIAL_TREE_SAMPLE = "HCC1395.html"
VIEWPORTS = {
    "desktop": {"width": 1440, "height": 1000},
    "mobile": {"width": 390, "height": 844},
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

    for viewport in VIEWPORTS:
        key = (INDEX_FILE, viewport, "full_page")
        plan[key] = {
            "step": step,
            "description": f"cohort index · {viewport} · initial full page",
            "filename": f"{step:02d}_{viewport}_index_full_page.png",
        }
        step += 1

    for sample_file in SAMPLE_FILES:
        sample = sample_file[:-5]
        for viewport in VIEWPORTS:
            for kind, description in (
                ("full_page", "initial full page"),
                ("region_view", "selected tree-bearing region browser"),
            ):
                key = (sample_file, viewport, kind)
                plan[key] = {
                    "step": step,
                    "description": f"{sample} · {viewport} · {description}",
                    "filename": f"{step:02d}_{viewport}_{slug(sample_file)}_{kind}.png",
                }
                step += 1
            if sample_file == SPECIAL_TREE_SAMPLE:
                key = (sample_file, viewport, "candidate_tree")
                plan[key] = {
                    "step": step,
                    "description": f"{sample} · {viewport} · expanded candidate tree/network detail",
                    "filename": f"{step:02d}_{viewport}_{slug(sample_file)}_candidate_tree.png",
                }
                step += 1

    if step != 33:
        raise AssertionError(f"Expected 32 screenshot steps, got {step - 1}")
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


def capture_page_screenshot(
    page: Page,
    path: Path,
    *,
    full_page: bool,
    reuse: bool,
) -> Dict[str, Any]:
    if reuse:
        result = png_metadata(path)
        result["reused"] = True
        result["status"] = "pass" if result.get("bytes", 0) > 0 else "fail"
        return result
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


def capture_locator_screenshot(
    locator: Locator,
    path: Path,
    *,
    reuse: bool,
) -> Dict[str, Any]:
    if reuse:
        result = png_metadata(path)
        result["reused"] = True
        result["status"] = "pass" if result.get("bytes", 0) > 0 else "fail"
        return result
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


def collect_dom_metrics(page: Page, page_kind: str) -> Dict[str, Any]:
    return page.evaluate(
        """pageKind => {
            const all = selector => Array.from(document.querySelectorAll(selector));
            const short = (value, n=240) => String(value || '').replace(/\s+/g, ' ').trim().slice(0, n);
            const visible = el => {
                const closed = el.closest('details:not([open])');
                if (closed && el !== closed && !el.matches('summary') && !el.closest('summary')) return false;
                const style = getComputedStyle(el);
                const rect = el.getBoundingClientRect();
                return style.display !== 'none' && style.visibility !== 'hidden' &&
                    Number(style.opacity || 1) !== 0 && rect.width > 0 && rect.height > 0;
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
            const nearestScroller = el => {
                let p = el.parentElement;
                while (p && p !== document.body) {
                    const s = getComputedStyle(p);
                    if (/(auto|scroll)/.test(s.overflowX) && p.scrollWidth > p.clientWidth + 1) return p;
                    p = p.parentElement;
                }
                return null;
            };
            const root = document.documentElement;
            const body = document.body;
            const viewportWidth = innerWidth;
            const headings = all('h1,h2,h3,h4,h5,h6').filter(visible).map((el, index) => ({
                index,
                level: Number(el.tagName.slice(1)),
                tag: el.tagName.toLowerCase(),
                id: el.id || null,
                text: short(el.innerText, 500),
                box: box(el)
            }));
            const links = all('a[href]').map((el, index) => {
                const href = el.getAttribute('href') || '';
                let targetExists = null;
                if (href.startsWith('#')) {
                    try { targetExists = href.length > 1 && !!document.querySelector(href); }
                    catch (_) { targetExists = false; }
                }
                return {
                    index,
                    text: short(el.innerText || el.getAttribute('aria-label')),
                    href,
                    resolved: el.href,
                    target_exists: targetExists,
                    visible: visible(el),
                    box: visible(el) ? box(el) : null
                };
            });
            const tables = all('table').filter(visible).map((el, index) => {
                const r = el.getBoundingClientRect();
                const scroller = nearestScroller(el);
                return {
                    index,
                    rows: el.rows ? el.rows.length : 0,
                    columns_max: Array.from(el.rows || []).reduce((n, row) => Math.max(n, row.cells.length), 0),
                    box: box(el),
                    exceeds_viewport: r.width > innerWidth + 1,
                    intentional_horizontal_scroller: !!scroller,
                    scroller: scroller ? {
                        tag: scroller.tagName.toLowerCase(),
                        id: scroller.id || null,
                        class: typeof scroller.className === 'string' ? short(scroller.className, 120) : null,
                        client_width: scroller.clientWidth,
                        scroll_width: scroller.scrollWidth
                    } : null
                };
            });
            const localScrollers = all('body *').filter(el => {
                if (!visible(el)) return false;
                const s = getComputedStyle(el);
                return /(auto|scroll)/.test(s.overflowX) && el.scrollWidth > el.clientWidth + 1;
            }).map((el, index) => ({
                index,
                tag: el.tagName.toLowerCase(),
                id: el.id || null,
                class: typeof el.className === 'string' ? short(el.className, 120) : null,
                client_width: el.clientWidth,
                scroll_width: el.scrollWidth,
                overflow_px: el.scrollWidth - el.clientWidth,
                box: box(el)
            }));
            const overflowOffenders = all('body *').filter(el => {
                if (!visible(el) || el.matches('.skip,.skip-link')) return false;
                const r = el.getBoundingClientRect();
                const offscreen = r.left < -1 || r.right > viewportWidth + 1;
                if (!offscreen) return false;
                const parent = el.parentElement;
                if (!parent || parent === document.body) return true;
                const pr = parent.getBoundingClientRect();
                const parentOffscreen = pr.left < -1 || pr.right > viewportWidth + 1;
                return !parentOffscreen;
            }).slice(0, 100).map(el => {
                const scroller = nearestScroller(el);
                return {
                    tag: el.tagName.toLowerCase(),
                    id: el.id || null,
                    class: typeof el.className === 'string' ? short(el.className, 120) : null,
                    box: box(el),
                    intentional_local_scroll: !!scroller,
                    scroller_id: scroller ? (scroller.id || null) : null,
                    scroller_class: scroller && typeof scroller.className === 'string' ? short(scroller.className, 120) : null
                };
            });
            const focusSelector = 'a[href],button,input,select,textarea,summary,[tabindex]:not([tabindex="-1"])';
            const focusables = all(focusSelector).filter(el => visible(el) && !el.disabled && el.tabIndex >= 0).map((el, index) => {
                const r = el.getBoundingClientRect();
                return {
                    index,
                    tag: el.tagName.toLowerCase(),
                    id: el.id || null,
                    role: el.getAttribute('role'),
                    name: short(el.innerText || el.getAttribute('aria-label') || el.value),
                    width: Math.round(r.width * 10) / 10,
                    height: Math.round(r.height * 10) / 10,
                    below_24px_target: r.width < 24 || r.height < 24,
                    below_44px_target: r.width < 44 || r.height < 44
                };
            });
            const treeSvgs = all('.tsw svg').filter(visible).map((svg, index) => {
                const r = svg.getBoundingClientRect();
                const iw = svg.width && svg.width.baseVal ? svg.width.baseVal.value : Number(svg.getAttribute('width')) || r.width;
                const ih = svg.height && svg.height.baseVal ? svg.height.baseVal.value : Number(svg.getAttribute('height')) || r.height;
                const scale = iw ? r.width / iw : 1;
                const stage = svg.closest('.tstage');
                return {
                    index,
                    box: box(svg),
                    intrinsic_width: Math.round(iw * 10) / 10,
                    intrinsic_height: Math.round(ih * 10) / 10,
                    scale_ratio: Math.round(scale * 1000) / 1000,
                    estimated_label_px: Math.round(9 * scale * 10) / 10,
                    text_nodes: svg.querySelectorAll('text').length,
                    edge_lines: svg.querySelectorAll('line').length,
                    nodes: svg.querySelectorAll('circle,rect').length,
                    view_box: svg.getAttribute('viewBox'),
                    stage_overflow_x: stage ? getComputedStyle(stage).overflowX : null,
                    stage_scrollable: !!stage && stage.scrollWidth > stage.clientWidth + 1,
                    clipped: !!stage && r.right > stage.getBoundingClientRect().right + 1
                };
            });
            const mainText = document.querySelector('main')?.innerText || '';
            const chrControl = document.getElementById('fchr');
            const chrScopeText = chrControl ? short(
                `${document.querySelector('label[for="fchr"]')?.innerText || ''} ${chrControl.selectedOptions?.[0]?.textContent || ''}`
            ) : '';
            const regionSignals = Array.from(mainText.matchAll(/Region-level denominator|regions? · 本頁主分母|主分母|eligible regions?|Region 瀏覽器/gi))
                .map(match => match[0]);
            const explicitGenomeSignals = Array.from(mainText.matchAll(/全基因組|whole[- ]genome|genome[- ]wide|全染色體/gi))
                .map(match => match[0]);
            return {
                page_kind: pageKind,
                document: {
                    title: document.title,
                    lang: root.lang || null,
                    ready_state: document.readyState,
                    viewport: {width: innerWidth, height: innerHeight, device_pixel_ratio: devicePixelRatio},
                    page: {width: root.scrollWidth, height: root.scrollHeight},
                    body_overflow_x_px: Math.max(0, Math.max(root.scrollWidth, body.scrollWidth) - innerWidth),
                    body_overflow_x_css: getComputedStyle(body).overflowX,
                    overflow_offenders: overflowOffenders,
                    unintended_overflow_offenders: overflowOffenders.filter(item => !item.intentional_local_scroll)
                },
                counts: {
                    elements: all('*').length,
                    h1: headings.filter(item => item.level === 1).length,
                    h2: headings.filter(item => item.level === 2).length,
                    h3: headings.filter(item => item.level === 3).length,
                    headings: headings.length,
                    main: all('main').filter(visible).length,
                    sections: all('section').filter(visible).length,
                    tables: tables.length,
                    table_rows: tables.reduce((sum, item) => sum + item.rows, 0),
                    links: links.filter(item => item.visible).length,
                    buttons: all('button').filter(visible).length,
                    inputs_selects: all('input,select,textarea').filter(visible).length,
                    details: all('details').filter(visible).length,
                    details_open: all('details[open]').filter(visible).length,
                    focusables: focusables.length,
                    svg: all('svg').filter(visible).length,
                    tree_svg: treeSvgs.length
                },
                headings,
                links,
                tables,
                local_horizontal_scrollers: localScrollers,
                focusables: {
                    items: focusables,
                    below_24px_target: focusables.filter(item => item.below_24px_target).length,
                    below_44px_target: focusables.filter(item => item.below_44px_target).length
                },
                scope_labels: {
                    region_signals: regionSignals,
                    region_scope_visible: regionSignals.length > 0,
                    chromosome_control_text: chrScopeText,
                    all_chromosomes_visible: /染色體\s*全部/.test(chrScopeText),
                    explicit_whole_genome_signals: explicitGenomeSignals,
                    explicit_whole_genome_visible: explicitGenomeSignals.length > 0
                },
                tree_svgs: treeSvgs,
                sample_switch_options: all('#sample-switch option').map(option => ({
                    text: short(option.textContent),
                    value: option.value,
                    selected: option.selected
                })),
                broken_images: all('img').filter(img => !img.complete || img.naturalWidth === 0).map(img => ({
                    src: img.getAttribute('src'),
                    alt: img.getAttribute('alt')
                }))
            };
        }""",
        page_kind,
    )


def validate_link_targets(dom: Dict[str, Any], source_dir: Path) -> Dict[str, Any]:
    results: List[Dict[str, Any]] = []
    for link in dom["links"]:
        href = link["href"]
        item = {"text": link["text"], "href": href, "resolved": link["resolved"]}
        if href.startswith("#"):
            item["kind"] = "internal_anchor"
            item["exists"] = bool(link["target_exists"])
        else:
            parsed = urlparse(link["resolved"])
            if parsed.scheme == "file":
                target = Path(unquote(parsed.path))
                item["kind"] = "local_file"
                item["target"] = str(target)
                item["exists"] = target.is_file()
            else:
                item["kind"] = parsed.scheme or "relative"
                item["exists"] = True
        results.append(item)

    for option in dom.get("sample_switch_options", []):
        target = (source_dir / option["value"]).resolve()
        results.append(
            {
                "text": f"sample-switch: {option['text']}",
                "href": option["value"],
                "kind": "select_option_local_file",
                "target": str(target),
                "exists": target.is_file(),
            }
        )
    return {
        "tested": len(results),
        "passed": sum(bool(item["exists"]) for item in results),
        "failed": sum(not bool(item["exists"]) for item in results),
        "results": results,
    }


def choose_candidate_region(page: Page) -> Dict[str, Any]:
    candidate = page.evaluate(
        """() => {
            const leafCount = edges => {
                const parents = new Set((edges || []).map(edge => edge[0]));
                const nodes = new Set();
                (edges || []).forEach(edge => { nodes.add(edge[0]); nodes.add(edge[1]); });
                return Array.from(nodes).filter(node => !parents.has(node)).length;
            };
            const candidates = R.map((region, index) => {
                let best = null;
                (region.lineages || []).forEach((lineage, lineageIndex) => {
                    const visualTrees = (lineage.trees || []).filter(tree => (tree.edges || []).length > 0);
                    if (!visualTrees.length) return;
                    const maxEdges = Math.max(...visualTrees.map(tree => tree.edges.length));
                    const maxLeaves = Math.max(...visualTrees.map(tree => leafCount(tree.edges)));
                    const treeCount = visualTrees.length;
                    const preferredCount = treeCount >= 2 && treeCount <= 8;
                    const score = (preferredCount ? 100000 : 0) + maxLeaves * 1000 + maxEdges * 50 + Math.min(treeCount, 20);
                    if (!best || score > best.score) {
                        best = {lineage_index: lineageIndex, tree_count: treeCount, max_edges: maxEdges, max_leaves: maxLeaves, score};
                    }
                });
                return best ? {
                    index,
                    region: region.region,
                    chrom: region.chrom,
                    hp_multiplicity: region.hp_multiplicity,
                    determinacy: region.region_determinacy,
                    n_sSNV: region.n_sSNV,
                    ...best
                } : null;
            }).filter(Boolean).sort((a, b) => b.score - a.score || a.index - b.index);
            return candidates[0] || null;
        }"""
    )
    if not candidate:
        raise RuntimeError("No tree-bearing region was found")

    page.locator("#fq").fill(candidate["region"])
    row = page.locator(f'#list .row[data-i="{candidate["index"]}"]')
    row.wait_for(state="visible", timeout=30_000)
    row.click(timeout=30_000)
    page.wait_for_function(
        "region => document.querySelector('#detail h3')?.textContent.includes(region)",
        arg=candidate["region"],
        timeout=30_000,
    )
    page.wait_for_selector("#detail .tsw svg", state="visible", timeout=30_000)
    page.wait_for_timeout(250)
    candidate["selection_verified"] = True
    candidate["rendered_tree_blocks"] = page.locator("#detail .tsw").count()
    candidate["rendered_tree_svgs"] = page.locator("#detail .tsw svg").count()
    return candidate


def mark_largest_tree_block(page: Page) -> Dict[str, Any]:
    return page.evaluate(
        """() => {
            const blocks = Array.from(document.querySelectorAll('#detail .tsw'));
            blocks.forEach(block => block.removeAttribute('data-audit-largest-tree'));
            const ranked = blocks.map((block, index) => {
                const svg = block.querySelector('.tslide:not([style*="display: none"]) svg, .tslide svg');
                const width = svg?.width?.baseVal?.value || Number(svg?.getAttribute('width')) || 0;
                const edges = svg ? svg.querySelectorAll('line').length : 0;
                return {block, index, width, edges, score: width + edges * 20};
            }).sort((a, b) => b.score - a.score);
            if (!ranked.length) return null;
            const selected = ranked[0];
            selected.block.setAttribute('data-audit-largest-tree', 'true');
            const thumbs = selected.block.querySelectorAll('.thumb');
            if (thumbs.length > 1) thumbs[1].click();
            const svg = selected.block.querySelector('.tslide:not([style*="display: none"]) svg, .tslide svg');
            return {
                block_index: selected.index,
                intrinsic_width: svg?.width?.baseVal?.value || Number(svg?.getAttribute('width')) || null,
                intrinsic_height: svg?.height?.baseVal?.value || Number(svg?.getAttribute('height')) || null,
                tree_number: thumbs.length > 1 ? 2 : 1,
                thumbs: thumbs.length
            };
        }"""
    )


def test_skip_link(page: Page) -> Dict[str, Any]:
    page.evaluate(
        """() => {
            document.body.setAttribute('tabindex', '-1');
            document.body.focus();
            document.documentElement.dataset.auditPreviousScrollBehavior = document.documentElement.style.scrollBehavior || '';
            document.documentElement.style.scrollBehavior = 'auto';
            scrollTo(0, 0);
        }"""
    )
    page.wait_for_timeout(80)
    page.keyboard.press("Tab")
    page.wait_for_timeout(100)
    focused = page.evaluate(
        """() => {
            const el = document.activeElement;
            if (!el) return null;
            const r = el.getBoundingClientRect();
            const s = getComputedStyle(el);
            return {
                tag: el.tagName.toLowerCase(),
                class: typeof el.className === 'string' ? el.className : null,
                text: (el.innerText || el.getAttribute('aria-label') || '').replace(/\s+/g, ' ').trim(),
                href: el.getAttribute('href'),
                in_view: r.top < innerHeight && r.bottom > 0 && r.left < innerWidth && r.right > 0,
                outline: `${s.outlineStyle} ${s.outlineWidth} ${s.outlineColor}`,
                box_shadow: s.boxShadow,
                indicator: (s.outlineStyle !== 'none' && parseFloat(s.outlineWidth || 0) > 0) || (s.boxShadow && s.boxShadow !== 'none')
            };
        }"""
    )
    target_state: Dict[str, Any] = {}
    if focused and focused.get("href", "").startswith("#"):
        page.keyboard.press("Enter")
        page.wait_for_timeout(160)
        target_state = page.evaluate(
            """href => {
                let target = null;
                try { target = document.querySelector(href); } catch (_) {}
                const r = target?.getBoundingClientRect();
                return {
                    hash: location.hash,
                    target_exists: !!target,
                    target_in_view: !!r && r.top < innerHeight && r.bottom > 0,
                    target_top: r ? Math.round(r.top * 10) / 10 : null
                };
            }""",
            focused["href"],
        )
    page.evaluate(
        """() => {
            document.body.removeAttribute('tabindex');
            document.activeElement?.blur();
            history.replaceState(null, '', location.pathname + location.search + location.hash.replace(location.hash, ''));
            scrollTo(0, 0);
            document.documentElement.style.scrollBehavior = document.documentElement.dataset.auditPreviousScrollBehavior || '';
            delete document.documentElement.dataset.auditPreviousScrollBehavior;
        }"""
    )
    passed = bool(
        focused
        and focused.get("href", "").startswith("#")
        and focused.get("in_view")
        and focused.get("indicator")
        and target_state.get("target_exists")
        and target_state.get("target_in_view")
    )
    return {"status": "pass" if passed else "fail", "focused": focused, "target": target_state}


def test_keyboard_focus(page: Page, limit: int = 120) -> Dict[str, Any]:
    available = page.evaluate(
        """() => Array.from(document.querySelectorAll(
            'a[href],button,input,select,textarea,summary,[tabindex]:not([tabindex="-1"])'
        )).filter(el => {
            const closed = el.closest('details:not([open])');
            if (closed && el !== closed && !el.matches('summary') && !el.closest('summary')) return false;
            const r = el.getBoundingClientRect();
            const s = getComputedStyle(el);
            return !el.disabled && el.tabIndex >= 0 && s.display !== 'none' && s.visibility !== 'hidden' && r.width > 0 && r.height > 0;
        }).length"""
    )
    page.evaluate(
        """() => {
            document.body.setAttribute('tabindex', '-1');
            document.body.focus();
            scrollTo(0, 0);
        }"""
    )
    sequence: List[Dict[str, Any]] = []
    for step in range(min(available + 1, limit)):
        page.keyboard.press("Tab")
        page.wait_for_timeout(30)
        state = page.evaluate(
            """step => {
                const el = document.activeElement;
                if (!el) return {step, missing: true};
                const r = el.getBoundingClientRect();
                const s = getComputedStyle(el);
                const name = (el.innerText || el.getAttribute('aria-label') || el.value || '').replace(/\s+/g, ' ').trim().slice(0, 220);
                const outline = s.outlineStyle !== 'none' && parseFloat(s.outlineWidth || 0) > 0;
                const shadow = s.boxShadow && s.boxShadow !== 'none';
                return {
                    step,
                    tag: el.tagName.toLowerCase(),
                    id: el.id || null,
                    class: typeof el.className === 'string' ? el.className.slice(0, 120) : null,
                    role: el.getAttribute('role'),
                    name,
                    href: el.getAttribute('href'),
                    visible: r.width > 0 && r.height > 0 && s.display !== 'none' && s.visibility !== 'hidden',
                    in_view: r.top < innerHeight && r.bottom > 0,
                    focus_indicator: !!(outline || shadow),
                    outline: `${s.outlineStyle} ${s.outlineWidth} ${s.outlineColor}`,
                    box_shadow: s.boxShadow,
                    box: {x: Math.round(r.x), y: Math.round(r.y), width: Math.round(r.width), height: Math.round(r.height)}
                };
            }""",
            step + 1,
        )
        sequence.append(state)
        if step > 0 and state.get("tag") == "body":
            break
    page.evaluate("document.body.removeAttribute('tabindex'); document.activeElement?.blur(); scrollTo(0, 0)")
    real = [item for item in sequence if item.get("tag") != "body" and not item.get("missing")]
    missing = [item for item in real if item.get("visible") and not item.get("focus_indicator")]
    identities = {(item.get("tag"), item.get("id"), item.get("class"), item.get("name")) for item in real}
    return {
        "available": available,
        "steps": len(sequence),
        "unique_targets": len(identities),
        "with_indicator": sum(bool(item.get("focus_indicator")) for item in real),
        "missing_indicator_count": len(missing),
        "missing_indicators": missing,
        "sequence": sequence,
    }


def test_details(page: Page) -> Dict[str, Any]:
    results: List[Dict[str, Any]] = []
    details = page.locator("details")
    for index in range(details.count()):
        detail = details.nth(index)
        summary = detail.locator("summary").first
        try:
            before = detail.evaluate("el => el.open")
            text = summary.inner_text(timeout=5_000).strip()[:240]
            summary.click(timeout=10_000)
            page.wait_for_timeout(60)
            after = detail.evaluate("el => el.open")
            summary.click(timeout=10_000)
            restored = detail.evaluate("el => el.open") == before
            passed = after != before and restored
            results.append({"index": index, "summary": text, "toggled": after != before, "restored": restored, "status": "pass" if passed else "fail"})
        except Exception as exc:
            results.append({"index": index, "status": "fail", "error": str(exc)})
    return {
        "tested": len(results),
        "passed": sum(item["status"] == "pass" for item in results),
        "failed": sum(item["status"] != "pass" for item in results),
        "results": results,
    }


def select_and_measure(page: Page, selector: str, value: str) -> Dict[str, Any]:
    before = page.locator(selector).input_value()
    page.locator(selector).select_option(value)
    page.wait_for_timeout(80)
    actual = page.locator(selector).input_value()
    count_text = page.locator("#cnt").inner_text().strip() if page.locator("#cnt").count() else None
    page.locator(selector).select_option(before)
    return {"selector": selector, "requested": value, "actual": actual, "count_text": count_text, "restored": page.locator(selector).input_value() == before, "status": "pass" if actual == value else "fail"}


def test_sample_controls(page: Page, candidate: Dict[str, Any]) -> Dict[str, Any]:
    tests: List[Dict[str, Any]] = []

    overview = page.locator("details.overview")
    before_open = overview.evaluate("el => el.open")
    overview.locator("summary").click()
    opened = overview.evaluate("el => el.open")
    overview.locator("summary").click()
    restored = overview.evaluate("el => el.open") == before_open
    tests.append({"control": "overview details", "opened": opened, "restored": restored, "status": "pass" if opened != before_open and restored else "fail"})

    copy_button = page.locator("#copy-link")
    before_text = copy_button.inner_text().strip()
    copy_button.click()
    page.wait_for_timeout(80)
    after_text = copy_button.inner_text().strip()
    copy_ok = after_text in {"已複製 Region 連結", "請複製網址列"}
    tests.append({"control": "copy region link", "before": before_text, "after": after_text, "status": "pass" if copy_ok else "fail"})

    tests.append({"control": "determinacy filter", **select_and_measure(page, "#fdet", candidate["determinacy"])})
    tests.append({"control": "HP multiplicity filter", **select_and_measure(page, "#fhp", str(candidate["hp_multiplicity"]))})
    tests.append({"control": "chromosome filter", **select_and_measure(page, "#fchr", candidate["chrom"])})
    tests.append({"control": "sort", **select_and_measure(page, "#fsort", "ntree")})

    query = page.locator("#fq")
    query.fill("")
    page.wait_for_timeout(80)
    cleared_count = page.locator("#cnt").inner_text().strip()
    query.fill(candidate["region"])
    page.wait_for_timeout(80)
    exact_rows = page.locator("#list .row").count()
    tests.append({"control": "region search", "cleared_count": cleared_count, "exact_rows": exact_rows, "status": "pass" if exact_rows >= 1 else "fail"})

    row = page.locator(f'#list .row[data-i="{candidate["index"]}"]')
    row.click()
    detail_title = page.locator("#detail h3").inner_text().strip()
    tests.append({"control": "region row selection", "detail_title": detail_title, "status": "pass" if candidate["region"] in detail_title else "fail"})

    sample_switch = page.locator("#sample-switch")
    selected_text = sample_switch.locator("option:checked").inner_text().strip()
    tests.append({"control": "sample switch", "options": sample_switch.locator("option").count(), "selected": selected_text, "status": "pass" if sample_switch.locator("option").count() == 7 else "fail"})

    next_button = page.locator('#detail .tnav[aria-label="下一棵等機率樹"]').first
    if next_button.count():
        block = next_button.locator("xpath=ancestor::div[contains(@class,'tsw')]")
        counter = block.locator(".tctr")
        before_counter = counter.inner_text().strip()
        next_button.click()
        after_counter = counter.inner_text().strip()
        prev_button = block.locator('.tnav[aria-label="上一棵等機率樹"]')
        prev_button.click()
        restored_counter = counter.inner_text().strip()
        tests.append({"control": "tree next/previous", "before": before_counter, "after": after_counter, "restored": restored_counter, "status": "pass" if after_counter != before_counter and restored_counter == before_counter else "fail"})
    else:
        tests.append({"control": "tree next/previous", "status": "not_applicable", "reason": "candidate has one rendered tree"})

    thumb = page.locator("#detail .tsw .thumb").nth(1)
    if thumb.count():
        block = thumb.locator("xpath=ancestor::div[contains(@class,'tsw')]")
        thumb.click()
        selected_second = thumb.evaluate("el => el.classList.contains('on')")
        first = block.locator(".thumb").first
        first.click()
        restored_first = first.evaluate("el => el.classList.contains('on')")
        tests.append({"control": "tree thumbnail", "selected_second": selected_second, "restored_first": restored_first, "status": "pass" if selected_second and restored_first else "fail"})
    else:
        tests.append({"control": "tree thumbnail", "status": "not_applicable", "reason": "candidate has one rendered tree"})

    details = test_details(page)
    tests.append({"control": "all details toggles", "tested": details["tested"], "failed": details["failed"], "status": "pass" if details["failed"] == 0 else "fail"})

    failures = [item for item in tests if item["status"] == "fail"]
    return {"tested": len(tests), "failed": len(failures), "tests": tests, "details": details}


def test_index_controls(page: Page) -> Dict[str, Any]:
    tests: List[Dict[str, Any]] = []
    scroller = page.locator(".table-scroll")
    if scroller.count():
        scrollable = scroller.evaluate("el => el.scrollWidth > el.clientWidth + 1")
        before = scroller.evaluate("el => el.scrollLeft")
        scroller.focus()
        page.keyboard.press("ArrowRight")
        page.wait_for_timeout(100)
        after = scroller.evaluate("el => el.scrollLeft")
        scroller.evaluate("el => { el.scrollLeft = 0; }")
        passed = (not scrollable) or after > before
        tests.append({"control": "cohort table horizontal keyboard scroll", "scrollable": scrollable, "before": before, "after": after, "status": "pass" if passed else "fail"})
    else:
        tests.append({"control": "cohort table horizontal keyboard scroll", "status": "fail", "reason": "missing .table-scroll"})
    sample_links = page.locator('a[href$=".html"]')
    tests.append({"control": "sample navigation links", "count": sample_links.count(), "status": "pass" if sample_links.count() >= 7 else "fail"})
    failures = [item for item in tests if item["status"] == "fail"]
    return {"tested": len(tests), "failed": len(failures), "tests": tests}


def attach_runtime_handlers(page: Page, runtime: Dict[str, Any]) -> None:
    page.on(
        "console",
        lambda message: runtime["console_errors"].append({"text": message.text, "location": message.location}) if message.type == "error" else None,
    )
    page.on("pageerror", lambda error: runtime["page_errors"].append(str(error)))
    page.on(
        "requestfailed",
        lambda request: runtime["request_failures"].append({"url": request.url, "failure": request.failure}),
    )
    page.on(
        "response",
        lambda response: runtime["response_errors"].append({"url": response.url, "status": response.status}) if response.status >= 400 else None,
    )


def audit_page_viewport(
    browser: Browser,
    source: Path,
    source_dir: Path,
    output_dir: Path,
    viewport_name: str,
    timeout_ms: int,
    reuse_screenshots: bool,
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
        "capture_order": [
            "01_load_and_wait",
            "02_capture_initial_full_page_before_interpretation",
            "03_select_tree_candidate_for_sample",
            "04_capture_region_and_tree_evidence",
            "05_collect_metrics_and_test_controls",
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
    page = context.new_page()
    page.set_default_timeout(timeout_ms)
    attach_runtime_handlers(page, result["runtime"])
    started = time.monotonic()

    try:
        nav_started = time.monotonic()
        page.goto(source.as_uri(), wait_until="load", timeout=timeout_ms)
        page.wait_for_load_state("networkidle", timeout=timeout_ms)
        if page_kind == "index":
            page.wait_for_function("document.querySelectorAll('#cohort-table tbody tr').length === 7", timeout=timeout_ms)
        else:
            page.wait_for_function("document.querySelectorAll('#list .row').length > 0", timeout=timeout_ms)
        page.wait_for_timeout(350)
        page.evaluate("scrollTo(0, 0)")
        result["load"] = {
            "status": "pass",
            "milliseconds": round((time.monotonic() - nav_started) * 1000, 1),
            "ready_state": page.evaluate("document.readyState"),
            "url": page.url,
        }

        full_key = (page_name, viewport_name, "full_page")
        full_plan = SCREENSHOT_PLAN[full_key]
        full_path = output_dir / SCREENSHOT_DIR_NAME / full_plan["filename"]
        full_meta = capture_page_screenshot(page, full_path, full_page=True, reuse=reuse_screenshots)
        full_record = register_screenshot(screenshots, full_key, full_meta)
        result["screenshots"].append(full_record)

        candidate: Optional[Dict[str, Any]] = None
        if page_kind == "sample":
            candidate = choose_candidate_region(page)
            result["candidate"] = candidate
            if viewport_name == "desktop":
                page.locator("#regions").scroll_into_view_if_needed(timeout=30_000)
                page.evaluate("scrollBy(0, -8)")
            else:
                page.locator("#detail").scroll_into_view_if_needed(timeout=30_000)
                page.evaluate("scrollBy(0, -8)")
            page.wait_for_timeout(120)
            region_key = (page_name, viewport_name, "region_view")
            region_plan = SCREENSHOT_PLAN[region_key]
            region_path = output_dir / SCREENSHOT_DIR_NAME / region_plan["filename"]
            region_meta = capture_page_screenshot(page, region_path, full_page=False, reuse=reuse_screenshots)
            region_record = register_screenshot(screenshots, region_key, region_meta)
            result["screenshots"].append(region_record)

            if page_name == SPECIAL_TREE_SAMPLE:
                tree_state = mark_largest_tree_block(page)
                result["special_tree_state"] = tree_state
                tree_key = (page_name, viewport_name, "candidate_tree")
                tree_plan = SCREENSHOT_PLAN[tree_key]
                tree_path = output_dir / SCREENSHOT_DIR_NAME / tree_plan["filename"]
                locator = page.locator('[data-audit-largest-tree="true"]')
                tree_meta = capture_locator_screenshot(locator, tree_path, reuse=reuse_screenshots)
                tree_record = register_screenshot(screenshots, tree_key, tree_meta)
                result["screenshots"].append(tree_record)

        result["dom"] = collect_dom_metrics(page, page_kind)
        result["dom"]["heading_issues"] = heading_issues(result["dom"]["headings"])
        result["link_targets"] = validate_link_targets(result["dom"], source_dir)
        result["interactions"] = {
            "skip_link": test_skip_link(page),
            "main_controls": test_index_controls(page) if page_kind == "index" else test_sample_controls(page, candidate or {}),
        }
        result["interactions"]["keyboard_focus"] = test_keyboard_focus(page)

        add_check(checks, "document_load", True, "load + networkidle + rendered UI", result["load"], page_name=page_name, viewport=viewport_name)
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
            add_check(checks, error_kind, not errors, {"count": 0}, {"count": len(errors), "items": errors}, page_name=page_name, viewport=viewport_name)

        document = result["dom"]["document"]
        add_check(checks, "body_horizontal_overflow", document["body_overflow_x_px"] <= 1, {"maximum_px": 1}, document["body_overflow_x_px"], page_name=page_name, viewport=viewport_name)
        add_check(checks, "unintended_offviewport_content", not document["unintended_overflow_offenders"], {"count": 0}, document["unintended_overflow_offenders"], page_name=page_name, viewport=viewport_name)
        add_check(checks, "heading_hierarchy", not result["dom"]["heading_issues"], "one H1 and no skipped visible levels", result["dom"]["heading_issues"], page_name=page_name, viewport=viewport_name)
        add_check(checks, "link_targets", result["link_targets"]["failed"] == 0, {"failed": 0}, {"tested": result["link_targets"]["tested"], "failed": result["link_targets"]["failed"]}, page_name=page_name, viewport=viewport_name)
        add_check(checks, "broken_images", not result["dom"]["broken_images"], {"count": 0}, result["dom"]["broken_images"], page_name=page_name, viewport=viewport_name)
        add_check(checks, "main_controls", result["interactions"]["main_controls"]["failed"] == 0, {"failed": 0}, result["interactions"]["main_controls"], page_name=page_name, viewport=viewport_name)
        add_check(checks, "skip_link_keyboard", result["interactions"]["skip_link"]["status"] == "pass", "first Tab exposes skip link and Enter reaches target", result["interactions"]["skip_link"], page_name=page_name, viewport=viewport_name)
        focus = result["interactions"]["keyboard_focus"]
        add_check(checks, "keyboard_focus_indicators", focus["missing_indicator_count"] == 0, {"missing": 0}, {"available": focus["available"], "steps": focus["steps"], "missing": focus["missing_indicator_count"], "targets": focus["missing_indicators"]}, page_name=page_name, viewport=viewport_name)
        add_check(checks, "region_scope_label", result["dom"]["scope_labels"]["region_scope_visible"], True, result["dom"]["scope_labels"], page_name=page_name, viewport=viewport_name)
        if page_kind == "sample":
            add_check(checks, "all_chromosomes_scope_control", result["dom"]["scope_labels"]["all_chromosomes_visible"], "染色體 全部 is visible", result["dom"]["scope_labels"], page_name=page_name, viewport=viewport_name)
            trees = result["dom"]["tree_svgs"]
            readable = bool(trees) and all(not tree["clipped"] and tree["estimated_label_px"] >= 8 for tree in trees)
            add_check(checks, "tree_svg_readability", readable, "rendered tree SVG, not clipped, estimated labels >= 8px", trees, page_name=page_name, viewport=viewport_name, severity="warning")
        else:
            add_check(
                checks,
                "explicit_whole_genome_scope_label",
                result["dom"]["scope_labels"]["explicit_whole_genome_visible"],
                "explicit 全基因組 / whole-genome label",
                result["dom"]["scope_labels"],
                page_name=page_name,
                viewport=viewport_name,
                severity="warning",
            )
    except Exception as exc:
        result["fatal_error"] = str(exc)
        add_check(checks, "page_viewport_audit", False, "completed without fatal exception", str(exc), page_name=page_name, viewport=viewport_name)
    finally:
        result["duration_seconds"] = round(time.monotonic() - started, 3)
        context.close()
    return result, checks


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Capture the complete 8-page layered workstation desktop/mobile before baseline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Exit codes: 0 all checks pass, 1 findings remain, 2 runner/configuration error.",
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR, help="Directory containing index.html and seven sample pages.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help="Audit output directory.")
    parser.add_argument("--metrics", type=Path, default=DEFAULT_METRICS, help="Metrics JSON output path.")
    parser.add_argument("--timeout-ms", type=int, default=90_000, help="Per-operation Playwright timeout.")
    parser.add_argument("--reuse-screenshots", action="store_true", help="Refresh metrics and interactions while reusing existing PNG evidence.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    source_dir = args.input_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    metrics_path = args.metrics.expanduser().resolve()

    if PLAYWRIGHT_IMPORT_ERROR:
        print(json.dumps({"status": "error", "error": PLAYWRIGHT_IMPORT_ERROR}, ensure_ascii=False))
        return 2
    if args.timeout_ms <= 0:
        print(json.dumps({"status": "error", "error": "--timeout-ms must be positive"}, ensure_ascii=False))
        return 2
    if output_dir == source_dir or source_dir in output_dir.parents:
        print(json.dumps({"status": "error", "error": "Output directory must not be inside the source HTML directory"}, ensure_ascii=False))
        return 2
    if metrics_path.parent != output_dir:
        print(json.dumps({"status": "error", "error": f"Metrics must be directly under {output_dir}"}, ensure_ascii=False))
        return 2

    sources = [source_dir / filename for filename in ALL_FILES]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        print(json.dumps({"status": "error", "error": "Missing required HTML", "missing": missing}, ensure_ascii=False, indent=2))
        return 2

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / SCREENSHOT_DIR_NAME).mkdir(parents=True, exist_ok=True)
    before = {path.name: file_record(path) for path in sources}
    started_at = utc_now()
    report: Dict[str, Any] = {
        "schema_version": "1.0",
        "audit_kind": "layered_workstation_complete_before_visual_audit",
        "task_type": "B_comprehensive_validation",
        "scope": {
            "pages_expected": 8,
            "pages": ALL_FILES,
            "sample_pages_expected": 7,
            "viewports": VIEWPORTS,
            "screenshots_expected": 32,
            "partial": False,
        },
        "generated_at": started_at,
        "input": {"directory": str(source_dir), "files_before": before},
        "output": {"directory": str(output_dir), "metrics": str(metrics_path), "screenshots": str(output_dir / SCREENSHOT_DIR_NAME)},
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
        browser = playwright.chromium.launch(headless=True, args=["--allow-file-access-from-files"])
        report["browser"] = {"engine": "chromium", "version": browser.version, "headless": True}
        for source in sources:
            for viewport_name in VIEWPORTS:
                page_result, page_checks = audit_page_viewport(
                    browser,
                    source,
                    source_dir,
                    output_dir,
                    viewport_name,
                    args.timeout_ms,
                    args.reuse_screenshots,
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

    after = {path.name: file_record(path) for path in sources}
    hash_changes = []
    for name in ALL_FILES:
        if before[name]["sha256"] != after[name]["sha256"]:
            hash_changes.append({"file": name, "before": before[name]["sha256"], "after": after[name]["sha256"]})
    report["input"]["files_after"] = after
    report["input"]["hash_changes"] = hash_changes
    report["input"]["unchanged"] = not hash_changes
    add_check(report["checks"], "source_html_hashes_unchanged", not hash_changes, {"changes": 0}, hash_changes)

    report["screenshots"].sort(key=lambda item: item["step"])
    screenshot_steps = [item["step"] for item in report["screenshots"]]
    add_check(
        report["checks"],
        "complete_screenshot_plan",
        screenshot_steps == list(range(1, 33)) and all(item.get("status") == "pass" for item in report["screenshots"]),
        {"steps": list(range(1, 33)), "failed": 0},
        {"steps": screenshot_steps, "failed": [item["step"] for item in report["screenshots"] if item.get("status") != "pass"]},
    )
    audited_pairs = {(item.get("page"), item.get("viewport")) for item in report["page_viewports"] if not item.get("fatal_error")}
    expected_pairs = {(name, viewport) for name in ALL_FILES for viewport in VIEWPORTS}
    add_check(report["checks"], "complete_page_viewport_scope", audited_pairs == expected_pairs, sorted(expected_pairs), sorted(audited_pairs))

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
