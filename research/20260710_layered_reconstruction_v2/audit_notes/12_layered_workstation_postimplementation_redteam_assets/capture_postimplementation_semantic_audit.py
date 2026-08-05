#!/usr/bin/env python3
"""Current-run Chromium/Playwright semantic and responsive audit."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from playwright.sync_api import sync_playwright


PAGES = [
    "index",
    "HCC1395",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
]

TARGETS = {
    "HCC1395": {"ftopo": "exact_and_topology_unique", "fevidence": "partial_only"},
    "COLO829": {"ftopo": "topology_unique_exact_multiple"},
    "H1437": {"ftopo": "incomplete"},
    "H2009": {
        "ftopo": "incomplete",
        "fsignal": "recurrence",
        "region": "chr8:79992384-80149786",
        "expected_filtered_before_region": 4,
    },
    "HCC1395_DORADO": {
        "ftopo": "no_primary_lineage",
        "region": "chr1:190064024-190196077",
        "expected_filtered_before_region": 720,
    },
    "HCC1937": {"ftopo": "exact_and_topology_unique"},
    "HCC1954": {"ftopo": "topology_multiple_exact_multiple"},
}

VIEWPORTS = {
    "desktop": {"width": 1440, "height": 1000},
    "mobile": {"width": 390, "height": 844},
}


def integer_status(text: str) -> int:
    return int(text.split()[0].replace(",", ""))


def visible_term_counts(page):
    text = page.locator("body").inner_text().lower()
    return {
        "softmax": text.count("softmax"),
        "posterior": text.count("posterior"),
        "vaf": len(re.findall(r"\bvaf\b", text)),
        "clone-c": len(re.findall(r"\bclone[- ]c\b", text)),
        "read_alt_fraction": text.count("read alt fraction"),
        "json": text.count("json"),
    }


def audit_index(page, viewport_name, shot_path):
    raw = page.locator("details.evidence-drawer")
    metrics = page.locator(".metric-card .metric-number").all_inner_texts()
    output = {
        "title": page.title(),
        "h1": page.locator("h1").inner_text(),
        "dataset_launchers": page.locator(".genome-launchers .genome-link").count(),
        "cohort_rows": page.locator("#cohort-table tbody tr").count(),
        "sample_topology_cards": page.locator(".sample-topology").count(),
        "metric_numbers": metrics,
        "topology_segments": page.locator(".topology-main .bar-segment").count(),
        "raw_drawer_open": raw.get_attribute("open") is not None,
        "raw_json_links_visible": sum(link.is_visible() for link in raw.locator('a[href$=".json"]').all()),
        "visible_terms": visible_term_counts(page),
    }
    if viewport_name == "desktop":
        page.screenshot(path=str(shot_path), full_page=True)
    else:
        page.screenshot(path=str(shot_path), full_page=True)
    return output


def audit_sample(page, sample, viewport_name, shot_path):
    target = TARGETS[sample]
    initial = {
        "chromosome_buttons": page.locator(".chrom-button").count(),
        "result_count": integer_status(page.locator("#result-count").inner_text()),
        "raw_drawer_open": page.locator("details.evidence-drawer").get_attribute("open") is not None,
        "raw_json_links_visible": sum(link.is_visible() for link in page.locator("details.evidence-drawer a").all()),
        "metrics": page.locator("#sample-metrics .metric .value").all_inner_texts(),
        "topology_legend_cards": page.locator("#topology-legend .legend-card").count(),
        "filters": page.locator(".filters select, .filters input").count(),
    }

    for control in ("ftopo", "fevidence", "fsignal"):
        if target.get(control):
            page.locator(f"#{control}").select_option(target[control])
    page.wait_for_timeout(80)
    filtered_before_region = integer_status(page.locator("#result-count").inner_text())
    if target.get("expected_filtered_before_region") is not None:
        assert filtered_before_region == target["expected_filtered_before_region"], (
            sample,
            filtered_before_region,
            target["expected_filtered_before_region"],
        )
    if target.get("region"):
        page.locator("#fq").fill(target["region"])
        page.wait_for_timeout(80)
    assert integer_status(page.locator("#result-count").inner_text()) >= 1
    first_result = page.locator(".result-row").first
    first_result.focus()
    focus_before_region_activation = page.evaluate(
        "({tag:document.activeElement.tagName,id:document.activeElement.id,cls:document.activeElement.className})"
    )
    first_result.press("Enter")
    page.locator("#detail-title").wait_for(state="visible")
    focus_after_region_activation = page.evaluate(
        "({tag:document.activeElement.tagName,id:document.activeElement.id,cls:document.activeElement.className})"
    )
    selected_region = page.locator("#detail-title").inner_text()
    if target.get("region"):
        assert selected_region == target["region"]

    page.locator("#detail").scroll_into_view_if_needed()
    page.wait_for_timeout(120)
    svg_titles = page.locator("#detail .network-scroll svg > title").evaluate_all(
        "els => els.map(el => el.textContent || '')"
    )
    edge_titles = page.locator("#detail .network-scroll svg line > title").evaluate_all(
        "els => els.map(el => el.textContent || '')"
    )
    network_aria = [
        page.locator("#detail .network-scroll svg").nth(index).aria_snapshot()
        for index in range(min(4, page.locator("#detail .network-scroll svg").count()))
    ]
    shape_tabs = page.locator("#detail .shape-tab")
    tree_tabs = page.locator("#detail .tree-tab")
    selected_shape = page.locator('#detail .shape-tab[aria-selected="true"]').first
    tab_keyboard = {"tested": False}
    if shape_tabs.count():
        selected_shape.focus()
        before_focus = page.evaluate("document.activeElement && document.activeElement.textContent")
        before_selected = selected_shape.inner_text()
        page.keyboard.press("ArrowRight")
        after_focus = page.evaluate("document.activeElement && document.activeElement.textContent")
        after_selected = page.locator('#detail .shape-tab[aria-selected="true"]').first.inner_text()
        tab_keyboard = {
            "tested": True,
            "before_focus": before_focus,
            "after_focus": after_focus,
            "before_selected": before_selected,
            "after_selected": after_selected,
            "arrow_right_changed_focus_or_selection": before_focus != after_focus or before_selected != after_selected,
        }

    tree_focus = {"tested": False}
    if tree_tabs.count() > 1:
        tree_tabs.nth(1).focus()
        before = page.evaluate(
            "({tag:document.activeElement.tagName,id:document.activeElement.id,cls:document.activeElement.className,text:(document.activeElement.innerText||'').slice(0,80)})"
        )
        tree_tabs.nth(1).evaluate("el => el.click()")
        page.wait_for_timeout(40)
        after = page.evaluate(
            "({tag:document.activeElement.tagName,id:document.activeElement.id,cls:document.activeElement.className,text:(document.activeElement.innerText||'').slice(0,80)})"
        )
        tree_focus = {"tested": True, "before": before, "after": after, "focus_preserved": before == after}

    caption = page.locator("#detail .site-table caption").inner_text() if page.locator("#detail .site-table caption").count() else None
    detail_text = page.locator("#detail").inner_text()
    output = {
        "initial": initial,
        "filtered_before_region_query": filtered_before_region,
        "filtered_after_region_query": integer_status(page.locator("#result-count").inner_text()),
        "selected_region": selected_region,
        "focus_before_region_activation": focus_before_region_activation,
        "focus_after_region_activation": focus_after_region_activation,
        "region_activation_focuses_detail": focus_after_region_activation.get("id") == "detail",
        "url_hash": page.evaluate("location.hash"),
        "verdict_text": page.locator("#detail .verdict p").nth(1).inner_text(),
        "verdict_badges": page.locator("#detail .verdict .badge").all_inner_texts(),
        "C_Topo": page.locator("#detail .ct-box strong").all_inner_texts(),
        "primary_lane_count": page.locator("#detail .lane:not(.auxiliary)").count(),
        "auxiliary_lane_count": page.locator("#detail .lane.auxiliary").count(),
        "candidate_complete_badges": page.locator("#detail .lane:not(.auxiliary) .lane-title .badge").all_inner_texts(),
        "scope_warning_count": page.locator("#detail .scope-warning").count(),
        "network_svg_count": page.locator("#detail .network-scroll svg").count(),
        "network_svg_role_img_count": page.locator('#detail .network-scroll svg[role="img"][aria-label]').count(),
        "network_region_count": page.locator('#detail .network-scroll[role="region"][tabindex="0"][aria-label]').count(),
        "network_svg_titles": svg_titles,
        "network_aria_snapshots": network_aria,
        "network_root_edge_count": sum(" edge ROOT>" in (title or "") for title in edge_titles),
        "network_acquisition_texts": page.locator("#detail .network-scroll svg text").evaluate_all(
            "els => els.map(el => el.textContent || '').filter(text => text.trim().startsWith('+'))"
        ),
        "edge_title_status_counts": {
            "forced": sum((title or "").startswith("forced edge ") for title in edge_titles),
            "variable": sum((title or "").startswith("variable edge ") for title in edge_titles),
            "unevaluated": sum((title or "").startswith("unevaluated edge ") for title in edge_titles),
        },
        "shape_tabs": shape_tabs.count(),
        "tree_tabs": tree_tabs.count(),
        "tab_aria_controls_count": page.locator('#detail [role="tab"][aria-controls]').count(),
        "tabpanel_count": page.locator('#detail [role="tabpanel"]').count(),
        "tab_explicit_tabindex_count": page.locator('#detail [role="tab"][tabindex]').count(),
        "tab_keyboard": tab_keyboard,
        "tree_candidate_activation_focus": tree_focus,
        "site_table_caption": caption,
        "site_table_text": page.locator("#detail .site-table").inner_text() if page.locator("#detail .site-table").count() else None,
        "detail_has_recurrence_facet": "recurrence facet" in detail_text,
        "detail_has_candidate_incomplete": "candidate incomplete" in detail_text,
        "detail_has_edge_scope_incomplete": "edge scope incomplete" in detail_text,
        "detail_visible_forced_legend_count": detail_text.count("forced across complete exact set"),
        "detail_visible_no_primary_warning": "此區沒有 primary mutation-bearing HP1／HP2" in detail_text,
        "raw_json_links_visible_after_selection": sum(
            link.is_visible() for link in page.locator("details.evidence-drawer a").all()
        ),
        "visible_terms": visible_term_counts(page),
    }
    page.screenshot(path=str(shot_path), full_page=False)
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-url", default="http://127.0.0.1:8765/layered_workstation")
    parser.add_argument(
        "--outdir",
        default=str(Path(__file__).resolve().parent),
    )
    args = parser.parse_args()
    outdir = Path(args.outdir).resolve()
    screenshot_dir = outdir / "screenshots"
    screenshot_dir.mkdir(parents=True, exist_ok=True)
    receipt = {
        "audit": "postimplementation_current_run_chromium_semantic_audit",
        "base_url": args.base_url,
        "browser": {},
        "viewports": VIEWPORTS,
        "pages": {},
        "console_errors": [],
        "page_errors": [],
        "request_failures": [],
    }
    executable = "/bip7_disk/liaoyoyo2001/.cache/ms-playwright/chromium-1217/chrome-linux64/chrome"
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True, executable_path=executable)
        receipt["browser"] = {
            "engine": "Chromium",
            "version": browser.version,
            "executable": executable,
        }
        for viewport_name, viewport in VIEWPORTS.items():
            context = browser.new_context(viewport=viewport, device_scale_factor=1)
            for name in PAGES:
                page = context.new_page()
                local_console = []
                local_page_errors = []
                local_request_failures = []
                page.on("console", lambda msg, bucket=local_console: bucket.append(msg.text) if msg.type == "error" else None)
                page.on("pageerror", lambda error, bucket=local_page_errors: bucket.append(str(error)))
                page.on("requestfailed", lambda request, bucket=local_request_failures: bucket.append({"url": request.url, "failure": request.failure}))
                url = f"{args.base_url}/{name}.html"
                page.goto(url, wait_until="domcontentloaded", timeout=120_000)
                if name == "index":
                    page.locator(".genome-launchers").wait_for(state="visible")
                else:
                    page.locator("#result-count").wait_for(state="visible")
                shot = screenshot_dir / f"{name}_{viewport_name}.png"
                item = audit_index(page, viewport_name, shot) if name == "index" else audit_sample(page, name, viewport_name, shot)
                overflow = page.evaluate(
                    "({innerWidth, bodyScrollWidth:document.body.scrollWidth, docScrollWidth:document.documentElement.scrollWidth, bodyScrollHeight:document.body.scrollHeight})"
                )
                item.update({
                    "url": page.url,
                    "http_status": 200,
                    "overflow": overflow,
                    "horizontal_page_overflow": max(overflow["bodyScrollWidth"], overflow["docScrollWidth"]) > overflow["innerWidth"] + 1,
                    "console_errors": local_console,
                    "page_errors": local_page_errors,
                    "request_failures": local_request_failures,
                    "screenshot": str(shot),
                    "screenshot_bytes": shot.stat().st_size,
                })
                receipt["pages"][f"{name}_{viewport_name}"] = item
                receipt["console_errors"].extend({"page": f"{name}_{viewport_name}", "message": value} for value in local_console)
                receipt["page_errors"].extend({"page": f"{name}_{viewport_name}", "message": value} for value in local_page_errors)
                receipt["request_failures"].extend({"page": f"{name}_{viewport_name}", **value} for value in local_request_failures)
                print(f"PASS capture {name} {viewport_name}: {shot.name}", flush=True)
                page.close()
            context.close()
        browser.close()
    receipt["summary"] = {
        "page_viewport_runs": len(receipt["pages"]),
        "screenshots": len(list(screenshot_dir.glob("*.png"))),
        "console_errors": len(receipt["console_errors"]),
        "page_errors": len(receipt["page_errors"]),
        "request_failures": len(receipt["request_failures"]),
        "horizontal_page_overflows": sum(item["horizontal_page_overflow"] for item in receipt["pages"].values()),
        "all_raw_drawers_initially_collapsed": all(
            not item.get("raw_drawer_open", item.get("initial", {}).get("raw_drawer_open", True))
            for item in receipt["pages"].values()
        ),
    }
    out = outdir / "playwright_semantic_receipt.json"
    out.write_text(json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"output": str(out), "summary": receipt["summary"]}, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
