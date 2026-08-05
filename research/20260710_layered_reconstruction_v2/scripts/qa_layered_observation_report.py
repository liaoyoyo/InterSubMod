#!/usr/bin/env python3
"""Browser QA for the standalone layered-reconstruction observation report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from playwright.sync_api import sync_playwright


def visible_count(page, selector: str) -> int:
    return page.locator(selector).evaluate_all(
        "els => els.filter(el => { const s=getComputedStyle(el); const r=el.getBoundingClientRect(); "
        "return s.display !== 'none' && s.visibility !== 'hidden' && r.width > 0 && r.height > 0; }).length"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--screenshot-dir", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    args.screenshot_dir.mkdir(parents=True, exist_ok=True)

    results = {"url": args.url, "scenarios": [], "errors": []}
    required_text = [
        "7/7",
        "read-AF",
        "not_evaluated",
        "7,928",
        "7,927",
        "25,639",
        "25,635",
        "NO-GO",
        "PARTIAL SCIENTIFIC VALIDATION",
        "S → W → HP → C → Topo",
        "C=1, Topo=1",
        "Most-likely topology：NOT EVALUATED",
    ]
    expected_chart_count = 14

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        for name, width, height, color_scheme in [
            ("desktop-light", 1440, 1000, "light"),
            ("desktop-dark", 1440, 1000, "dark"),
            ("tablet-light", 768, 1000, "light"),
            ("mobile-light", 390, 844, "light"),
            ("mobile-360", 360, 800, "light"),
            ("mobile-320", 320, 720, "light"),
        ]:
            context = browser.new_context(
                viewport={"width": width, "height": height},
                color_scheme=color_scheme,
            )
            page = context.new_page()
            console_errors = []
            page_errors = []
            page.on("console", lambda message, sink=console_errors: sink.append(message.text) if message.type == "error" else None)
            page.on("pageerror", lambda error, sink=page_errors: sink.append(str(error)))
            page.goto(args.url, wait_until="networkidle")
            page.wait_for_timeout(1500)

            body_text = page.locator("body").inner_text()
            unnamed_figures = page.locator("figure:not([aria-labelledby])").count()
            tables_without_caption = page.locator("table:not(:has(caption))").count()
            th_without_scope = page.locator("th:not([scope])").count()
            summary_labels = page.locator("summary").all_inner_texts()
            duplicate_summary_labels = sorted({label for label in summary_labels if summary_labels.count(label) > 1})
            source_title_overlap = page.locator(".source-figure").evaluate_all(
                "figs => figs.filter(fig => { const b=fig.querySelector(':scope > button.source-tooltip'); const h=fig.querySelector('.card-head h3'); if(!b||!h) return false; const x=b.getBoundingClientRect(), y=h.getBoundingClientRect(); return !(x.right<y.left||x.left>y.right||x.bottom<y.top||x.top>y.bottom); }).length"
            )
            hp_chart_text = " ".join(page.locator("[data-recharts-chart='hp-multiplicity'] [data-recharts-live] svg text").evaluate_all("els => els.map(el => el.textContent || '')"))
            scenario = {
                "name": name,
                "viewport": [width, height],
                "color_scheme": color_scheme,
                "chart_hosts": page.locator("[data-recharts-chart]").count(),
                "live_svg": page.locator("[data-recharts-live] svg").count(),
                "visible_fallbacks": visible_count(page, "[data-recharts-fallback]"),
                "figure_source_buttons": page.locator(".source-figure > button.source-tooltip").count(),
                "source_tooltips": page.locator(".source-tooltip").count(),
                "remote_scripts": page.locator("script[src^='http']").count(),
                "remote_stylesheets": page.locator("link[rel='stylesheet'][href^='http']").count(),
                "missing_required_text": [text for text in required_text if text not in body_text],
                "console_errors": console_errors,
                "page_errors": page_errors,
                "unnamed_figures": unnamed_figures,
                "tables_without_caption": tables_without_caption,
                "th_without_scope": th_without_scope,
                "duplicate_summary_labels": duplicate_summary_labels,
                "skip_link_present": page.locator("a.skip-link[href='#report-main']").count() == 1,
                "inline_source_tabindex": page.locator(".inline-source[tabindex]").count(),
                "source_trigger_title_overlap": source_title_overlap,
                "hp_chart_dataset_labels": sum(1 for label in ["COLO829", "H1437", "H2009", "HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954"] if label in hp_chart_text),
            }

            page.keyboard.press("Escape")
            page.evaluate("document.activeElement && document.activeElement.blur()")
            page.keyboard.press("Tab")
            scenario["first_tab_is_skip_link"] = page.locator(":focus").evaluate("el => el.matches('a.skip-link')")

            first_source = page.locator(".source-figure > button.source-tooltip").first
            if first_source.count():
                first_source.focus()
                page.wait_for_timeout(150)
                scenario["focused_source_expanded"] = first_source.get_attribute("aria-expanded")
                scenario["focused_source_tooltip_visible"] = visible_count(
                    page, ".source-figure > button.source-tooltip:first-of-type .source-tooltip-content"
                )
            else:
                scenario["focused_source_expanded"] = None
                scenario["focused_source_tooltip_visible"] = 0

            screenshot = args.screenshot_dir / f"{name}.png"
            page.screenshot(path=str(screenshot), full_page=True)
            scenario["screenshot"] = str(screenshot)
            results["scenarios"].append(scenario)
            context.close()

        for name, width in [("no-js-desktop", 1440), ("no-js-mobile", 390)]:
            context = browser.new_context(
                viewport={"width": width, "height": 900},
                java_script_enabled=False,
                color_scheme="light",
            )
            page = context.new_page()
            page.goto(args.url, wait_until="networkidle")
            body_text = page.locator("body").inner_text()
            scenario = {
                "name": name,
                "viewport": [width, 900],
                "javascript": False,
                "chart_hosts": page.locator("[data-recharts-chart]").count(),
                "live_svg": page.locator("[data-recharts-live] svg").count(),
                "visible_fallbacks": visible_count(page, "[data-recharts-fallback]"),
                "fallback_tables": page.locator("[data-recharts-fallback] table").count(),
                "missing_required_text": [text for text in required_text if text not in body_text],
            }
            screenshot = args.screenshot_dir / f"{name}.png"
            page.screenshot(path=str(screenshot), full_page=True)
            scenario["screenshot"] = str(screenshot)
            results["scenarios"].append(scenario)
            context.close()
        browser.close()

    for scenario in results["scenarios"]:
        if scenario["chart_hosts"] != expected_chart_count:
            results["errors"].append(f"{scenario['name']}: expected {expected_chart_count} chart hosts")
        if scenario.get("javascript", True):
            if scenario["live_svg"] != expected_chart_count:
                results["errors"].append(f"{scenario['name']}: expected {expected_chart_count} live SVG charts")
            if scenario["visible_fallbacks"] != 0:
                results["errors"].append(f"{scenario['name']}: live mode must hide fallbacks")
            if scenario["figure_source_buttons"] != expected_chart_count:
                results["errors"].append(f"{scenario['name']}: expected {expected_chart_count} chart source buttons")
            if scenario["remote_scripts"] or scenario["remote_stylesheets"]:
                results["errors"].append(f"{scenario['name']}: remote runtime dependency found")
            if scenario["console_errors"] or scenario["page_errors"]:
                results["errors"].append(f"{scenario['name']}: browser console/page error")
            if scenario["focused_source_tooltip_visible"] != 1:
                results["errors"].append(f"{scenario['name']}: keyboard source tooltip not visible")
            if scenario["unnamed_figures"] or scenario["tables_without_caption"] or scenario["th_without_scope"]:
                results["errors"].append(f"{scenario['name']}: figure/table accessibility semantics failed")
            if scenario["duplicate_summary_labels"]:
                results["errors"].append(f"{scenario['name']}: duplicate details summary labels")
            if not scenario["skip_link_present"] or not scenario["first_tab_is_skip_link"]:
                results["errors"].append(f"{scenario['name']}: skip-link keyboard contract failed")
            if scenario["inline_source_tabindex"]:
                results["errors"].append(f"{scenario['name']}: inline source numbers pollute tab order")
            if scenario["source_trigger_title_overlap"]:
                results["errors"].append(f"{scenario['name']}: source trigger overlaps chart title")
            if width <= 800 and scenario["hp_chart_dataset_labels"] != 7:
                results["errors"].append(f"{scenario['name']}: HP chart does not expose all 7 dataset labels")
        else:
            if scenario["live_svg"] != 0 or scenario["visible_fallbacks"] != expected_chart_count or scenario["fallback_tables"] != expected_chart_count:
                results["errors"].append(f"{scenario['name']}: static chart fallback contract failed")
        if scenario["missing_required_text"]:
            results["errors"].append(f"{scenario['name']}: missing required text {scenario['missing_required_text']}")

    results["status"] = "pass" if not results["errors"] else "fail"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(results, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results, ensure_ascii=False, indent=2))
    raise SystemExit(0 if results["status"] == "pass" else 1)


if __name__ == "__main__":
    main()
