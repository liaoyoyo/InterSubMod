#!/usr/bin/env python3
"""Report elements that escape the portable dashboard viewport."""

import argparse
import json
from pathlib import Path

from playwright.sync_api import sync_playwright


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("html", type=Path)
    parser.add_argument("--stable-scrollbar", action="store_true")
    parser.add_argument("--browser", type=Path)
    parser.add_argument("--width", type=int, default=1440)
    parser.add_argument("--height", type=int, default=1000)
    parser.add_argument("--topbar-fix", action="store_true")
    parser.add_argument("--open-sample-filter", action="store_true")
    parser.add_argument("--select-sample")
    args = parser.parse_args()

    with sync_playwright() as playwright:
        launch_kwargs = {"headless": True}
        if args.browser:
            launch_kwargs["executable_path"] = str(args.browser.resolve())
        browser = playwright.chromium.launch(**launch_kwargs)
        page = browser.new_page(viewport={"width": args.width, "height": args.height})
        page.goto(args.html.resolve().as_uri())
        page.wait_for_selector("#data-analytics-portable-reader-root", timeout=8000)
        if args.stable_scrollbar:
            page.add_style_tag(
                content="html{overflow-y:scroll!important;scrollbar-gutter:stable!important}"
            )
        if args.topbar_fix:
            page.add_style_tag(
                content="""
                .analytics-top-bar {
                  width: calc(100% + var(--ds-gutter) + var(--ds-gutter)) !important;
                  margin-right: calc(0px - var(--ds-gutter)) !important;
                  margin-left: calc(0px - var(--ds-gutter)) !important;
                }
                @media (max-width: 600px) {
                  .chart-legend {
                    width: 100% !important;
                    max-width: 100% !important;
                    min-width: 0 !important;
                    flex: 1 1 100% !important;
                  }
                }
                """
            )
        page.wait_for_timeout(600)
        if args.open_sample_filter or args.select_sample:
            page.locator("button", has_text="Sample filter").first.click(
                force=True, timeout=5000
            )
            page.wait_for_timeout(250)
        if args.select_sample:
            option = page.locator("button").filter(has_text=args.select_sample)
            exact_option = None
            for index in range(option.count()):
                candidate = option.nth(index)
                if candidate.inner_text().strip() == args.select_sample:
                    exact_option = candidate
                    break
            if exact_option is None:
                raise RuntimeError(f"Sample option not found: {args.select_sample}")
            exact_option.click(force=True, timeout=5000)
            page.wait_for_timeout(450)
        result = page.evaluate(
            """
            () => {
              const de = document.documentElement;
              const elements = [...document.querySelectorAll('*')]
                .map((element) => {
                  const rect = element.getBoundingClientRect();
                  return {
                    tag: element.tagName,
                    className: String(element.className).slice(0, 180),
                    id: element.id,
                    left: rect.left,
                    right: rect.right,
                    width: rect.width,
                    scrollWidth: element.scrollWidth,
                    clientWidth: element.clientWidth,
                    overflowX: getComputedStyle(element).overflowX,
                  };
                })
                .filter((item) => item.right > de.clientWidth + 1 || item.left < -1)
                .sort((a, b) => b.right - a.right)
                .slice(0, 60);
              const unclippedElements = [...document.querySelectorAll('*')]
                .filter((element) => !element.closest('.table-wrap'))
                .map((element) => {
                  const rect = element.getBoundingClientRect();
                  return {
                    tag: element.tagName,
                    className: String(element.className).slice(0, 180),
                    id: element.id,
                    left: rect.left,
                    right: rect.right,
                    width: rect.width,
                    scrollWidth: element.scrollWidth,
                    clientWidth: element.clientWidth,
                    overflowX: getComputedStyle(element).overflowX,
                  };
                })
                .filter((item) => item.right > de.clientWidth + 1 || item.left < -1)
                .sort((a, b) => b.right - a.right)
                .slice(0, 60);
              const widest = [...document.querySelectorAll('.table-scroll-content')]
                .sort((a, b) => b.getBoundingClientRect().right - a.getBoundingClientRect().right)[0];
              const ancestry = [];
              for (let node = widest; node; node = node.parentElement) {
                const rect = node.getBoundingClientRect();
                ancestry.push({
                  tag: node.tagName,
                  className: String(node.className).slice(0, 180),
                  id: node.id,
                  left: rect.left,
                  right: rect.right,
                  width: rect.width,
                  scrollWidth: node.scrollWidth,
                  clientWidth: node.clientWidth,
                  overflowX: getComputedStyle(node).overflowX,
                  minWidth: getComputedStyle(node).minWidth,
                  maxWidth: getComputedStyle(node).maxWidth,
                });
              }
              const legend = [...document.querySelectorAll('.chart-legend')]
                .sort((a, b) => b.getBoundingClientRect().width - a.getBoundingClientRect().width)[0];
              const legendAncestry = [];
              for (let node = legend; node; node = node.parentElement) {
                const rect = node.getBoundingClientRect();
                const style = getComputedStyle(node);
                legendAncestry.push({
                  tag: node.tagName,
                  className: String(node.className).slice(0, 180),
                  id: node.id,
                  left: rect.left,
                  right: rect.right,
                  width: rect.width,
                  display: style.display,
                  position: style.position,
                  flexWrap: style.flexWrap,
                  justifyContent: style.justifyContent,
                  overflowX: style.overflowX,
                  widthCss: style.width,
                  minWidth: style.minWidth,
                  maxWidth: style.maxWidth,
                  gap: style.gap,
                  padding: style.padding,
                  margin: style.margin,
                });
              }
              const dashboardRules = [];
              for (const sheet of [...document.styleSheets]) {
                try {
                  for (const rule of [...sheet.cssRules]) {
                    if (
                      rule.cssText?.includes('.dashboard-shell') ||
                      rule.cssText?.includes('.analytics-top-bar') ||
                      rule.cssText?.includes('.chart-legend')
                    ) {
                      dashboardRules.push(rule.cssText.slice(0, 1000));
                    }
                  }
                } catch (_error) {
                  // Ignore inaccessible stylesheets; portable HTML should not have any.
                }
              }
              const controls = [...document.querySelectorAll('select, button')].map((node) => ({
                tag: node.tagName,
                ariaLabel: node.getAttribute('aria-label'),
                name: node.getAttribute('name'),
                text: node.textContent?.trim().slice(0, 160),
                value: node.value,
                options: node.tagName === 'SELECT'
                  ? [...node.options].map((option) => ({ text: option.textContent, value: option.value }))
                  : undefined,
              }));
              const metricCards = [
                ...document.querySelectorAll('.metric-card, .report-metric-card'),
              ].map((node) => ({
                className: String(node.className),
                text: node.textContent?.trim().replace(/\s+/g, ' ').slice(0, 500),
              }));
              const renderedTables = [...document.querySelectorAll('section.table-panel')].map((node) => ({
                id: node.id,
                rows: node.querySelectorAll('tbody tr').length,
                text: node.textContent?.trim().replace(/\s+/g, ' ').slice(0, 500),
              }));
              return {
                innerWidth,
                clientWidth: de.clientWidth,
                scrollWidth: de.scrollWidth,
                bodyScrollWidth: document.body.scrollWidth,
                elements,
                unclippedElements,
                widestTableAncestry: ancestry,
                legendAncestry,
                dashboardRules,
                controls,
                metricCards,
                renderedTables,
              };
            }
            """
        )
        print(json.dumps(result, indent=2))
        browser.close()


if __name__ == "__main__":
    main()
