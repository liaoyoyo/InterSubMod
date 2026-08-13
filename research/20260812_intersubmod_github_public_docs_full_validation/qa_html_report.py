#!/usr/bin/env python3
"""Render and mechanically QA the standalone audit report."""

import argparse
import json
from pathlib import Path
from urllib.parse import unquote, urlparse

from playwright.sync_api import sync_playwright


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("html", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    html_path = args.html.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    findings = {
        "input": str(html_path),
        "screenshots": {},
        "print_pdf": "",
        "viewports": {},
        "console_errors": [],
        "page_errors": [],
        "broken_internal_links": [],
        "missing_anchor_targets": [],
        "svg": {},
    }

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        for name, width, height in (("desktop", 1440, 1000), ("mobile", 390, 844)):
            page = browser.new_page(viewport={"width": width, "height": height})
            page.on("console", lambda msg: findings["console_errors"].append(msg.text) if msg.type == "error" else None)
            page.on("pageerror", lambda exc: findings["page_errors"].append(str(exc)))
            page.goto(html_path.as_uri(), wait_until="networkidle")
            screenshot = out_dir / f"{name}.png"
            page.screenshot(path=str(screenshot), full_page=True)

            metrics = page.evaluate(
                """() => ({
                    viewport: document.documentElement.clientWidth,
                    documentScroll: document.documentElement.scrollWidth,
                    bodyScroll: document.body.scrollWidth,
                    h1: document.querySelectorAll('h1').length,
                    h2: document.querySelectorAll('h2').length,
                    tables: document.querySelectorAll('table').length,
                    details: document.querySelectorAll('details').length,
                    wideElements: [...document.querySelectorAll('body *')]
                      .filter(e => {
                        const r = e.getBoundingClientRect();
                        return r.right > document.documentElement.clientWidth + 1 || r.left < -1;
                      })
                      .slice(0, 20)
                      .map(e => ({
                        tag: e.tagName,
                        cls: e.className?.toString().slice(0, 80) || '',
                        text: (e.textContent || '').trim().replace(/\s+/g, ' ').slice(0, 100),
                        left: Math.round(e.getBoundingClientRect().left),
                        right: Math.round(e.getBoundingClientRect().right),
                        scrollWidth: e.scrollWidth
                      })),
                    uncontainedWide: [...document.querySelectorAll('body *')]
                      .filter(e => {
                        const r = e.getBoundingClientRect();
                        const wide = r.right > document.documentElement.clientWidth + 1 || r.left < -1;
                        return wide && !e.closest('.side,.table-wrap');
                      })
                      .slice(0, 20)
                      .map(e => ({
                        tag: e.tagName,
                        cls: e.className?.toString().slice(0, 80) || '',
                        text: (e.textContent || '').trim().replace(/\s+/g, ' ').slice(0, 100),
                        left: Math.round(e.getBoundingClientRect().left),
                        right: Math.round(e.getBoundingClientRect().right)
                      })),
                    localAnchors: [...document.querySelectorAll('a[href^="#"]')].map(a => a.getAttribute('href')),
                    svgs: [...document.querySelectorAll('svg')].map(s => ({
                        viewBox: s.getAttribute('viewBox'),
                        role: s.getAttribute('role'),
                        title: s.querySelector('title')?.textContent || '',
                        desc: s.querySelector('desc')?.textContent || ''
                    }))
                })"""
            )
            findings["screenshots"][name] = str(screenshot)
            findings["viewports"][name] = metrics
            for anchor in metrics["localAnchors"]:
                if anchor and not page.locator(anchor).count():
                    findings["missing_anchor_targets"].append(anchor)

            if name == "desktop":
                findings["svg"] = {
                    "count": len(metrics["svgs"]),
                    "all_have_viewbox": all(item["viewBox"] for item in metrics["svgs"]),
                    "all_have_role": all(item["role"] == "img" for item in metrics["svgs"]),
                    "all_have_title": all(item["title"] for item in metrics["svgs"]),
                    "all_have_desc": all(item["desc"] for item in metrics["svgs"]),
                }

                for href in page.locator('a[href]:not([href^="#"]):not([href^="http"]):not([href^="mailto:"])').evaluate_all(
                    "els => els.map(e => e.getAttribute('href'))"
                ):
                    parsed = urlparse(href)
                    target = (html_path.parent / unquote(parsed.path)).resolve()
                    if not target.exists():
                        findings["broken_internal_links"].append({"href": href, "resolved": str(target)})
                page.emulate_media(media="print")
                print_pdf = out_dir / "print.pdf"
                page.pdf(path=str(print_pdf), format="A4", print_background=True)
                findings["print_pdf"] = str(print_pdf)
            page.close()
        browser.close()

    prohibited = {
        "external_http_assets": any(token in html_path.read_text(encoding="utf-8") for token in ('src="http://', 'src="https://', 'href="https://fonts')),
        "gradient": "gradient" in html_path.read_text(encoding="utf-8").lower(),
        "text_shadow": "text-shadow" in html_path.read_text(encoding="utf-8").lower(),
        "backdrop_filter": "backdrop-filter" in html_path.read_text(encoding="utf-8").lower(),
    }
    findings["prohibited_patterns"] = prohibited
    findings["pass"] = (
        not findings["console_errors"]
        and not findings["page_errors"]
        and not findings["broken_internal_links"]
        and not findings["missing_anchor_targets"]
        and bool(findings["print_pdf"])
        and all(v["documentScroll"] <= v["viewport"] for v in findings["viewports"].values())
        and all(not v["uncontainedWide"] for v in findings["viewports"].values())
        and findings["svg"].get("count") == 4
        and all(value for key, value in findings["svg"].items() if key != "count")
        and not any(prohibited.values())
    )

    receipt = out_dir / "qa_receipt.json"
    receipt.write_text(json.dumps(findings, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(findings, ensure_ascii=False, indent=2))
    return 0 if findings["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
