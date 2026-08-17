#!/usr/bin/env python3
"""Browser and static QA for the standalone P0 correction report."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlparse

from playwright.sync_api import sync_playwright


VIEWPORTS = (("desktop", 1440, 1000), ("mobile", 390, 844))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("html", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    html_path = args.html.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    source = html_path.read_text(encoding="utf-8")

    findings: dict[str, object] = {
        "input": str(html_path),
        "size_bytes": len(source.encode("utf-8")),
        "screenshots": {},
        "print_pdf": "",
        "viewports": {},
        "console_errors": [],
        "page_errors": [],
        "external_requests": [],
        "broken_internal_links": [],
        "missing_anchor_targets": [],
        "svg": {},
    }

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        for name, width, height in VIEWPORTS:
            page = browser.new_page(viewport={"width": width, "height": height})
            page.on(
                "console",
                lambda msg: findings["console_errors"].append(msg.text)
                if msg.type == "error"
                else None,
            )
            page.on("pageerror", lambda exc: findings["page_errors"].append(str(exc)))
            page.on(
                "request",
                lambda request: findings["external_requests"].append(request.url)
                if request.url.startswith(("http://", "https://"))
                else None,
            )
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
                  tocLinks: document.querySelectorAll('aside.toc a[href^="#"]').length,
                  details: document.querySelectorAll('details').length,
                  tables: document.querySelectorAll('table').length,
                  localAnchors: [...document.querySelectorAll('a[href^="#"]')]
                    .map(a => a.getAttribute('href')),
                  uncontainedWide: [...document.querySelectorAll('body *')]
                    .filter(e => {
                      const r = e.getBoundingClientRect();
                      const wide = r.right > document.documentElement.clientWidth + 1 || r.left < -1;
                      return wide && !e.closest('.table-wrap, pre');
                    })
                    .slice(0, 20)
                    .map(e => ({
                      tag: e.tagName,
                      cls: e.className?.toString().slice(0, 80) || '',
                      text: (e.textContent || '').trim().replace(/\s+/g, ' ').slice(0, 100),
                      left: Math.round(e.getBoundingClientRect().left),
                      right: Math.round(e.getBoundingClientRect().right)
                    })),
                  svgs: [...document.querySelectorAll('svg')].map(s => ({
                    viewBox: s.getAttribute('viewBox'),
                    role: s.getAttribute('role'),
                    title: s.querySelector('title')?.textContent || '',
                    desc: s.querySelector('desc')?.textContent || '',
                    labelledBy: s.getAttribute('aria-labelledby') || ''
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
                    "all_have_aria_labelledby": all(
                        item["labelledBy"] for item in metrics["svgs"]
                    ),
                }
                hrefs = page.locator(
                    'a[href]:not([href^="#"]):not([href^="http"]):not([href^="mailto:"])'
                ).evaluate_all("els => els.map(e => e.getAttribute('href'))")
                for href in hrefs:
                    parsed = urlparse(href)
                    target = (html_path.parent / unquote(parsed.path)).resolve()
                    if not target.exists():
                        findings["broken_internal_links"].append(
                            {"href": href, "resolved": str(target)}
                        )

                page.emulate_media(media="print")
                hidden_print_bodies = page.locator(
                    "details:not([open]) > .card-body"
                ).evaluate_all(
                    "els => els.filter(e => getComputedStyle(e).display === 'none').length"
                )
                findings["hidden_print_detail_bodies"] = hidden_print_bodies
                print_pdf = out_dir / "print.pdf"
                page.pdf(path=str(print_pdf), format="A4", print_background=True)
                findings["print_pdf"] = str(print_pdf)
            page.close()
        browser.close()

    static = {
        "doctype": source.lstrip().lower().startswith("<!doctype html>"),
        "lang_zh_tw": '<html lang="zh-TW">' in source,
        "viewport_meta": 'name="viewport"' in source,
        "external_asset_markup": bool(
            re.search(
                r'(?:src|href)=["\']https?://|<script[^>]+src=|<link[^>]+stylesheet',
                source,
                re.I,
            )
        ),
        "gradient_count": len(re.findall(r"(?:linear|radial)-gradient", source, re.I)),
        "backdrop_filter": "backdrop-filter" in source.lower(),
        "text_shadow": "text-shadow" in source.lower(),
        "glow_shadow": bool(re.search(r"box-shadow\s*:[^;]*0\s+0\s+\d+px", source, re.I)),
        "svg_filter": bool(re.search(r"<(?:filter|feGaussianBlur)\b", source, re.I)),
        "svg_animation": bool(re.search(r"<(?:animate|animateTransform)\b", source, re.I)),
        "placeholder": "{{" in source or "}}" in source,
        "oversize": len(source.encode("utf-8")) >= 500_000,
        "emoji_heading": bool(
            re.search(r"<h[1-6][^>]*>[^<]*(?:✅|❌|⭐|🔴|🟠|🟡|🟢|🔍|📌|🧪)", source)
        ),
    }
    findings["static"] = static

    viewports = findings["viewports"].values()
    svg = findings["svg"]
    findings["pass"] = (
        not findings["console_errors"]
        and not findings["page_errors"]
        and not findings["external_requests"]
        and not findings["broken_internal_links"]
        and not findings["missing_anchor_targets"]
        and bool(findings["print_pdf"])
        and findings.get("hidden_print_detail_bodies") == 0
        and all(v["documentScroll"] <= v["viewport"] for v in viewports)
        and all(not v["uncontainedWide"] for v in findings["viewports"].values())
        and findings["viewports"]["desktop"]["h1"] == 1
        and findings["viewports"]["desktop"]["h2"]
        == findings["viewports"]["desktop"]["tocLinks"]
        and 1 <= svg.get("count", 0) <= 5
        and all(value for key, value in svg.items() if key != "count")
        and static["doctype"]
        and static["lang_zh_tw"]
        and static["viewport_meta"]
        and not static["external_asset_markup"]
        and static["gradient_count"] <= 2
        and not any(
            static[key]
            for key in (
                "backdrop_filter",
                "text_shadow",
                "glow_shadow",
                "svg_filter",
                "svg_animation",
                "placeholder",
                "oversize",
                "emoji_heading",
            )
        )
    )

    receipt = out_dir / "qa_receipt.json"
    receipt.write_text(
        json.dumps(findings, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(findings, ensure_ascii=False, indent=2))
    return 0 if findings["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
