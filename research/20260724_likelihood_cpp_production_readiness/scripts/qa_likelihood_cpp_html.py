#!/usr/bin/env python3
"""Deterministic browser QA for the standalone likelihood/C++ explainer."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from playwright.sync_api import Browser, sync_playwright


REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
WORK = REPO / "research/20260724_likelihood_cpp_production_readiness"
HTML = WORK / "20260724_Likelihood與C++後續流程一步步教學_01.standalone.html"
OUTPUT = WORK / "results/20260724_likelihood_cpp教學HTML_QA_01"
RECEIPT = OUTPUT / "qa_receipt.json"
CHROMIUM = Path(
    "/bip7_disk/liaoyoyo2001/.cache/ms-playwright/"
    "chromium-1223/chrome-linux64/chrome"
)

VIEWPORTS = (
    ("desktop", 1440, 1000),
    ("tablet", 1024, 900),
    ("mobile", 390, 844),
    ("narrow", 320, 700),
)

EXPECTED_SECTIONS = (
    "overview",
    "encoding",
    "emission",
    "mixture",
    "ranking",
    "edge",
    "cap",
    "hardcase",
    "july",
    "cpp",
    "decision",
    "sources",
)

REQUIRED_TEXT = (
    "122,281,152",
    "0.000209%",
    "max_sets=256",
    "candidate_generation_complete=true",
    "MAX_SETS_REACHED / ABSTAIN",
    "Production directed topology",
    "0/7",
    "96.154%",
    "July raw-all",
    "cap單獨不保證",
    "ISOLATED_PERFECT_EVENT_H_EQ_M_MINUS_D_FAST_PATH_NOT_PRODUCTION",
    "validation_evidence_eligible=false",
    "big7↔big8",
    "Gate A · content equivalence",
    "Gate B · production authority",
    "Gate C · adapter/software readiness",
    "117,760 edge rows",
    "39,846 components",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def collect_viewport(browser: Browser, name: str, width: int, height: int) -> dict:
    console_errors: list[str] = []
    page_errors: list[str] = []
    external_requests: list[str] = []
    context = browser.new_context(
        viewport={"width": width, "height": height},
        color_scheme="light",
        reduced_motion="reduce",
    )
    page = context.new_page()
    page.on(
        "console",
        lambda message: console_errors.append(message.text)
        if message.type == "error"
        else None,
    )
    page.on("pageerror", lambda error: page_errors.append(str(error)))

    def record_request(request) -> None:
        if request.url.startswith(("http://", "https://", "ws://", "wss://")):
            external_requests.append(request.url)

    page.on("request", record_request)
    page.goto(HTML.as_uri(), wait_until="load")
    page.wait_for_timeout(150)

    checks = page.evaluate(
        """({expectedSections, requiredText}) => {
          const bodyText = document.body.textContent;
          const tables = [...document.querySelectorAll('table')];
          const svgs = [...document.querySelectorAll('svg')];
          const links = [...document.querySelectorAll('.toc a')];
          const unresolvedLinks = links
            .map(link => link.getAttribute('href'))
            .filter(href => !href || !document.querySelector(href));
          return {
            lang: document.documentElement.lang,
            title: document.title,
            h1Count: document.querySelectorAll('h1').length,
            h2Count: document.querySelectorAll('h2').length,
            sectionIds: expectedSections.map(id => ({
              id,
              present: Boolean(document.getElementById(id)),
            })),
            missingText: requiredText.filter(text => !bodyText.includes(text)),
            tableCount: tables.length,
            tablesWithoutCaption: tables.filter(table => !table.querySelector('caption')).length,
            headersWithoutScope: [...document.querySelectorAll('th')]
              .filter(th => !th.hasAttribute('scope')).length,
            svgCount: svgs.length,
            svgsWithoutTitle: svgs.filter(svg => !svg.querySelector('title')).length,
            svgsWithoutDesc: svgs.filter(svg => !svg.querySelector('desc')).length,
            formulaCount: document.querySelectorAll('[role="math"]').length,
            detailsCount: document.querySelectorAll('details').length,
            unresolvedLinks,
            externalAssets: document.querySelectorAll(
              'script[src], link[rel="stylesheet"], img[src], iframe[src]'
            ).length,
            mainTextLength: (document.querySelector('main')?.innerText || '').length,
            clientWidth: document.documentElement.clientWidth,
            scrollWidth: document.documentElement.scrollWidth,
            horizontalOverflow:
              document.documentElement.scrollWidth > document.documentElement.clientWidth + 1,
          };
        }""",
        {"expectedSections": EXPECTED_SECTIONS, "requiredText": REQUIRED_TEXT},
    )

    interactions = {"run": name == "desktop"}
    if name == "desktop":
        page.locator("#bqRange").fill("10")
        interactions["bq10"] = {
            "match": page.locator("#matchProb").inner_text(),
            "flip": page.locator("#flipProb").inner_text(),
        }
        page.locator("#candidateT2").click()
        interactions["candidateT2"] = {
            "states": page.locator("#candidateStates").inner_text(),
            "ll": page.locator("#candidateLL").inner_text(),
        }
        page.locator("#capRange").fill("8192")
        interactions["cap8192"] = {
            "count": page.locator("#capCount").inner_text(),
            "percent": page.locator("#capPercent").inner_text(),
        }
        details_count = page.locator("details").count()
        page.locator("#expandAll").click()
        interactions["expandedDetails"] = page.locator("details[open]").count()
        page.locator("#collapseAll").click()
        interactions["collapsedDetails"] = details_count - page.locator("details[open]").count()
        before_theme = page.locator("html").get_attribute("data-theme")
        page.locator("#themeToggle").click()
        after_theme = page.locator("html").get_attribute("data-theme")
        interactions["theme"] = {"before": before_theme, "after": after_theme}
        page.locator("#themeToggle").click()

    screenshot = OUTPUT / f"{name}_{width}x{height}.png"
    page.evaluate("window.scrollTo(0, 0)")
    page.wait_for_timeout(100)
    page.screenshot(path=str(screenshot), full_page=False)
    context.close()

    return {
        "name": name,
        "viewport": {"width": width, "height": height},
        "checks": checks,
        "interactions": interactions,
        "console_errors": console_errors,
        "page_errors": page_errors,
        "external_requests": external_requests,
        "screenshot": str(screenshot),
    }


def collect_no_js(browser: Browser) -> dict:
    context = browser.new_context(
        java_script_enabled=False,
        viewport={"width": 390, "height": 844},
    )
    page = context.new_page()
    page.goto(HTML.as_uri(), wait_until="load")
    result = page.evaluate(
        """() => ({
          title: document.title,
          mainTextLength: (document.querySelector('main')?.innerText || '').length,
          keyConclusionPresent: document.body.innerText.includes(
            'cap單獨不保證'
          ),
          horizontalOverflow:
            document.documentElement.scrollWidth > document.documentElement.clientWidth + 1,
        })"""
    )
    context.close()
    return result


def collect_print(browser: Browser) -> dict:
    context = browser.new_context(viewport={"width": 1240, "height": 900})
    page = context.new_page()
    page.goto(HTML.as_uri(), wait_until="load")
    page.emulate_media(media="print")
    page.evaluate("document.querySelectorAll('details').forEach(node => { node.open = true; })")
    result = page.evaluate(
        """() => ({
          navDisplay: getComputedStyle(document.querySelector('.toc')).display,
          scopeBarDisplay: getComputedStyle(document.querySelector('.scope-bar')).display,
          openDetails: document.querySelectorAll('details[open]').length,
          totalDetails: document.querySelectorAll('details').length,
          bodyBackground: getComputedStyle(document.body).backgroundColor,
        })"""
    )
    pdf = OUTPUT / "print_A4.pdf"
    page.pdf(path=str(pdf), format="A4", print_background=True)
    result["pdf"] = str(pdf)
    context.close()
    return result


def result_is_valid(result: dict) -> bool:
    checks = result["checks"]
    expected_section_count = len(EXPECTED_SECTIONS)
    interactions = result["interactions"]
    structural = (
        checks["lang"] == "zh-Hant"
        and checks["h1Count"] == 1
        and checks["h2Count"] >= expected_section_count
        and all(item["present"] for item in checks["sectionIds"])
        and not checks["missingText"]
        and checks["tableCount"] >= 3
        and checks["tablesWithoutCaption"] == 0
        and checks["headersWithoutScope"] == 0
        and checks["svgCount"] >= 2
        and checks["svgsWithoutTitle"] == 0
        and checks["svgsWithoutDesc"] == 0
        and checks["formulaCount"] >= 2
        and checks["detailsCount"] >= 4
        and not checks["unresolvedLinks"]
        and checks["externalAssets"] == 0
        and checks["mainTextLength"] > 8_000
        and not checks["horizontalOverflow"]
        and not result["console_errors"]
        and not result["page_errors"]
        and not result["external_requests"]
    )
    if not interactions["run"]:
        return structural
    return structural and (
        interactions["bq10"] == {"match": "0.964286", "flip": "0.035714"}
        and interactions["candidateT2"] == {
            "states": "{RR, RA, AA}",
            "ll": "−427.745135",
        }
        and interactions["cap8192"] == {"count": "8,192", "percent": "0.006699%"}
        and interactions["expandedDetails"] == checks["detailsCount"]
        and interactions["collapsedDetails"] == checks["detailsCount"]
        and interactions["theme"] == {"before": "paper", "after": "night"}
    )


def main() -> int:
    if not HTML.is_file():
        raise FileNotFoundError(HTML)
    if not CHROMIUM.is_file():
        raise FileNotFoundError(CHROMIUM)
    OUTPUT.mkdir(parents=True, exist_ok=True)

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            headless=True,
            executable_path=str(CHROMIUM),
            args=["--no-sandbox", "--disable-dev-shm-usage", "--disable-gpu"],
        )
        try:
            viewports = [
                collect_viewport(browser, name, width, height)
                for name, width, height in VIEWPORTS
            ]
            no_js = collect_no_js(browser)
            print_result = collect_print(browser)
        finally:
            browser.close()

    passed = (
        all(result_is_valid(result) for result in viewports)
        and no_js["mainTextLength"] > 8_000
        and no_js["keyConclusionPresent"]
        and not no_js["horizontalOverflow"]
        and print_result["navDisplay"] == "none"
        and print_result["scopeBarDisplay"] == "none"
        and print_result["openDetails"] == print_result["totalDetails"]
    )
    receipt = {
        "schema_version": "1.0.0",
        "created_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "task_type": "B_COMPREHENSIVE_VALIDATION",
        "html": str(HTML),
        "html_sha256": sha256(HTML),
        "browser": str(CHROMIUM),
        "viewports": viewports,
        "no_javascript": no_js,
        "print": print_result,
        "pass": passed,
    }
    RECEIPT.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
