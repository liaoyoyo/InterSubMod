#!/usr/bin/env python3
"""Chromium desktop/mobile visual and interaction audit for exact-PS pages."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from playwright.sync_api import Browser, Page, sync_playwright


HERE = Path(__file__).resolve().parent
TOPIC_ROOT = HERE.parent
REPO_ROOT = HERE.parents[2]
WORKSTATION = REPO_ROOT / "docs" / "methodology" / "_assets" / "layered_workstation"
DEFAULT_OUTPUT = TOPIC_ROOT / "artifacts" / "exact_ps_workstation_audit"
SAMPLES = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
VIEWPORTS = {
    "desktop": {"width": 1440, "height": 1000},
    "mobile": {"width": 390, "height": 844},
}
AUTHORITY = "20260724-exact-ps-hp-strict-read-linkage"
UI_CONTRACT = "layered-workstation-exact-ps-v1"


class AuditFailure(RuntimeError):
    """Raised on a browser contract violation."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditFailure(message)


def launch_browser(playwright: Any) -> Browser:
    try:
        return playwright.chromium.launch(headless=True)
    except Exception:
        for candidate in (
            Path("/usr/bin/chromium"),
            Path("/usr/bin/chromium-browser"),
            Path("/usr/bin/google-chrome"),
        ):
            if candidate.is_file():
                return playwright.chromium.launch(
                    headless=True,
                    executable_path=str(candidate),
                    args=["--no-sandbox"],
                )
        raise


def page_observers(page: Page) -> tuple[list[str], list[str]]:
    console_errors: list[str] = []
    external_requests: list[str] = []
    page.on(
        "console",
        lambda message: console_errors.append(message.text)
        if message.type == "error"
        else None,
    )
    page.on("pageerror", lambda error: console_errors.append(str(error)))
    page.on(
        "request",
        lambda request: external_requests.append(request.url)
        if request.url.startswith(("http://", "https://"))
        else None,
    )
    return console_errors, external_requests


def geometry(page: Page) -> dict[str, Any]:
    return page.evaluate(
        """() => ({
          viewport: {w: innerWidth, h: innerHeight},
          html: {scrollWidth: document.documentElement.scrollWidth, clientWidth: document.documentElement.clientWidth},
          body: {scrollWidth: document.body.scrollWidth, clientWidth: document.body.clientWidth},
          scroll: {x: scrollX, y: scrollY},
          title: document.title,
          authority: document.querySelector('meta[name="intersubmod-authority"]')?.content,
          contract: document.querySelector('meta[name="intersubmod-ui-contract"]')?.content
        })"""
    )


def audit_index(page: Page, output: Path, mode: str) -> dict[str, Any]:
    console_errors, external_requests = page_observers(page)
    page.goto((WORKSTATION / "index.html").as_uri(), wait_until="domcontentloaded")
    page.wait_for_selector('a[href="HCC1395.html"]', timeout=30_000)
    info = geometry(page)
    require(info["authority"] == AUTHORITY, f"index {mode}: authority mismatch")
    require(info["contract"] == UI_CONTRACT, f"index {mode}: UI contract mismatch")
    require(
        info["html"]["scrollWidth"] <= info["html"]["clientWidth"] + 1,
        f"index {mode}: horizontal overflow {info['html']}",
    )
    require(
        page.locator('a.sample-link').count() == 7,
        f"index {mode}: expected seven sample launchers",
    )
    page.screenshot(path=str(output / f"index_{mode}.png"), full_page=False)
    require(not console_errors, f"index {mode}: console errors {console_errors}")
    require(not external_requests, f"index {mode}: external requests {external_requests}")
    return {
        **info,
        "sample_links": page.locator("a.sample-link").count(),
        "console_errors": console_errors,
        "external_requests": external_requests,
    }


def audit_sample(page: Page, sample: str, output: Path, mode: str) -> dict[str, Any]:
    console_errors, external_requests = page_observers(page)
    page.goto((WORKSTATION / f"{sample}.html").as_uri(), wait_until="domcontentloaded")
    page.wait_for_selector('[data-testid="genome-canvas"]', timeout=60_000)
    page.wait_for_function(
        "() => document.getElementById('filterStatus')?.textContent.includes('regions')",
        timeout=60_000,
    )
    info = geometry(page)
    require(info["authority"] == AUTHORITY, f"{sample} {mode}: authority mismatch")
    require(info["contract"] == UI_CONTRACT, f"{sample} {mode}: UI contract mismatch")
    require(
        info["html"]["scrollWidth"] <= info["html"]["clientWidth"] + 1,
        f"{sample} {mode}: horizontal overflow {info['html']}",
    )
    require(
        page.locator('[data-testid="sample-metrics"] .metric').count() == 5,
        f"{sample} {mode}: KPI card count mismatch",
    )

    page.locator("#genome").scroll_into_view_if_needed()
    page.wait_for_timeout(250)
    legend_buttons = page.locator("#legend .legend-btn")
    require(legend_buttons.count() >= 3, f"{sample} {mode}: legend too small")
    legend_buttons.nth(0).click()
    legend_buttons.nth(1).click()
    require(
        legend_buttons.nth(0).get_attribute("aria-pressed") == "true"
        and legend_buttons.nth(1).get_attribute("aria-pressed") == "true",
        f"{sample} {mode}: multiselect union failed",
    )
    legend_buttons.nth(0).click()
    require(
        legend_buttons.nth(0).get_attribute("aria-pressed") == "false",
        f"{sample} {mode}: second-click deselect failed",
    )

    page.locator("#clearSearch").click()
    first_region = page.locator("#regionList .region-button").first
    require(first_region.count() == 1, f"{sample} {mode}: no region list item")
    first_region.click()
    page.wait_for_selector("#regionDetail .detail-id", timeout=10_000)
    require(
        page.locator("#regionDetail .detail-id").inner_text().startswith("chr"),
        f"{sample} {mode}: region detail did not render",
    )

    if sample == "HCC1395":
        for button in page.locator(
            '#legend .legend-btn[aria-pressed="true"]'
        ).all():
            button.click()
        page.locator("#regionSearch").fill("chr10:87818272-87928739")
        page.locator("#searchForm button[type=submit]").click()
        page.wait_for_timeout(250)
        require(
            "87818272" in page.locator("#regionDetail").inner_text()
            and "87840023" in page.locator("#regionDetail").inner_text(),
            f"{sample} {mode}: target exact group not found",
        )
        require(
            "87888228" not in page.locator("#regionDetail .detail-id").inner_text(),
            f"{sample} {mode}: legacy four-site region leaked into detail key",
        )

    page.evaluate(
        """() => {
          document.documentElement.style.scrollBehavior='auto';
          const target=document.querySelector('#genome');
          window.scrollTo({left:0,top:target.offsetTop,behavior:'auto'});
        }"""
    )
    page.wait_for_timeout(150)
    page.screenshot(path=str(output / f"{sample}_{mode}_genome.png"), full_page=False)
    page.locator("#regionDetail").screenshot(
        path=str(output / f"{sample}_{mode}_detail.png")
    )
    page.evaluate("window.scrollTo({left:0,top:0,behavior:'auto'})")
    page.wait_for_timeout(200)
    page.screenshot(path=str(output / f"{sample}_{mode}_top.png"), full_page=False)
    require(not console_errors, f"{sample} {mode}: console errors {console_errors}")
    require(not external_requests, f"{sample} {mode}: external requests {external_requests}")
    return {
        **info,
        "legend_buttons": legend_buttons.count(),
        "filter_status": page.locator("#filterStatus").inner_text(),
        "detail_id": page.locator("#regionDetail .detail-id").inner_text(),
        "console_errors": console_errors,
        "external_requests": external_requests,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--sample", choices=SAMPLES, action="append")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    samples = tuple(args.sample) if args.sample else SAMPLES
    receipt: dict[str, Any] = {
        "schema_name": "intersubmod.exact_ps_layered_workstation.browser_audit",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "authority": AUTHORITY,
        "ui_contract": UI_CONTRACT,
        "samples": list(samples),
        "viewports": VIEWPORTS,
        "checks": {},
        "pages": {},
    }
    try:
        with sync_playwright() as playwright:
            browser = launch_browser(playwright)
            try:
                for mode, viewport in VIEWPORTS.items():
                    context = browser.new_context(viewport=viewport)
                    try:
                        page = context.new_page()
                        receipt["pages"][f"index:{mode}"] = audit_index(
                            page, output, mode
                        )
                        page.close()
                        for sample in samples:
                            page = context.new_page()
                            receipt["pages"][f"{sample}:{mode}"] = audit_sample(
                                page, sample, output, mode
                            )
                            page.close()
                    finally:
                        context.close()
            finally:
                browser.close()
        receipt["checks"] = {
            "all_pages_audited": len(receipt["pages"]) == 2 * (1 + len(samples)),
            "no_console_errors": all(
                not row["console_errors"] for row in receipt["pages"].values()
            ),
            "no_external_requests": all(
                not row["external_requests"] for row in receipt["pages"].values()
            ),
            "no_horizontal_overflow": all(
                row["html"]["scrollWidth"] <= row["html"]["clientWidth"] + 1
                for row in receipt["pages"].values()
            ),
        }
        receipt["all_pass"] = all(receipt["checks"].values())
    except Exception as exc:
        receipt["all_pass"] = False
        receipt["failure"] = str(exc)
    receipt_path = output / "receipt.json"
    receipt_path.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"INPUT={WORKSTATION}")
    print(f"OUTPUT={receipt_path}")
    print(json.dumps({"all_pass": receipt["all_pass"], **receipt["checks"]}))
    if not receipt["all_pass"]:
        print(f"FAIL: {receipt.get('failure', receipt['checks'])}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
