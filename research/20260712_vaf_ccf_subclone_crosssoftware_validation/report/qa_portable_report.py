#!/usr/bin/env python3
"""Render and smoke-test the portable validation report with Playwright."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from playwright.sync_api import sync_playwright


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    html = args.html.resolve()
    if not html.is_file():
        raise FileNotFoundError(html)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    console_errors: list[str] = []
    page_errors: list[str] = []
    checks: dict[str, bool] = {}
    rendered: dict[str, object] = {}

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        for name, width, height in (
            ("desktop", 1440, 1000),
            ("mobile", 390, 844),
        ):
            page = browser.new_page(viewport={"width": width, "height": height})
            page.on(
                "console",
                lambda message: console_errors.append(message.text)
                if message.type == "error"
                else None,
            )
            page.on("pageerror", lambda error: page_errors.append(str(error)))
            page.goto(html.as_uri(), wait_until="networkidle")
            page.wait_for_timeout(1000)
            body_text = page.locator("body").inner_text()
            visible_h1_texts = page.locator("h1").evaluate_all(
                "elements => elements.filter(element => { const rect = element.getBoundingClientRect(); "
                "return rect.width > 0 && rect.height > 0; }).map(element => element.textContent.trim())"
            )
            screenshot = output_dir / f"report_{name}.png"
            page.screenshot(path=str(screenshot), full_page=True)
            rendered[name] = {
                "viewport": {"width": width, "height": height},
                "body_text_chars": len(body_text),
                "h1_count": page.locator("h1").count(),
                "visible_h1_count": len(visible_h1_texts),
                "visible_h1_texts": visible_h1_texts,
                "table_count": page.locator("table").count(),
                "svg_count": page.locator("svg").count(),
                "canvas_count": page.locator("canvas").count(),
                "screenshot": str(screenshot),
                "screenshot_sha256": sha256_file(screenshot),
            }
            checks[f"{name}_body_nonempty"] = len(body_text) > 5000
            checks[f"{name}_h1_consistent"] = (
                len(visible_h1_texts) >= 1 and len(set(visible_h1_texts)) == 1
            )
            checks[f"{name}_title_present"] = "HCC1395" in body_text and "VAF" in body_text
            checks[f"{name}_pyclone_present"] = "PyClone-VI" in body_text
            checks[f"{name}_topology_present"] = "5,720" in body_text and "69.39%" in body_text
            checks[f"{name}_minor_metric_present"] = "0.381" in body_text
            checks[f"{name}_regional_precision_present"] = (
                "15,713/15,713" in body_text
                and "96.75%" in body_text
                and "54.62%" in body_text
            )
            checks[f"{name}_clone_region_bridge_present"] = (
                "14,369" in body_text
                and "61.76%" in body_text
                and "subclonal union" in body_text
            )
            checks[f"{name}_cause_audit_present"] = (
                "703 regions" in body_text
                and "identifiability" in body_text
                and "OR=0.921" in body_text
            )
            checks[f"{name}_claim_ceiling_present"] = (
                "不支持「每區域已準確回復唯一真實演化樹」" in body_text
            )
            checks[f"{name}_tables_present"] = page.locator("table").count() >= 12
            checks[f"{name}_visuals_present"] = (
                page.locator("svg").count() + page.locator("canvas").count()
            ) >= 6
            page.close()
        browser.close()

    checks["console_errors_absent"] = not console_errors
    checks["page_errors_absent"] = not page_errors
    receipt = {
        "schema_name": "intersubmod.portable_report_browser_qa",
        "schema_version": "1.0.0",
        "html": str(html),
        "html_sha256": sha256_file(html),
        "checks": checks,
        "console_errors": console_errors,
        "page_errors": page_errors,
        "rendered": rendered,
        "status": "PASS" if all(checks.values()) else "FAIL",
    }
    receipt_path = output_dir / "browser_qa_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"status={receipt['status']} checks={sum(checks.values())}/{len(checks)}")
    print(f"receipt={receipt_path}")
    if receipt["status"] != "PASS":
        for name, passed in checks.items():
            if not passed:
                print(f"FAIL {name}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
