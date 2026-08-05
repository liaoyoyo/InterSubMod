#!/usr/bin/env python3
"""Browser-level QA for the self-contained computation-complexity explainer."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List

from playwright.sync_api import sync_playwright


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check(name: str, observed: Any, expected: Any) -> Dict[str, Any]:
    return {
        "name": name,
        "observed": observed,
        "expected": expected,
        "pass": observed == expected,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--screenshot-dir", required=True, type=Path)
    args = parser.parse_args()

    html_path = args.html.resolve()
    receipt_path = args.receipt.resolve()
    screenshot_dir = args.screenshot_dir.resolve()
    if not html_path.is_file():
        raise SystemExit(f"HTML not found: {html_path}")
    if receipt_path.exists():
        raise SystemExit(f"Refusing to overwrite receipt: {receipt_path}")
    screenshot_dir.mkdir(parents=True, exist_ok=True)

    checks: List[Dict[str, Any]] = []
    console_errors: List[str] = []
    page_errors: List[str] = []
    external_requests: List[str] = []
    screenshots: List[str] = []

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1440, "height": 1000})
        page.on(
            "console",
            lambda message: console_errors.append(message.text)
            if message.type == "error"
            else None,
        )
        page.on("pageerror", lambda error: page_errors.append(str(error)))
        page.on(
            "request",
            lambda request: external_requests.append(request.url)
            if request.url.startswith(("http://", "https://"))
            else None,
        )

        response = page.goto(html_path.as_uri(), wait_until="networkidle")
        checks.append(
            check(
                "file_navigation_status",
                None if response is None else response.status,
                200,
            )
        )
        checks.append(
            check(
                "verified_evidence_ready",
                page.locator("html").get_attribute("data-evidence-ready"),
                "true",
            )
        )
        checks.append(check("single_h1", page.locator("h1").count(), 1))
        checks.append(check("major_sections", page.locator("main > section").count(), 11))
        checks.append(
            check("formula_rows_q0_to_q12", page.locator("#formula-body tr").count(), 13)
        )
        checks.append(
            check(
                "observed_rows_q0_to_q12",
                page.locator("#observed-q-body tr").count(),
                13,
            )
        )
        checks.append(
            check(
                "completion_chart_has_13_complete_bars",
                page.locator("#q-completion-chart rect").count() >= 13,
                True,
            )
        )
        body_text = page.locator("body").inner_text()
        for token in (
            "98,955",
            "35,754,927",
            "85,621",
            "99,966",
            "83,252",
            "10,717",
            "39,648",
            "23,858",
            "8,449",
            "63,506",
            "88.26%",
        ):
            checks.append(check(f"key_text_{token}", token in body_text, True))
        checks.append(
            check(
                "coarse_shape_claim_is_explicit",
                "單一 rooted-unlabeled 粗形狀" in body_text,
                True,
            )
        )
        checks.append(
            check(
                "forbidden_exact_topology_aggregate_claim_absent",
                "可確認一個 exact topology" in body_text,
                False,
            )
        )
        checks.append(
            check(
                "current_flow_does_not_claim_timeout",
                "cap / timeout / family" in body_text,
                False,
            )
        )
        checks.append(
            check(
                "extreme_partial_joint_exponent",
                page.locator("#extreme-partial-copy sup").inner_text(),
                "581",
            )
        )

        desktop_overflow = page.evaluate(
            "() => document.documentElement.scrollWidth - window.innerWidth"
        )
        checks.append(check("desktop_horizontal_overflow_px", desktop_overflow <= 2, True))
        desktop_shot = screenshot_dir / "01_desktop_full.png"
        page.screenshot(path=str(desktop_shot), full_page=False)
        screenshots.append(str(desktop_shot))

        page.locator("#q-number").fill("12")
        page.locator("#u-list").fill("11,10")
        page.locator("#cv-input").fill("360")
        page.locator("#parent-list").fill("1,2,3")
        checks.append(check("calculator_cube_q12", page.locator("#calc-cube").inner_text(), "4,096 / 24,576"))
        checks.append(check("calculator_partial_sum", page.locator("#calc-partial").inner_text(), "3,072"))
        checks.append(check("calculator_fixed_tree_count", page.locator("#calc-tree").inner_text(), "6"))
        checks.append(
            check(
                "calculator_family_tree_illustration",
                page.locator("#calc-family-tree").inner_text(),
                "2,160",
            )
        )
        checks.append(
            check(
                "calculator_q12_empirical_zero_complete",
                "0%" in page.locator("#calc-empirical").inner_text(),
                True,
            )
        )
        page.locator("#u-list").fill("13")
        checks.append(check("calculator_invalid_u_visible", page.locator("#calc-error").is_visible(), True))
        checks.append(
            check(
                "calculator_invalid_u_message",
                "不可大於 q=12" in page.locator("#calc-error").inner_text(),
                True,
            )
        )

        page.set_viewport_size({"width": 390, "height": 844})
        page.reload(wait_until="networkidle")
        mobile_overflow = page.evaluate(
            "() => document.documentElement.scrollWidth - window.innerWidth"
        )
        checks.append(check("mobile_horizontal_overflow_px", mobile_overflow <= 2, True))
        checks.append(
            check(
                "mobile_verified_evidence_ready",
                page.locator("html").get_attribute("data-evidence-ready"),
                "true",
            )
        )
        mobile_shot = screenshot_dir / "02_mobile_full.png"
        page.screenshot(path=str(mobile_shot), full_page=False)
        screenshots.append(str(mobile_shot))

        page.emulate_media(media="print")
        checks.append(
            check(
                "print_footer_present",
                "Task Type F" in page.locator("footer").inner_text(),
                True,
            )
        )
        print_shot = screenshot_dir / "03_print_first_view.png"
        page.screenshot(path=str(print_shot), full_page=False)
        screenshots.append(str(print_shot))
        browser.close()

    checks.append(check("console_errors", console_errors, []))
    checks.append(check("page_errors", page_errors, []))
    checks.append(check("external_network_requests", external_requests, []))
    all_pass = all(item["pass"] for item in checks)
    receipt = {
        "schema_version": "1.0",
        "purpose": "Browser QA for the self-contained Hypercube computation explainer",
        "html": {
            "path": str(html_path),
            "sha256": sha256(html_path),
            "size_bytes": html_path.stat().st_size,
        },
        "viewports": {
            "desktop": {"width": 1440, "height": 1000},
            "mobile": {"width": 390, "height": 844},
        },
        "screenshots": screenshots,
        "checks": checks,
        "all_pass": all_pass,
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "all_pass": all_pass,
                "checks_passed": sum(item["pass"] for item in checks),
                "checks_total": len(checks),
                "html": str(html_path),
                "receipt": str(receipt_path),
                "screenshots": screenshots,
            },
            ensure_ascii=False,
        )
    )
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
