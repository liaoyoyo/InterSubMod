#!/usr/bin/env python3
"""Independent structural, numerical, portability, and browser QA."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from playwright.sync_api import sync_playwright


REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
WORK = REPO / "research/20260724_methyl_linked_ssnv_g1_visual_report"
HTML = WORK / "20260724_methyl_linked_sSNV_G1全樣本關聯圖譜_01.standalone.html"
DATA = WORK / "data/report_data.json"
RESULTS = WORK / "results"
RECEIPT = RESULTS / "20260724_G1全樣本HTML_QA收據_01.json"
TOP_SCREENSHOT = RESULTS / "01_G1_atlas_top_desktop.png"
CASE_SCREENSHOT = RESULTS / "02_G1_HCC1395_case_desktop.png"
MOBILE_SCREENSHOT = RESULTS / "03_G1_atlas_mobile.png"
CHROMIUM = Path(
    "/bip7_disk/liaoyoyo2001/.cache/ms-playwright/chromium-1223/chrome-linux64/chrome"
)


def structural_qa(data: dict) -> dict:
    cases = data["cases"]
    assert data["task_type"] == "B_comprehensive_validation"
    assert data["scope"]["partial"] is False
    assert data["scope"]["technical_datasets"] == 7
    assert data["scope"]["biological_samples"] == 6
    assert len(cases) == 7
    assert data["funnel"]["all_sites"] == 469849
    assert data["funnel"]["m1"] == 102842
    assert data["funnel"]["m2_evaluable"] == 1867
    assert data["funnel"]["m2_pass"] == 919
    assert data["funnel"]["exact_pairs"] == 147
    assert data["funnel"]["by"] == 10
    assert data["funnel"]["g1"] == 7
    assert data["funnel"]["g2_family_sites"] == 58
    assert data["headline"]["g2_global_by"] == 0
    assert data["headline"]["identified_mutation_orders"] == 0
    assert data["headline"]["strict_w_containers"] == 10
    assert data["headline"]["strict_direct_support"] == 790
    assert sum(case["methylation"]["core_read_rows"] for case in cases) == 826
    assert sum(len(case["strict_w"]) for case in cases) == 10
    assert sum(
        evidence["direct_support"]
        for case in cases
        for evidence in case["strict_w"]
    ) == 790
    assert all(case["callability"]["pass"] for case in cases)
    assert all(case["statistics"]["conditional_p"] == 0.001 for case in cases)
    assert all(case["endpoint_b"]["n_compatible_relation_models"] == 0 for case in cases)
    assert all(len(case["methylation"]["rows"]) == case["methylation"]["core_read_rows"] for case in cases)
    assert all(len(case["methylation"]["display_cpg_positions"]) <= 40 for case in cases)
    assert all(
        all(group["group"].startswith("MG-") for group in case["groups"])
        for case in cases
    )
    assert sum(case["partner"]["truth"] == "FP" for case in cases) == 1
    assert sum(row["g1"] for row in data["dataset_rows"]) == 7

    html = HTML.read_text(encoding="utf-8")
    assert "<script src=" not in html
    assert "<link rel=" not in html
    assert "http://" not in html
    assert "https://" not in html
    required = [
        "可以確認「哪個 linked 位點和甲基群有關」；不能確認誰造成誰",
        "0 IDENTIFIED ORDERS",
        "MG 是甲基群，不是 HP",
        "非演化樹",
        "NOT_SIGNED_TASK_B_FINAL_RELEASE",
    ]
    # The release status is embedded as JSON rather than visible prose.
    assert all(text in html for text in required)
    forbidden = [
        "已證明是 subclone",
        "confirmed evolutionary tree",
        "已證明甲基造成突變",
        "已確認父子 clone",
    ]
    assert not any(text in html for text in forbidden)
    return {
        "cases": len(cases),
        "core_reads": sum(case["methylation"]["core_read_rows"] for case in cases),
        "strict_w_containers": sum(len(case["strict_w"]) for case in cases),
        "standalone_no_external_assets": True,
        "claim_guardrail_terms_present": True,
    }


def browser_qa(browser, name: str, width: int, height: int, capture: bool) -> dict:
    console_errors: list[str] = []
    page_errors: list[str] = []
    external_requests: list[str] = []
    context = browser.new_context(
        viewport={"width": width, "height": height},
        color_scheme="light",
        reduced_motion="reduce",
    )
    page = context.new_page()
    page.set_default_timeout(60_000)
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
        if request.url.startswith(("http://", "https://", "ws://", "wss://"))
        else None,
    )
    page.goto(HTML.as_uri(), wait_until="load", timeout=120_000)
    page.add_style_tag(content="html{scroll-behavior:auto!important}")
    page.wait_for_timeout(1000)

    title = page.title()
    case_count = page.locator(".case").count()
    index_count = page.locator("#case-index tbody tr").count()
    canvas_count = page.locator("canvas.heatmap").count()
    if case_count != 7 or index_count != 7 or canvas_count != 7:
        raise AssertionError(
            f"{name}: cases={case_count}, index={index_count}, canvases={canvas_count}"
        )
    heading = " ".join(page.locator("h1").first.inner_text().split())
    if heading != "甲基群 × linked-sSNV 全樣本 G1 關聯圖譜":
        raise AssertionError(f"{name}: title heading mismatch: {heading!r}")

    canvas_checks = page.locator("canvas.heatmap").evaluate_all(
        """canvases => canvases.map(canvas => {
          const rect = canvas.getBoundingClientRect();
          const ctx = canvas.getContext('2d');
          const sample = ctx.getImageData(0, 0, canvas.width, canvas.height).data;
          let nonTransparent = 0;
          for (let i = 3; i < sample.length; i += Math.max(4, Math.floor(sample.length / 5000 / 4) * 4)) {
            if (sample[i] > 0) nonTransparent++;
          }
          return {
            cssWidth: rect.width,
            cssHeight: rect.height,
            bitmapWidth: canvas.width,
            bitmapHeight: canvas.height,
            nonTransparent,
          };
        })"""
    )
    if any(
        item["cssWidth"] < 250
        or item["cssHeight"] < 200
        or item["bitmapWidth"] < 250
        or item["bitmapHeight"] < 200
        or item["nonTransparent"] < 25
        for item in canvas_checks
    ):
        raise AssertionError(f"{name}: heatmap canvas failed {canvas_checks}")

    geometry = page.evaluate(
        """() => ({
          clientWidth: document.documentElement.clientWidth,
          scrollWidth: document.documentElement.scrollWidth,
          zeroCases: [...document.querySelectorAll('.case')].filter(node => {
            const r=node.getBoundingClientRect(); return r.width <= 0 || r.height <= 0;
          }).length,
          zeroPanes: [...document.querySelectorAll('.pane')].filter(node => {
            const r=node.getBoundingClientRect(); return r.width <= 0 || r.height <= 0;
          }).length
        })"""
    )
    if geometry["scrollWidth"] > geometry["clientWidth"] + 1:
        raise AssertionError(f"{name}: page horizontal overflow {geometry}")
    if geometry["zeroCases"] or geometry["zeroPanes"]:
        raise AssertionError(f"{name}: zero-size sections {geometry}")

    filter_result = None
    if name == "desktop":
        page.locator('button[data-filter="HCC1395"]').click()
        visible_cases = page.locator(".case:visible").count()
        if visible_cases != 1:
            raise AssertionError(f"Dataset filter left {visible_cases} cases")
        filter_result = "PASS"
        page.locator('button[data-filter="全部"]').click()

    if capture:
        page.evaluate("window.scrollTo(0,0)")
        page.wait_for_timeout(100)
        page.screenshot(path=str(TOP_SCREENSHOT), full_page=False)
        first_case = page.locator(".case").first
        page.evaluate(
            "node => window.scrollTo(0, Math.max(0, node.offsetTop - 64))",
            first_case.element_handle(),
        )
        page.wait_for_timeout(100)
        page.screenshot(path=str(CASE_SCREENSHOT), full_page=False)
    if name == "mobile":
        page.evaluate("window.scrollTo(0,0)")
        page.screenshot(path=str(MOBILE_SCREENSHOT), full_page=False)

    if console_errors or page_errors or external_requests:
        raise AssertionError(
            f"{name}: console={console_errors}, page={page_errors}, external={external_requests}"
        )
    result = {
        "name": name,
        "viewport": {"width": width, "height": height},
        "title": title,
        "case_count": case_count,
        "index_count": index_count,
        "canvas_count": canvas_count,
        "canvas_checks": canvas_checks,
        "geometry": geometry,
        "filter_interaction": filter_result,
        "console_errors": console_errors,
        "page_errors": page_errors,
        "external_requests": external_requests,
    }
    context.close()
    return result


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    if not HTML.exists() or not DATA.exists() or not CHROMIUM.exists():
        raise FileNotFoundError("HTML, report data, or Chromium is missing")
    data = json.loads(DATA.read_text(encoding="utf-8"))
    structural = structural_qa(data)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            headless=True,
            executable_path=str(CHROMIUM),
            args=["--no-sandbox", "--disable-dev-shm-usage", "--disable-gpu"],
        )
        try:
            viewports = [
                browser_qa(browser, "desktop", 1440, 1050, True),
                browser_qa(browser, "mobile", 390, 844, False),
            ]
        finally:
            browser.close()
    receipt = {
        "schema_name": "intersubmod.g1_html_qa_receipt",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS",
        "html": str(HTML),
        "data": str(DATA),
        "browser": str(CHROMIUM),
        "structural": structural,
        "viewports": viewports,
        "screenshots": [
            str(TOP_SCREENSHOT),
            str(CASE_SCREENSHOT),
            str(MOBILE_SCREENSHOT),
        ],
    }
    RECEIPT.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": receipt["status"],
                "receipt": str(RECEIPT),
                "desktop_cases": viewports[0]["case_count"],
                "mobile_cases": viewports[1]["case_count"],
                "console_errors": 0,
                "external_requests": 0,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
