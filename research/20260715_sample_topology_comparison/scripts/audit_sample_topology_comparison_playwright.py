#!/usr/bin/env python3
"""Chromium/Playwright audit for the sample topology comparison dashboard."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from playwright.sync_api import Browser, Page, sync_playwright


SCRIPT_PATH = Path(__file__).resolve()
TOPIC_ROOT = SCRIPT_PATH.parents[1]
REPO_ROOT = SCRIPT_PATH.parents[3]
DEFAULT_INDEX = REPO_ROOT / "docs/methodology/_assets/layered_workstation/index.html"
DEFAULT_COMPARISON = TOPIC_ROOT / "artifacts/sample_topology_comparison.json"
DEFAULT_OUTPUT = TOPIC_ROOT / "qa/playwright"
VIEWPORTS = {
    "desktop": {"width": 1440, "height": 1000},
    "mobile": {"width": 390, "height": 844},
}
DIMENSIONS = ("structural", "read_af_selection", "morphology")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def add_check(run: dict[str, Any], name: str, passed: bool, expected: Any, actual: Any) -> None:
    run["checks"].append(
        {
            "name": name,
            "status": "pass" if passed else "fail",
            "expected": expected,
            "actual": actual,
        }
    )


def scroll_to(page: Page, selector: str) -> None:
    # The workstation intentionally uses smooth scrolling for readers. Disable it
    # inside the audit so a viewport screenshot cannot capture an intermediate
    # animation frame and falsely suggest that rows are missing.
    page.evaluate(
        """selector => {
          const element = document.querySelector(selector);
          if (!element) return;
          document.documentElement.style.scrollBehavior = 'auto';
          window.scrollTo(0, Math.max(0, element.getBoundingClientRect().top + window.scrollY - 12));
        }""",
        selector,
    )
    page.wait_for_function(
        "selector => { const element=document.querySelector(selector); return element && Math.abs(element.getBoundingClientRect().top - 12) <= 2; }",
        arg=selector,
    )
    page.wait_for_timeout(80)


def screenshot(page: Page, path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    page.screenshot(path=str(path), full_page=False)
    return str(path.resolve())


def dimension_state(page: Page) -> dict[str, Any]:
    return page.evaluate(
        """() => {
          const embedded = JSON.parse(document.getElementById('sample-comparison-data').textContent);
          const active = document.querySelector('.comparison-tab[aria-pressed="true"]');
          const profileRows = [...document.querySelectorAll('#comparison-profile-chart .comparison-profile-row')];
          const matrix = document.querySelector('.distance-matrix');
          const matrixRows = [...(matrix?.querySelectorAll('tbody tr') || [])];
          const matrixValues = matrixRows.map(row => [...row.querySelectorAll('td')].map(cell => {
            const button = cell.querySelector('button');
            return button ? Number(button.textContent) : 0;
          }));
          const confusion = document.querySelector('.confusion-table');
          const confusionValues = [...(confusion?.querySelectorAll('tbody td') || [])]
            .map(cell => Number((cell.textContent || '').replaceAll(',', '')));
          const ledger = document.querySelector('.comparison-table');
          const segmentSums = profileRows.map(row => [...row.querySelectorAll('.bar-segment')]
            .reduce((sum, segment) => sum + Number.parseFloat(segment.style.width || '0'), 0));
          return {
            embedded_dataset_count: embedded.datasets.length,
            embedded_biological_ids: [...new Set(embedded.datasets.map(row => row.biological_id))].sort(),
            active_dimension: active?.dataset.comparisonDimension || null,
            active_tab_count: document.querySelectorAll('.comparison-tab[aria-pressed="true"]').length,
            profile_count: profileRows.length,
            profile_segment_sums: segmentSums,
            ledger_rows: ledger?.querySelectorAll('tbody tr').length || 0,
            ledger_columns: ledger?.querySelectorAll('thead th').length || 0,
            matrix_rows: matrixRows.length,
            matrix_columns: matrixRows[0]?.querySelectorAll('td').length || 0,
            matrix_values: matrixValues,
            hcc_matrix_value: Number(document.querySelector('.matrix-cell[data-left="HCC1395"][data-right="HCC1395_DORADO"]')?.textContent || NaN),
            confusion_rows: confusion?.querySelectorAll('tbody tr').length || 0,
            confusion_columns: confusion?.querySelectorAll('thead th').length || 0,
            confusion_sum: confusionValues.reduce((sum, value) => sum + value, 0),
            comparison_title: document.getElementById('comparison-profile-title')?.textContent || '',
            confusion_title: document.getElementById('hcc-confusion-title')?.textContent || ''
          };
        }"""
    )


def audit_viewport(
    browser: Browser,
    index_path: Path,
    comparison: dict[str, Any],
    comparison_sha256: str,
    output_dir: Path,
    viewport_name: str,
    timeout_ms: int,
) -> dict[str, Any]:
    started = time.monotonic()
    run: dict[str, Any] = {
        "viewport": viewport_name,
        "viewport_size": VIEWPORTS[viewport_name],
        "input": str(index_path.resolve()),
        "checks": [],
        "console_errors": [],
        "page_errors": [],
        "screenshots": {},
    }
    context = browser.new_context(
        viewport=VIEWPORTS[viewport_name],
        locale="zh-TW",
        color_scheme="light",
        device_scale_factor=1,
    )
    page = context.new_page()
    page.set_default_timeout(timeout_ms)
    page.on("console", lambda message: run["console_errors"].append(message.text) if message.type == "error" else None)
    page.on("pageerror", lambda error: run["page_errors"].append(str(error)))
    try:
        page.goto(index_path.resolve().as_uri(), wait_until="load", timeout=timeout_ms)
        page.wait_for_load_state("networkidle", timeout=timeout_ms)
        add_check(run, "document_load", True, "load + networkidle", page.evaluate("document.readyState"))

        meta = page.locator('meta[name="intersubmod-sample-topology-comparison-sha256"]').get_attribute("content")
        add_check(run, "comparison_sha256_meta_binding", meta == comparison_sha256, comparison_sha256, meta)

        initial = dimension_state(page)
        add_check(run, "dataset_and_biological_counts", initial["embedded_dataset_count"] == 7 and len(initial["embedded_biological_ids"]) == 6, {"datasets": 7, "biological_ids": 6}, initial)

        dimension_receipts: dict[str, Any] = {}
        order = comparison["pairwise_composition"]["matrices"]["order"]
        for dimension in DIMENSIONS:
            page.locator(f'.comparison-tab[data-comparison-dimension="{dimension}"]').click()
            page.wait_for_timeout(80)
            state = dimension_state(page)
            expected_categories = comparison["dimension_contracts"][dimension]["categories"]
            expected_hcc = comparison["hcc1395_technical_pair"]["marginal_profiles"]["dimensions"][dimension]["tvd"]
            expected_matched = comparison["hcc1395_technical_pair"]["matched_primary_dimensions"][dimension]["n"]
            symmetric = all(
                abs(state["matrix_values"][i][j] - state["matrix_values"][j][i]) <= 0.00051
                for i in range(7)
                for j in range(7)
            )
            diagonal_zero = all(abs(state["matrix_values"][i][i]) <= 1e-12 for i in range(7))
            passed = (
                state["active_dimension"] == dimension
                and state["active_tab_count"] == 1
                and state["profile_count"] == 7
                and all(abs(total - 100.0) <= 0.001 for total in state["profile_segment_sums"])
                and state["ledger_rows"] == 7
                and state["ledger_columns"] == len(expected_categories) + 2
                and state["matrix_rows"] == 7
                and state["matrix_columns"] == 7
                and symmetric
                and diagonal_zero
                and abs(state["hcc_matrix_value"] - expected_hcc) <= 0.00051
                and state["confusion_rows"] == len(expected_categories)
                and state["confusion_columns"] == len(expected_categories) + 1
                and state["confusion_sum"] == expected_matched
            )
            dimension_receipts[dimension] = state
            add_check(
                run,
                f"dimension_sync_{dimension}",
                passed,
                {
                    "active": dimension,
                    "profiles": 7,
                    "ledger_rows": 7,
                    "matrix": "7x7 symmetric diagonal=0",
                    "hcc_tvd": expected_hcc,
                    "confusion_sum": expected_matched,
                },
                state,
            )

        page.locator('.comparison-tab[data-comparison-dimension="structural"]').click()
        page.locator('.matrix-cell[data-left="COLO829"][data-right="H1437"]').click()
        regular_text = page.locator("#pair-inspector").inner_text()
        regular_upper = regular_text.upper()
        add_check(
            run,
            "regular_pair_inspector_has_no_technical_claim",
            "DIFFERENT BIOLOGICAL SAMPLES" in regular_upper and "TECHNICAL PAIR" not in regular_upper,
            "different biological samples; composition only",
            regular_text,
        )
        page.locator('.matrix-cell[data-left="HCC1395"][data-right="HCC1395_DORADO"]').click()
        hcc_text = page.locator("#pair-inspector").inner_text()
        hcc_upper = hcc_text.upper()
        add_check(
            run,
            "hcc_pair_inspector_relationship",
            "SAME BIOLOGICAL SAMPLE" in hcc_upper and "TECHNICAL PAIR" in hcc_upper,
            "same biological sample technical pair",
            hcc_text,
        )

        touch_sizes = page.evaluate(
            """() => ({
              tabs: [...document.querySelectorAll('.comparison-tab')].map(element => element.getBoundingClientRect().height),
              cells: [...document.querySelectorAll('.matrix-cell')].slice(0, 8).map(element => element.getBoundingClientRect().height)
            })"""
        )
        add_check(
            run,
            "interactive_targets_at_least_44px",
            all(value >= 44 for value in touch_sizes["tabs"] + touch_sizes["cells"]),
            {"minimum_height": 44},
            touch_sizes,
        )

        layout = page.evaluate(
            """() => {
              const matrix = document.querySelector('.heatmap-scroll');
              const links = [...document.querySelectorAll('a[href]')].filter(link => (link.getAttribute('href') || '').toLowerCase().includes('.json'));
              return {
                body_overflow_px: Math.max(0, document.documentElement.scrollWidth - document.documentElement.clientWidth),
                matrix_client_width: matrix?.clientWidth || 0,
                matrix_scroll_width: matrix?.scrollWidth || 0,
                json_links: links.length,
                json_links_inside_details: links.filter(link => link.closest('details')).length,
                json_links_in_open_details: links.filter(link => link.closest('details')?.open).length,
                hcc_metric_cards: document.querySelectorAll('.hcc-metric-grid article').length,
                exact_edge_visible: (document.querySelector('#hcc1395-technical-validation')?.innerText || '').includes('35.3%')
              };
            }"""
        )
        add_check(run, "body_and_internal_overflow", layout["body_overflow_px"] <= 1 and layout["matrix_scroll_width"] >= layout["matrix_client_width"], {"body_overflow_max": 1, "matrix_internal_scroll_allowed": True}, layout)
        add_check(run, "json_links_hidden_in_closed_details", layout["json_links"] > 0 and layout["json_links_inside_details"] == layout["json_links"] and layout["json_links_in_open_details"] == 0, "all JSON links inside closed details", layout)
        add_check(run, "hcc_exact_signature_visible", layout["hcc_metric_cards"] == 4 and layout["exact_edge_visible"], {"metric_cards": 4, "exact_edge": "35.3%"}, layout)

        profile_geometry = page.evaluate(
            """() => {
              const list = document.querySelector('#comparison-profile-chart');
              const listRect = list?.getBoundingClientRect();
              return [...document.querySelectorAll('#comparison-profile-chart .comparison-profile-row')].map(row => {
                const rowRect = row.getBoundingClientRect();
                const barRect = row.querySelector('.topology-bar')?.getBoundingClientRect();
                const nameRect = row.querySelector('.comparison-profile-name')?.getBoundingClientRect();
                const denominatorRect = row.querySelector('.comparison-profile-n')?.getBoundingClientRect();
                return {
                  dataset: row.dataset.dataset,
                  row_height: rowRect.height,
                  bar_width: barRect?.width || 0,
                  bar_height: barRect?.height || 0,
                  name_visible: (nameRect?.width || 0) > 0 && (nameRect?.height || 0) > 0,
                  denominator_visible: (denominatorRect?.width || 0) > 0 && (denominatorRect?.height || 0) > 0,
                  within_list: Boolean(listRect) && rowRect.top >= listRect.top - 1 && rowRect.bottom <= listRect.bottom + 1
                };
              });
            }"""
        )
        geometry_ok = len(profile_geometry) == 7 and all(
            row["bar_width"] > 0
            and row["bar_height"] >= 16
            and row["name_visible"]
            and row["denominator_visible"]
            and row["within_list"]
            and row["row_height"] <= 140
            for row in profile_geometry
        )
        add_check(
            run,
            "all_profile_rows_visually_present",
            geometry_ok,
            {"rows": 7, "bar_height_min": 16, "row_height_max": 140, "all_children_visible": True},
            profile_geometry,
        )

        scroll_to(page, "#cohort-comparison")
        comparison_top = page.locator("#cohort-comparison").evaluate("element => element.getBoundingClientRect().top")
        run["screenshots"]["comparison"] = screenshot(page, output_dir / f"comparison__{viewport_name}.png")
        scroll_to(page, "#hcc1395-technical-validation")
        hcc_top = page.locator("#hcc1395-technical-validation").evaluate("element => element.getBoundingClientRect().top")
        run["screenshots"]["hcc_validation"] = screenshot(page, output_dir / f"hcc_validation__{viewport_name}.png")
        add_check(
            run,
            "screenshot_scroll_alignment",
            abs(comparison_top - 12) <= 2 and abs(hcc_top - 12) <= 2,
            {"comparison_top": 12, "hcc_top": 12, "tolerance": 2},
            {"comparison_top": comparison_top, "hcc_top": hcc_top},
        )
        add_check(run, "console_and_page_errors", not run["console_errors"] and not run["page_errors"], {"console_errors": 0, "page_errors": 0}, {"console_errors": run["console_errors"], "page_errors": run["page_errors"]})
        run["dimension_receipts"] = dimension_receipts
    except Exception as exc:
        add_check(run, "audit_exception", False, "no exception", str(exc))
    finally:
        run["duration_seconds"] = round(time.monotonic() - started, 3)
        context.close()
    return run


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    parser.add_argument("--comparison", type=Path, default=DEFAULT_COMPARISON)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout-ms", type=int, default=120_000)
    args = parser.parse_args()

    index_path = args.index.resolve()
    comparison_path = args.comparison.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    comparison = json.loads(comparison_path.read_text(encoding="utf-8"))
    comparison_sha256 = sha256_file(comparison_path)
    report: dict[str, Any] = {
        "schema_name": "intersubmod.sample_topology_comparison_playwright_audit",
        "schema_version": "1.0.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "index": str(index_path),
            "index_sha256": sha256_file(index_path),
            "comparison": str(comparison_path),
            "comparison_sha256": comparison_sha256,
        },
        "runs": [],
    }
    started = time.monotonic()
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True, args=["--allow-file-access-from-files"])
        report["browser"] = {"engine": "chromium", "version": browser.version, "headless": True}
        for viewport_name in VIEWPORTS:
            print(f"[RUN] {viewport_name}", flush=True)
            run = audit_viewport(browser, index_path, comparison, comparison_sha256, output_dir, viewport_name, args.timeout_ms)
            report["runs"].append(run)
            failed = [check["name"] for check in run["checks"] if check["status"] != "pass"]
            print(f"[{'PASS' if not failed else 'FAIL'}] {viewport_name} failures={','.join(failed) if failed else 'none'}", flush=True)
        browser.close()

    checks = [check for run in report["runs"] for check in run["checks"]]
    failures = [check for check in checks if check["status"] != "pass"]
    screenshots = [path for run in report["runs"] for path in run["screenshots"].values()]
    report["summary"] = {
        "status": "PASS" if not failures else "FAIL",
        "runs": len(report["runs"]),
        "checks_passed": len(checks) - len(failures),
        "checks_total": len(checks),
        "screenshots": len(screenshots),
        "elapsed_seconds": round(time.monotonic() - started, 3),
    }
    receipt_path = output_dir / "validation_receipt.json"
    receipt_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        f"RESULT status={report['summary']['status']} runs={report['summary']['runs']} "
        f"checks={report['summary']['checks_passed']}/{report['summary']['checks_total']} "
        f"screenshots={report['summary']['screenshots']}",
        flush=True,
    )
    print(f"RECEIPT {receipt_path}", flush=True)
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
