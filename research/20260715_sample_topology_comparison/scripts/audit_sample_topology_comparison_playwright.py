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
    "narrow": {"width": 320, "height": 800},
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
          const navHeight = document.querySelector('.local-nav')?.getBoundingClientRect().height || 0;
          window.scrollTo(0, Math.max(0, element.getBoundingClientRect().top + window.scrollY - navHeight - 12));
        }""",
        selector,
    )
    page.wait_for_function(
        "selector => { const element=document.querySelector(selector); const navHeight=document.querySelector('.local-nav')?.getBoundingClientRect().height || 0; return element && Math.abs(element.getBoundingClientRect().top - navHeight - 12) <= 2; }",
        arg=selector,
    )
    page.wait_for_timeout(80)


def screenshot(page: Page, path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    page.screenshot(path=str(path), full_page=False)
    return str(path.resolve())


def text_contrast_failures(page: Page) -> list[dict[str, Any]]:
    """Return visible direct-text nodes below WCAG 2.x AA contrast thresholds."""

    return page.evaluate(
        r"""() => {
          const parse = value => {
            const match = value && value.match(/rgba?\(([-.\d]+)[, ]+([-.\d]+)[, ]+([-.\d]+)(?:[, /]+([-.\d]+))?\)/);
            return match
              ? [Number(match[1]), Number(match[2]), Number(match[3]), match[4] == null ? 1 : Number(match[4])]
              : [0, 0, 0, 0];
          };
          const blend = (foreground, background) => {
            const alpha = foreground[3] + background[3] * (1 - foreground[3]);
            if (!alpha) return [255, 255, 255, 1];
            return [
              (foreground[0] * foreground[3] + background[0] * background[3] * (1 - foreground[3])) / alpha,
              (foreground[1] * foreground[3] + background[1] * background[3] * (1 - foreground[3])) / alpha,
              (foreground[2] * foreground[3] + background[2] * background[3] * (1 - foreground[3])) / alpha,
              alpha
            ];
          };
          const backgroundFor = element => {
            const layers = [];
            for (let current = element; current; current = current.parentElement) {
              const color = parse(getComputedStyle(current).backgroundColor);
              if (color[3]) layers.push(color);
            }
            return layers.reverse().reduce((background, layer) => blend(layer, background), [255, 255, 255, 1]);
          };
          const luminance = color => {
            const channels = color.slice(0, 3).map(value => {
              const normalized = value / 255;
              return normalized <= 0.04045
                ? normalized / 12.92
                : ((normalized + 0.055) / 1.055) ** 2.4;
            });
            return 0.2126 * channels[0] + 0.7152 * channels[1] + 0.0722 * channels[2];
          };
          const failures = [];
          for (const element of document.querySelectorAll('body *')) {
            const text = [...element.childNodes]
              .filter(node => node.nodeType === Node.TEXT_NODE)
              .map(node => node.textContent.trim())
              .filter(Boolean)
              .join(' ');
            if (!text || element.closest('.sr-only')) continue;
            const style = getComputedStyle(element);
            const rect = element.getBoundingClientRect();
            if (style.display === 'none' || style.visibility === 'hidden' || Number(style.opacity) < 0.1 || rect.width < 1 || rect.height < 1) continue;
            const foreground = parse(style.color);
            const background = backgroundFor(element);
            const foregroundLuminance = luminance(foreground);
            const backgroundLuminance = luminance(background);
            const ratio = (Math.max(foregroundLuminance, backgroundLuminance) + 0.05)
              / (Math.min(foregroundLuminance, backgroundLuminance) + 0.05);
            const fontSize = Number.parseFloat(style.fontSize);
            const fontWeight = Number.parseInt(style.fontWeight, 10) || 400;
            const large = fontSize >= 24 || (fontSize >= 18.66 && fontWeight >= 700);
            const required = large ? 3 : 4.5;
            if (ratio + 1e-6 < required) {
              failures.push({
                ratio: Number(ratio.toFixed(2)),
                required,
                text: text.slice(0, 100),
                tag: element.tagName.toLowerCase(),
                class_name: String(element.className || '').slice(0, 100)
              });
            }
          }
          return failures.slice(0, 40);
        }"""
    )


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
          const confusionValues = [...(confusion?.querySelectorAll('tbody td:not(.confusion-total)') || [])]
            .map(cell => Number((cell.textContent || '').replaceAll(',', '')));
          const ledger = document.querySelector('.comparison-table');
          const segmentSums = profileRows.map(row => [...row.querySelectorAll('.bar-segment')]
            .reduce((sum, segment) => sum + Number.parseFloat(segment.style.width || '0'), 0));
          return {
            embedded_dataset_count: embedded.datasets.length,
            embedded_biological_ids: [...new Set(embedded.datasets.map(row => row.biological_id))].sort(),
            active_dimension: active?.dataset.comparisonDimension || null,
            active_tab_count: document.querySelectorAll('.comparison-tab[aria-pressed="true"]').length,
            active_profile_scope: document.querySelector('[data-profile-scope][aria-pressed="true"]')?.dataset.profileScope || null,
            active_distance_scope: document.querySelector('[data-distance-scope][aria-pressed="true"]')?.dataset.distanceScope || null,
            active_confusion_mode: document.querySelector('[data-confusion-mode][aria-pressed="true"]')?.dataset.confusionMode || null,
            aggregate_profile_count: document.querySelectorAll('#aggregate-profile-chart .aggregate-profile').length,
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
            confusion_row_totals: [...(confusion?.querySelectorAll('tbody .confusion-total') || [])].map(cell => (cell.textContent || '').trim()),
            table_captions: document.querySelectorAll('table > caption').length,
            selected_matrix_cells: document.querySelectorAll('.matrix-cell[aria-pressed="true"]').length,
            comparison_title: document.getElementById('comparison-profile-title')?.textContent || '',
            matrix_title: document.getElementById('comparison-matrix-title')?.textContent || '',
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
                and state["active_profile_scope"] == "dataset"
                and state["active_distance_scope"] == "full"
                and state["active_confusion_mode"] == "count"
                and state["aggregate_profile_count"] == 2
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
                and state["confusion_columns"] == len(expected_categories) + 2
                and state["confusion_sum"] == expected_matched
                and state["selected_matrix_cells"] == 2
                and state["table_captions"] >= 5
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

            page.locator('[data-distance-scope="conditional"]').click()
            page.wait_for_timeout(80)
            conditional_state = dimension_state(page)
            expected_conditional_hcc = comparison["hcc1395_technical_pair"]["marginal_profiles"]["dimensions"][dimension]["conditional_evaluable_tvd"]
            expected_conditional: dict[tuple[str, str], float] = {}
            for record in comparison["pairwise_composition"]["records"]:
                expected_conditional[(record["left"], record["right"])] = record["dimensions"][dimension]["conditional_evaluable_tvd"]
                expected_conditional[(record["right"], record["left"])] = record["dimensions"][dimension]["conditional_evaluable_tvd"]
            conditional_values_match = all(
                (
                    abs(conditional_state["matrix_values"][i][j]) <= 1e-12
                    if left == right
                    else abs(conditional_state["matrix_values"][i][j] - expected_conditional[(left, right)]) <= 0.00051
                )
                for i, left in enumerate(order)
                for j, right in enumerate(order)
            )
            conditional_passed = (
                conditional_state["active_distance_scope"] == "conditional"
                and conditional_state["matrix_rows"] == 7
                and conditional_state["matrix_columns"] == 7
                and conditional_values_match
                and abs(conditional_state["hcc_matrix_value"] - expected_conditional_hcc) <= 0.00051
                and "Conditional evaluable" in conditional_state["matrix_title"]
            )
            add_check(
                run,
                f"conditional_matrix_{dimension}",
                conditional_passed,
                {"scope": "conditional", "pairs": 21, "hcc_tvd": expected_conditional_hcc},
                conditional_state,
            )
            page.locator('[data-distance-scope="full"]').click()
            page.wait_for_timeout(40)

        page.locator('.comparison-tab[data-comparison-dimension="structural"]').click()
        page.locator('[data-profile-scope="biological"]').click()
        page.wait_for_timeout(80)
        macro_state = page.evaluate(
            """() => ({
              active: document.querySelector('[data-profile-scope][aria-pressed="true"]')?.dataset.profileScope,
              ids: [...document.querySelectorAll('#comparison-profile-chart [data-biological-id]')].map(row => row.dataset.biologicalId),
              segment_sums: [...document.querySelectorAll('#comparison-profile-chart .comparison-profile-row')].map(row => [...row.querySelectorAll('.bar-segment')].reduce((sum, segment) => sum + Number.parseFloat(segment.style.width || '0'), 0)),
              invented_counts: [...document.querySelectorAll('#comparison-profile-chart .comparison-profile-n')].some(node => /W_primary/i.test(node.textContent || ''))
            })"""
        )
        add_check(
            run,
            "biological_macro_profiles",
            macro_state["active"] == "biological"
            and len(macro_state["ids"]) == 6
            and set(macro_state["ids"]) == {"COLO829", "H1437", "H2009", "HCC1395", "HCC1937", "HCC1954"}
            and all(abs(value - 100.0) <= 0.001 for value in macro_state["segment_sums"])
            and not macro_state["invented_counts"],
            {"biological_ids": 6, "equal_weight_profiles": True, "invented_counts": False},
            macro_state,
        )
        scroll_to(page, "#comparison-profile-chart")
        run["screenshots"]["biological_macro"] = screenshot(
            page, output_dir / f"biological_macro__{viewport_name}.png"
        )
        page.locator('[data-profile-scope="dataset"]').click()

        page.locator('[data-distance-scope="conditional"]').click()
        page.wait_for_timeout(80)
        scroll_to(page, "#comparison-matrix-title")
        run["screenshots"]["conditional_matrix"] = screenshot(
            page, output_dir / f"conditional_matrix__{viewport_name}.png"
        )
        page.locator('[data-distance-scope="full"]').click()

        page.locator('.hcc-deep-dive > summary').click()
        page.locator('[data-confusion-mode="row-percent"]').click()
        page.wait_for_timeout(80)
        row_percent_state = page.evaluate(
            """() => ({
              active: document.querySelector('[data-confusion-mode][aria-pressed="true"]')?.dataset.confusionMode,
              rows: [...document.querySelectorAll('.confusion-table tbody tr')].map(row => [...row.querySelectorAll('td:not(.confusion-total)')].map(cell => Number.parseFloat(cell.textContent || '0'))),
              totals: [...document.querySelectorAll('.confusion-table .confusion-total')].map(cell => (cell.textContent || '').trim()),
              first_cell: document.querySelector('.confusion-table tbody td:not(.confusion-total)')?.textContent?.trim()
            })"""
        )
        row_sums = [sum(row) for row in row_percent_state["rows"]]
        expected_first = (
            100
            * comparison["hcc1395_technical_pair"]["matched_primary_dimensions"]["structural"]["matrix"]["exact_and_topology_unique"]["exact_and_topology_unique"]
            / comparison["hcc1395_technical_pair"]["matched_primary_dimensions"]["structural"]["left_counts"]["exact_and_topology_unique"]
        )
        add_check(
            run,
            "confusion_row_percent",
            row_percent_state["active"] == "row-percent"
            and all(abs(total - 100.0) <= 0.11 for total in row_sums)
            and all(value == "100.0%" for value in row_percent_state["totals"])
            and abs(float(row_percent_state["first_cell"].rstrip("%")) - expected_first) <= 0.051,
            {"row_sums": "100.0% ±0.1", "first_cell": expected_first},
            {**row_percent_state, "row_sums": row_sums},
        )
        scroll_to(page, "#hcc-confusion-title")
        run["screenshots"]["confusion_row_percent"] = screenshot(
            page, output_dir / f"confusion_row_percent__{viewport_name}.png"
        )
        contrast_failures = text_contrast_failures(page)
        add_check(
            run,
            "visible_text_wcag_aa_contrast",
            not contrast_failures,
            {"failures": 0, "normal_text_ratio": 4.5, "large_text_ratio": 3.0},
            {"failures": len(contrast_failures), "examples": contrast_failures},
        )
        page.locator('[data-confusion-mode="count"]').click()
        page.locator('.hcc-deep-dive > summary').click()

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
              toggles: [...document.querySelectorAll('.inline-toggle')].map(element => element.getBoundingClientRect().height),
              nav: [...document.querySelectorAll('.local-nav a')].map(element => element.getBoundingClientRect().height),
              cells: [...document.querySelectorAll('.matrix-cell')].slice(0, 8).map(element => element.getBoundingClientRect().height)
            })"""
        )
        add_check(
            run,
            "interactive_targets_at_least_44px",
            all(value >= 44 for value in touch_sizes["tabs"] + touch_sizes["toggles"] + touch_sizes["nav"] + touch_sizes["cells"]),
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
                exact_edge_visible: (document.querySelector('#hcc1395-technical-validation')?.innerText || '').includes('35.3%'),
                allele_site_label_visible: (document.querySelector('#hcc1395-technical-validation')?.innerText || '').includes('Exact retained allele-site'),
                answer_before_hcc: document.getElementById('answer-strip')?.compareDocumentPosition(document.getElementById('hcc1395-technical-validation')) & Node.DOCUMENT_POSITION_FOLLOWING,
                hcc_before_comparison: document.getElementById('hcc1395-technical-validation')?.compareDocumentPosition(document.getElementById('cohort-comparison')) & Node.DOCUMENT_POSITION_FOLLOWING,
                local_nav_visible: document.querySelector('.local-nav')?.getBoundingClientRect().height >= 44,
                table_captions: document.querySelectorAll('table > caption').length,
                selected_pair_aria_pressed: document.querySelectorAll('.matrix-cell[aria-pressed="true"]').length,
                hcc_details_open: document.querySelector('.hcc-deep-dive')?.open === true
              };
            }"""
        )
        add_check(run, "body_and_internal_overflow", layout["body_overflow_px"] <= 1 and layout["matrix_scroll_width"] >= layout["matrix_client_width"], {"body_overflow_max": 1, "matrix_internal_scroll_allowed": True}, layout)
        add_check(run, "json_links_hidden_in_closed_details", layout["json_links"] > 0 and layout["json_links_inside_details"] == layout["json_links"] and layout["json_links_in_open_details"] == 0, "all JSON links inside closed details", layout)
        add_check(run, "hcc_exact_signature_visible", layout["hcc_metric_cards"] == 4 and layout["exact_edge_visible"], {"metric_cards": 4, "exact_edge": "35.3%"}, layout)
        add_check(
            run,
            "answer_first_navigation_and_semantics",
            bool(layout["answer_before_hcc"])
            and bool(layout["hcc_before_comparison"])
            and layout["local_nav_visible"]
            and layout["table_captions"] >= 5
            and layout["selected_pair_aria_pressed"] == 2
            and layout["allele_site_label_visible"]
            and not layout["hcc_details_open"],
            {"order": "answer → HCC → comparison", "sticky_nav": True, "captions": ">=5", "selected_pair_aria": 2, "allele_site_label": True},
            layout,
        )

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
        nav_height = page.locator(".local-nav").evaluate("element => element.getBoundingClientRect().height")
        expected_screenshot_top = nav_height + 12
        add_check(
            run,
            "screenshot_scroll_alignment",
            abs(comparison_top - expected_screenshot_top) <= 2
            and abs(hcc_top - expected_screenshot_top) <= 2,
            {
                "comparison_top": expected_screenshot_top,
                "hcc_top": expected_screenshot_top,
                "tolerance": 2,
            },
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
