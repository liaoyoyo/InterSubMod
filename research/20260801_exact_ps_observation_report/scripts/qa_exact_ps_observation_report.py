#!/usr/bin/env python3
"""Browser, no-JS, responsive and print QA for the Exact-PS observation report."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import sys
from typing import Any

try:
    from playwright.sync_api import Browser, Page, sync_playwright
except ImportError as exc:  # pragma: no cover - exercised by deployment environment
    raise SystemExit(f"Playwright is required: {exc}")


SCHEMA_NAME = "intersubmod.exact_ps_observation_report.browser_qa"
SCHEMA_VERSION = "1.0.0"
REPORT_DATA_SCHEMA = "intersubmod.exact_ps_observation_report.data"
REPORT_DATA_VERSION = "1.0.0"
EXPECTED_SAMPLES = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
VIEWPORTS = (
    ("desktop", 1440, 1000),
    ("laptop", 1024, 900),
    ("mobile", 390, 844),
    ("narrow", 320, 700),
)
SENTINELS = (
    "98,955",
    "10,717",
    "71,955",
    "39,648",
    "8,449",
    "NOT_INTEGRATED",
    "association-only",
)


class QaError(RuntimeError):
    """A release-blocking browser or document-contract failure."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise QaError(message)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing output: {path}")
    require(not path.is_symlink(), f"symlink is not allowed: {path}")
    resolved = path.resolve(strict=True)
    require(stat.S_ISREG(resolved.stat().st_mode), f"not a regular file: {path}")
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_path(resolved),
    }


def read_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing JSON: {path}")
    try:
        payload = json.loads(
            path.read_text(encoding="utf-8"),
            parse_constant=lambda value: (_ for _ in ()).throw(
                ValueError(f"non-finite JSON constant: {value}")
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise QaError(f"cannot read JSON {path}: {exc}") from exc
    require(isinstance(payload, dict), f"JSON root is not an object: {path}")
    return payload


def validate_report_data(report_data: dict[str, Any]) -> None:
    require(
        report_data.get("schema_name") == REPORT_DATA_SCHEMA,
        "report-data schema mismatch",
    )
    require(
        report_data.get("schema_version") == REPORT_DATA_VERSION,
        "report-data version mismatch",
    )
    require(
        report_data.get("report_status") == "validated-derived-observation",
        "report-data status mismatch",
    )
    scope = report_data.get("scope")
    require(isinstance(scope, dict), "report-data scope missing")
    require(
        scope.get("technical_datasets") == list(EXPECTED_SAMPLES),
        "report-data technical dataset scope mismatch",
    )
    require(scope.get("technical_dataset_count") == 7, "technical dataset count drift")
    require(scope.get("biological_id_count") == 6, "biological ID count drift")
    require(scope.get("chromosome_scope") == "chr1–22", "chromosome scope drift")
    rows = report_data.get("per_sample")
    require(isinstance(rows, list) and len(rows) == 7, "per-sample rows missing")
    require(
        [row.get("sample") for row in rows] == list(EXPECTED_SAMPLES),
        "per-sample dataset identity/order drift",
    )
    require(
        report_data.get("cn_loh_readiness", {}).get("status") == "NOT_INTEGRATED",
        "CN/LOH status drift",
    )
    require(
        report_data.get("cn_loh_readiness", {}).get("availability") is None,
        "CN/LOH unknown must be null",
    )
    require(
        report_data.get("methyl_auxiliary", {}).get("status") == "association-only",
        "methyl role drift",
    )
    require(
        report_data.get("methyl_auxiliary", {}).get("topology_rescoring") is False,
        "methyl topology-rescoring claim drift",
    )
    require(
        report_data.get("tree_decision", {}).get("materialization_status")
        == "POLICY_ONLY_NOT_MATERIALIZED",
        "tree-decision materialization claim drift",
    )
    checks = report_data.get("conservation_checks")
    require(
        isinstance(checks, list)
        and checks
        and all(check.get("status") == "PASS" for check in checks),
        "report-data conservation checks failed",
    )
    provenance = report_data.get("provenance", {})
    require(provenance.get("authority_artifact_count") == 13, "authority count drift")
    require(provenance.get("strict_nested_bundle_count") == 9, "nested bundle count drift")


def write_json_exclusive(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(
        payload,
        ensure_ascii=False,
        allow_nan=False,
        indent=2,
        sort_keys=True,
    ) + "\n"
    try:
        with path.open("x", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
    except FileExistsError as exc:
        raise QaError(f"refusing to overwrite QA receipt: {path}") from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--html", type=Path, required=True)
    parser.add_argument("--report-data", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--chromium",
        type=Path,
        help="Optional Chromium executable; Playwright bundled Chromium is used otherwise.",
    )
    return parser.parse_args()


def attach_error_collectors(
    page: Page,
    console_errors: list[str],
    page_errors: list[str],
    external_requests: list[str],
    allowed_document_uri: str,
) -> None:
    page.on(
        "console",
        lambda message: (
            console_errors.append(message.text) if message.type == "error" else None
        ),
    )
    page.on("pageerror", lambda error: page_errors.append(str(error)))
    page.on(
        "request",
        lambda request: (
            external_requests.append(request.url)
            if request.url != allowed_document_uri
            and not request.url.startswith(("data:", "blob:", "about:"))
            else None
        ),
    )


def assert_common_document_contract(
    page: Page, report_data: dict[str, Any]
) -> dict[str, Any]:
    page.wait_for_load_state("domcontentloaded")
    require(page.locator("html").get_attribute("lang") == "zh-Hant", "html lang drift")
    require(page.locator("h1").count() == 1, "document must contain exactly one h1")
    require(
        page.locator(
            'body[data-report-status="validated-derived-observation"]'
        ).count()
        == 1,
        "report status marker is missing",
    )
    for selector in (
        "#overview",
        "#funnel",
        "#topology",
        "#samples",
        "#cn-loh",
        "#methyl",
        "#provenance",
    ):
        section = page.locator(selector)
        require(section.count() == 1, f"missing section: {selector}")
        require(section.is_visible(), f"section is not visible: {selector}")
        require(len(section.inner_text().strip()) >= 30, f"section is empty: {selector}")

    body_text = page.locator("body").inner_text()
    for sentinel in SENTINELS:
        require(sentinel in body_text, f"missing reader sentinel: {sentinel}")
    sample_rows = page.locator("#sample-table tbody tr[data-sample]")
    require(sample_rows.count() == 7, "sample table must contain seven datasets")
    observed_samples = sample_rows.evaluate_all(
        "rows => rows.map(row => row.dataset.sample)"
    )
    require(
        observed_samples == list(EXPECTED_SAMPLES),
        f"sample table identity/order drift: {observed_samples}",
    )
    embedded_text = page.locator("#report-data").text_content()
    require(embedded_text is not None, "embedded report-data is missing")
    try:
        embedded_data = json.loads(embedded_text)
    except json.JSONDecodeError as exc:
        raise QaError(f"embedded report-data is malformed: {exc}") from exc
    require(embedded_data == report_data, "HTML embedded data does not match report_data")

    headline = report_data["headline"]
    expected_metrics = {
        "final_groups": headline["final_groups"],
        "ranked": headline["ranked"],
        "one_rooted_unlabeled_topology": headline[
            "one_rooted_unlabeled_topology"
        ],
        "resource_abstain": headline["resource_abstain"],
    }
    for metric, expected in expected_metrics.items():
        target = page.locator(f'[data-metric="{metric}"]')
        require(target.count() == 1, f"missing metric binding: {metric}")
        require(
            target.get_attribute("data-value") == str(expected),
            f"metric DOM/data mismatch: {metric}",
        )

    sample_data = {row["sample"]: row for row in report_data["per_sample"]}
    field_paths = {
        "strict_read_linked_w": ("input", "strict_read_linked_w"),
        "final_groups": ("solver", "final_groups"),
        "resource_abstain": ("solver", "resource_abstain"),
        "ranked": ("solver", "ranked"),
        "unique_tree": ("topology", "resolution", "UNIQUE_TREE"),
        "tied_same_topology": (
            "topology",
            "resolution",
            "TIED_SAME_TOPOLOGY",
        ),
        "tied_cross_topology": (
            "topology",
            "resolution",
            "TIED_CROSS_TOPOLOGY",
        ),
        "one_exact_topology": ("topology", "one_exact_topology"),
        "methyl_formal": ("methyl", "formal_units"),
        "methyl_robust": ("methyl", "robust_associations"),
    }
    for sample in EXPECTED_SAMPLES:
        row_locator = page.locator(
            f'#sample-table tbody tr[data-sample="{sample}"]'
        )
        source = sample_data[sample]
        for field, path in field_paths.items():
            expected: Any = source
            for key in path:
                expected = expected[key]
            cell = row_locator.locator(f'[data-field="{field}"]')
            require(cell.count() == 1, f"{sample} missing bound field: {field}")
            require(
                cell.get_attribute("data-value") == str(expected),
                f"{sample} DOM/data mismatch: {field}",
            )
    require(page.locator("svg.chart-svg").count() >= 3, "fewer than three charts")
    for index in range(page.locator("svg.chart-svg").count()):
        chart = page.locator("svg.chart-svg").nth(index)
        require(
            chart.locator(":scope > title").count() == 1,
            f"chart {index} lacks one direct title",
        )
        require(
            chart.locator(":scope > desc").count() == 1,
            f"chart {index} lacks one direct desc",
        )
    require(page.locator("table caption").count() >= 1, "tables lack captions")
    require(
        page.locator("th[scope]").count() >= 5,
        "table headers lack explicit scope attributes",
    )
    overflow = page.evaluate(
        """() => ({
          scrollWidth: document.documentElement.scrollWidth,
          innerWidth: window.innerWidth,
          bodyWidth: document.body.getBoundingClientRect().width
        })"""
    )
    require(
        int(overflow["scrollWidth"]) <= int(overflow["innerWidth"]) + 2,
        f"document horizontal overflow: {overflow}",
    )
    return {
        "title": page.title(),
        "sample_rows": sample_rows.count(),
        "sample_ids": observed_samples,
        "chart_count": page.locator("svg.chart-svg").count(),
        "overflow": overflow,
        "embedded_data_equal": True,
        "metric_bindings_pass": True,
        "sample_bindings_pass": True,
    }


def exercise_javascript(page: Page) -> dict[str, Any]:
    button = page.locator('[data-sample-filter="H2009"]')
    require(button.count() == 1, "H2009 sample filter is missing")
    button.click()
    row = page.locator('#sample-table tbody tr[data-sample="H2009"]')
    require(row.count() == 1, "H2009 sample row is missing")
    active = row.evaluate(
        """element => (
          element.classList.contains('is-focus') ||
          element.getAttribute('data-highlighted') === 'true'
        )"""
    )
    require(bool(active), "sample filter did not focus H2009")
    require(
        button.get_attribute("aria-pressed") == "true",
        "sample filter aria state did not update",
    )
    return {"sample_filter": "H2009", "focused": True}


def run_viewport(
    browser: Browser,
    html_uri: str,
    output_dir: Path,
    name: str,
    width: int,
    height: int,
    report_data: dict[str, Any],
    *,
    javascript_enabled: bool = True,
) -> tuple[dict[str, Any], Path, list[str], list[str], list[str]]:
    context = browser.new_context(
        viewport={"width": width, "height": height},
        java_script_enabled=javascript_enabled,
        color_scheme="light",
        reduced_motion="reduce",
    )
    page = context.new_page()
    console_errors: list[str] = []
    page_errors: list[str] = []
    external_requests: list[str] = []
    attach_error_collectors(
        page,
        console_errors,
        page_errors,
        external_requests,
        html_uri,
    )
    page.goto(html_uri, wait_until="domcontentloaded")
    observed = assert_common_document_contract(page, report_data)
    if javascript_enabled and name == "desktop":
        observed["interaction"] = exercise_javascript(page)
    screenshot = output_dir / f"{name}.png"
    page.screenshot(path=str(screenshot), full_page=True)
    context.close()
    return observed, screenshot, console_errors, page_errors, external_requests


def render_print_pdf(
    browser: Browser,
    html_uri: str,
    output_dir: Path,
    report_data: dict[str, Any],
) -> tuple[Path, dict[str, Any], list[str], list[str], list[str]]:
    context = browser.new_context(
        viewport={"width": 1240, "height": 900},
        color_scheme="light",
        reduced_motion="reduce",
    )
    page = context.new_page()
    console_errors: list[str] = []
    page_errors: list[str] = []
    external_requests: list[str] = []
    attach_error_collectors(
        page,
        console_errors,
        page_errors,
        external_requests,
        html_uri,
    )
    page.goto(html_uri, wait_until="domcontentloaded")
    page.emulate_media(media="print")
    observed = assert_common_document_contract(page, report_data)
    page.locator("details").evaluate_all(
        "elements => elements.forEach(element => element.open = true)"
    )
    print_overflow = page.evaluate(
        """() => ({
          scrollWidth: document.documentElement.scrollWidth,
          innerWidth: window.innerWidth
        })"""
    )
    require(
        int(print_overflow["scrollWidth"]) <= int(print_overflow["innerWidth"]) + 2,
        f"print-media horizontal overflow: {print_overflow}",
    )
    pdf = output_dir / "report_A4.pdf"
    page.pdf(
        path=str(pdf),
        format="A4",
        print_background=True,
        margin={"top": "12mm", "right": "10mm", "bottom": "12mm", "left": "10mm"},
    )
    require(pdf.is_file() and pdf.stat().st_size > 10_000, "print PDF is empty")
    pdf_bytes = pdf.read_bytes()
    require(pdf_bytes.startswith(b"%PDF-"), "print output lacks PDF header")
    page_count = len(re.findall(rb"/Type\s*/Page\b", pdf_bytes))
    require(page_count >= 1, "print PDF contains no pages")
    context.close()
    observed["print_overflow"] = print_overflow
    observed["pdf_page_count"] = page_count
    return (
        pdf,
        observed,
        console_errors,
        page_errors,
        external_requests,
    )


def main() -> int:
    args = parse_args()
    html = args.html.expanduser().resolve()
    report_data_path = args.report_data.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    receipt_path = output_dir / "browser_qa_receipt.json"
    receipt: dict[str, Any] = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "inputs": {},
        "viewports": {},
        "no_js": {},
        "print": {},
        "artifacts": {},
        "checks": {},
        "all_pass": False,
    }
    try:
        require(html.is_file(), f"missing HTML: {html}")
        report_data = read_json(report_data_path)
        validate_report_data(report_data)
        require(
            not output_dir.exists() or not any(output_dir.iterdir()),
            f"QA output directory is not empty: {output_dir}",
        )
        require(not output_dir.is_symlink(), f"QA output directory is a symlink: {output_dir}")
        output_dir.mkdir(parents=True, exist_ok=True)
        receipt["inputs"] = {
            "html": identity(html),
            "report_data": identity(report_data_path),
        }
        launch_args: dict[str, Any] = {"headless": True}
        if args.chromium is not None:
            executable = args.chromium.expanduser().resolve()
            require(executable.is_file(), f"missing Chromium executable: {executable}")
            launch_args["executable_path"] = str(executable)

        all_console_errors: list[str] = []
        all_page_errors: list[str] = []
        all_external_requests: list[str] = []
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(**launch_args)
            for name, width, height in VIEWPORTS:
                observed, screenshot, console, page_errors, external = run_viewport(
                    browser,
                    html.as_uri(),
                    output_dir,
                    name,
                    width,
                    height,
                    report_data,
                )
                receipt["viewports"][name] = {
                    "width": width,
                    "height": height,
                    "observed": observed,
                    "screenshot": identity(screenshot),
                }
                all_console_errors.extend(console)
                all_page_errors.extend(page_errors)
                all_external_requests.extend(external)

            observed, screenshot, console, page_errors, external = run_viewport(
                browser,
                html.as_uri(),
                output_dir,
                "no_js",
                1024,
                900,
                report_data,
                javascript_enabled=False,
            )
            receipt["no_js"] = {
                "observed": observed,
                "screenshot": identity(screenshot),
            }
            all_console_errors.extend(console)
            all_page_errors.extend(page_errors)
            all_external_requests.extend(external)
            pdf, print_observed, console, page_errors, external = render_print_pdf(
                browser,
                html.as_uri(),
                output_dir,
                report_data,
            )
            receipt["print"] = {
                "pdf": identity(pdf),
                "observed": print_observed,
            }
            all_console_errors.extend(console)
            all_page_errors.extend(page_errors)
            all_external_requests.extend(external)
            browser.close()

        receipt["checks"] = {
            "report_data_schema_pass": True,
            "report_data_contract_pass": True,
            "html_data_binding_pass": all(
                viewport["observed"].get("embedded_data_equal") is True
                for viewport in receipt["viewports"].values()
            ),
            "sample_identity_and_binding_pass": all(
                viewport["observed"].get("sample_bindings_pass") is True
                for viewport in receipt["viewports"].values()
            ),
            "four_viewports_pass": len(receipt["viewports"]) == 4,
            "no_js_core_content_pass": bool(receipt["no_js"]),
            "print_pdf_pass": bool(receipt["print"]),
            "print_page_and_overflow_pass": bool(
                receipt["print"].get("observed", {}).get("pdf_page_count")
            ),
            "console_errors_zero": not all_console_errors,
            "page_errors_zero": not all_page_errors,
            "external_requests_zero": not all_external_requests,
            "sample_interaction_pass": bool(
                receipt["viewports"]["desktop"]["observed"].get("interaction")
            ),
        }
        receipt["diagnostics"] = {
            "console_errors": all_console_errors,
            "page_errors": all_page_errors,
            "external_requests": sorted(set(all_external_requests)),
        }
        receipt["all_pass"] = all(receipt["checks"].values())
        require(receipt["all_pass"], f"browser QA checks failed: {receipt['checks']}")
    except Exception as exc:  # noqa: BLE001 - receipt must preserve the blocking cause
        receipt["failure"] = f"{type(exc).__name__}: {exc}"
        receipt["all_pass"] = False
    try:
        write_json_exclusive(receipt_path, receipt)
    except QaError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2
    print(f"[INPUT] html={html}")
    print(f"[INPUT] report_data={report_data_path}")
    print(f"[OUTPUT] receipt={receipt_path}")
    print(
        json.dumps(
            {"all_pass": receipt["all_pass"], "checks": receipt["checks"]},
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if receipt["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
