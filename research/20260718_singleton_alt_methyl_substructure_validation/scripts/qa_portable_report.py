#!/usr/bin/env python3
"""Capture deterministic desktop/mobile QA evidence for the portable report."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

from playwright.sync_api import sync_playwright


TOPIC_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CHROMIUM = Path(
    "/bip7_disk/liaoyoyo2001/.cache/ms-playwright/chromium-1223/chrome-linux64/chrome"
)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_json(path: Path, payload: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def ensure_topic_path(path: Path) -> None:
    try:
        path.resolve().relative_to(TOPIC_ROOT.resolve())
    except ValueError as exc:
        raise ValueError(f"QA output must remain within {TOPIC_ROOT}") from exc


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", type=Path, required=True)
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--chromium", type=Path, default=DEFAULT_CHROMIUM)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    ensure_topic_path(args.output_dir)
    if args.output_dir.exists():
        raise FileExistsError(args.output_dir)
    args.output_dir.mkdir(parents=True)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    expected = {
        "blocks": len(artifact["manifest"]["blocks"]),
        "charts": len(artifact["manifest"].get("charts", [])),
        "html_blocks": sum(block["type"] == "html" for block in artifact["manifest"]["blocks"]),
        "tables": len(artifact["manifest"].get("tables", [])),
    }
    console_errors: list[str] = []
    page_errors: list[str] = []
    failed_requests: list[str] = []
    results: list[dict[str, Any]] = []
    receipt_path = args.output_dir / "browser_qa_receipt.json"
    try:
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(
                headless=True,
                executable_path=str(args.chromium.resolve()),
            )
            for name, width, height in (("desktop", 1440, 1000), ("mobile", 390, 844)):
                context = browser.new_context(
                    viewport={"width": width, "height": height},
                    color_scheme="light",
                    device_scale_factor=1,
                )
                page = context.new_page()
                page.on(
                    "console",
                    lambda message, viewport=name: (
                        console_errors.append(f"{viewport}: {message.text}")
                        if message.type == "error"
                        else None
                    ),
                )
                page.on("pageerror", lambda error, viewport=name: page_errors.append(f"{viewport}: {error}"))
                page.on(
                    "requestfailed",
                    lambda request, viewport=name: failed_requests.append(
                        f"{viewport}: {request.url} {request.failure}"
                    ),
                )
                page.goto(args.html.resolve().as_uri(), wait_until="load", timeout=60_000)
                page.wait_for_selector("#data-analytics-portable-reader", state="visible", timeout=30_000)
                page.wait_for_function(
                    """() => {
                        const fallback = document.querySelector('[data-portable-fallback]');
                        const reader = document.querySelector('#data-analytics-portable-reader');
                        return reader && reader.getBoundingClientRect().height > 0 &&
                               fallback && getComputedStyle(fallback).display === 'none';
                    }""",
                    timeout=30_000,
                )
                page.wait_for_timeout(500)
                text = page.locator("body").inner_text()
                for token in (
                    "50,432",
                    "8,279",
                    "observed operational M2-PASS yield",
                    "confirmed cellular subclone=0",
                    "linear ancestry=0",
                    "PS 不在 M2 measured axes",
                ):
                    if token not in text:
                        raise AssertionError(f"{name}: missing visible token {token!r}")
                geometry = page.evaluate(
                    """() => {
                        const root = document.documentElement;
                        const reader = document.querySelector('#data-analytics-portable-reader');
                        const charts = [...document.querySelectorAll(
                            '#data-analytics-portable-reader section[data-artifact-kind="chart"]'
                        )];
                        const tables = [...document.querySelectorAll(
                            '#data-analytics-portable-reader section[data-artifact-kind="table"]'
                        )];
                        const frames = [...document.querySelectorAll('iframe.report-html-frame')];
                        return {
                            clientWidth: root.clientWidth,
                            scrollWidth: root.scrollWidth,
                            readerWidth: reader?.getBoundingClientRect().width || 0,
                            readerHeight: reader?.getBoundingClientRect().height || 0,
                            chartCount: charts.length,
                            tableCount: tables.length,
                            htmlFrameCount: frames.length,
                            zeroGeometryCharts: charts.filter(node => {
                                const rect = node.getBoundingClientRect();
                                return rect.width <= 0 || rect.height <= 0;
                            }).length,
                            zeroGeometryTables: tables.filter(node => {
                                const rect = node.getBoundingClientRect();
                                return rect.width <= 0 || rect.height <= 0;
                            }).length,
                        };
                    }"""
                )
                if geometry["scrollWidth"] > geometry["clientWidth"] + 1:
                    raise AssertionError(f"{name}: horizontal overflow {geometry}")
                if geometry["readerWidth"] <= 0 or geometry["readerHeight"] <= 0:
                    raise AssertionError(f"{name}: reader has zero geometry")
                if geometry["chartCount"] != expected["charts"]:
                    raise AssertionError(f"{name}: chart count {geometry['chartCount']} != {expected['charts']}")
                if geometry["tableCount"] != expected["tables"]:
                    raise AssertionError(f"{name}: table count {geometry['tableCount']} != {expected['tables']}")
                if geometry["htmlFrameCount"] != expected["html_blocks"]:
                    raise AssertionError(
                        f"{name}: HTML frame count {geometry['htmlFrameCount']} != {expected['html_blocks']}"
                    )
                if geometry["zeroGeometryCharts"] or geometry["zeroGeometryTables"]:
                    raise AssertionError(f"{name}: zero-geometry visual {geometry}")
                frame_checks = []
                for frame in page.frames:
                    if frame == page.main_frame:
                        continue
                    owner = frame.frame_element()
                    owner_class = owner.get_attribute("class") or ""
                    if "report-html-frame" not in owner_class.split():
                        continue
                    image_count = frame.locator("img").count()
                    if image_count:
                        image_data = frame.locator("img").first.evaluate(
                            """img => ({
                                complete: img.complete,
                                naturalWidth: img.naturalWidth,
                                naturalHeight: img.naturalHeight,
                                width: img.getBoundingClientRect().width,
                                height: img.getBoundingClientRect().height,
                            })"""
                        )
                        if not image_data["complete"] or image_data["naturalWidth"] < 100 or image_data["naturalHeight"] < 100:
                            raise AssertionError(f"{name}: embedded methyl heatmap failed {image_data}")
                        frame_checks.append(image_data)
                if len(frame_checks) != expected["html_blocks"]:
                    raise AssertionError(
                        f"{name}: embedded image frames {len(frame_checks)} != {expected['html_blocks']}"
                    )
                screenshot_path = args.output_dir / f"{name}.png"
                page.screenshot(path=str(screenshot_path), full_page=True)
                results.append(
                    {
                        "viewport": name,
                        "width": width,
                        "height": height,
                        "geometry": geometry,
                        "embedded_images": frame_checks,
                        "screenshot": str(screenshot_path.resolve()),
                    }
                )
                context.close()
            browser.close()
        if console_errors or page_errors or failed_requests:
            raise AssertionError(
                f"browser errors: console={console_errors}, page={page_errors}, requests={failed_requests}"
            )
        receipt = {
            "schema_name": "intersubmod.singleton_sidecar_browser_qa_receipt",
            "schema_version": "1.0.0",
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "ok": True,
            "html": {
                "path": str(args.html.resolve()),
                "sha256": sha256_path(args.html),
            },
            "artifact": {
                "path": str(args.artifact.resolve()),
                "sha256": sha256_path(args.artifact),
            },
            "chromium": str(args.chromium.resolve()),
            "expected": expected,
            "results": results,
            "console_errors": console_errors,
            "page_errors": page_errors,
            "failed_requests": failed_requests,
        }
        safe_json(receipt_path, receipt)
        for path in args.output_dir.rglob("*"):
            if path.is_file():
                path.chmod(0o444)
        print(json.dumps(receipt, ensure_ascii=False, sort_keys=True))
        return 0
    except Exception as exc:
        if not receipt_path.exists():
            safe_json(
                receipt_path,
                {
                    "schema_name": "intersubmod.singleton_sidecar_browser_qa_receipt",
                    "schema_version": "1.0.0",
                    "created_at_utc": datetime.now(timezone.utc).isoformat(),
                    "ok": False,
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "traceback": traceback.format_exc(),
                    "partial_results": results,
                    "console_errors": console_errors,
                    "page_errors": page_errors,
                    "failed_requests": failed_requests,
                },
            )
        for path in args.output_dir.rglob("*"):
            if path.is_file():
                path.chmod(0o444)
        print(json.dumps({"ok": False, "receipt": str(receipt_path), "error": str(exc)}), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
