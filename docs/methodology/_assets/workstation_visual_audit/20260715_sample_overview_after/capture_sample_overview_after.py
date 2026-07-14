#!/usr/bin/env python3
"""Capture and validate the canonical layered workstation overview redesign.

Inputs:
  docs/methodology/_assets/layered_workstation/index.html
  docs/methodology/_assets/layered_workstation/{7 datasets}.html

Outputs:
  this directory / metrics.json
  this directory / screenshots/*.png

The audit is intentionally independent from the renderer's embedded GRCh38
length table.  A mismatch in coordinate scaling, panel closure, runtime errors,
offline behavior, keyboard interaction, or body overflow makes the command
exit non-zero.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from playwright.sync_api import Page, sync_playwright


ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
WORKSTATION = ROOT / "docs/methodology/_assets/layered_workstation"
RENDERER = (
    ROOT
    / "docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py"
)
EXPECTED_RENDERER_SHA256 = hashlib.sha256(RENDERER.read_bytes()).hexdigest()
OUT = ROOT / "docs/methodology/_assets/workstation_visual_audit/20260715_sample_overview_after"
SCREENSHOTS = OUT / "screenshots"
SAMPLES = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)
VIEWPORTS = {
    "desktop1440": {"width": 1440, "height": 1000},
    "tablet768": {"width": 768, "height": 1024},
    "mobile390": {"width": 390, "height": 844},
    "narrow320": {"width": 320, "height": 720},
}
GRCH38 = {
    "chr1": 248_956_422,
    "chr2": 242_193_529,
    "chr3": 198_295_559,
    "chr4": 190_214_555,
    "chr5": 181_538_259,
    "chr6": 170_805_979,
    "chr7": 159_345_973,
    "chr8": 145_138_636,
    "chr9": 138_394_717,
    "chr10": 133_797_422,
    "chr11": 135_086_622,
    "chr12": 133_275_309,
    "chr13": 114_364_328,
    "chr14": 107_043_718,
    "chr15": 101_991_189,
    "chr16": 90_338_345,
    "chr17": 83_257_441,
    "chr18": 80_373_285,
    "chr19": 58_617_616,
    "chr20": 64_444_167,
    "chr21": 46_709_983,
    "chr22": 50_818_468,
}
MODES = ("determinacy", "evidence", "primary-hp", "n-ssnv", "cn-sidecar")
PANEL_IDS = ("topology-count", "candidate-count", "determinacy", "hp-h3", "region-size")


def capture_section(page: Page, selector: str, target: Path) -> None:
    locator = page.locator(selector)
    locator.scroll_into_view_if_needed()
    page.evaluate("document.activeElement && document.activeElement.blur()")
    page.wait_for_timeout(80)
    # Playwright locator screenshots rebase fixed-position elements to the
    # clipped section.  A correctly off-canvas skip link can therefore appear
    # inside the capture even though it is invisible in the real viewport.
    # Hide it only while recording the clean visual-reference state.
    page.evaluate(
        """() => document.querySelectorAll('.skip-link').forEach(link => {
            link.dataset.captureVisibility = link.style.visibility;
            link.style.visibility = 'hidden';
        })"""
    )
    try:
        locator.screenshot(path=str(target), animations="disabled")
    finally:
        page.evaluate(
            """() => document.querySelectorAll('.skip-link').forEach(link => {
                link.style.visibility = link.dataset.captureVisibility || '';
                delete link.dataset.captureVisibility;
            })"""
        )


def audit_sample(page: Page, sample: str, viewport: str) -> dict:
    page.wait_for_function(
        "document.querySelectorAll('.ideogram-mark').length === JSON.parse(document.getElementById('workstation-data').textContent).canonical_sample.W_tree",
        timeout=120_000,
    )
    metrics = page.evaluate(
        """({lengths, panelIds}) => {
            const data = JSON.parse(document.getElementById('workstation-data').textContent);
            const overviewById = Object.fromEntries(data.sample_overview.panels.map(panel => [panel.id, panel]));
            const panelMetrics = panelIds.map(id => {
                const panel = document.querySelector(`[data-overview-panel="${id}"]`);
                const expected = overviewById[id];
                const bins = [...panel.querySelectorAll('[data-overview-bin]')].map(bin => ({
                    key: bin.dataset.binKey,
                    count: Number(bin.dataset.count),
                    denominator: Number(bin.dataset.denominator),
                    width: parseFloat(bin.querySelector('[data-overview-bar]').style.width)
                }));
                return {
                    id,
                    denominator_key: panel.dataset.denominatorKey,
                    total: Number(panel.dataset.total),
                    sum: bins.reduce((sum, bin) => sum + bin.count, 0),
                    bins,
                    expected_bins: expected.bins.map(bin => ({key: bin.key, count: Number(bin.count)})),
                    visible_denominator: (panel.innerText || '').includes(panel.dataset.denominatorKey) &&
                        (panel.innerText || '').includes(Number(panel.dataset.total).toLocaleString())
                };
            });
            const tracks = [...document.querySelectorAll('.ideogram-track')];
            const marks = [...document.querySelectorAll('.ideogram-mark')];
            const tracksByChrom = Object.fromEntries(tracks.map(track => [track.dataset.chrom, track]));
            const regionRows = Object.fromEntries(data.region_index.map(row => [row.region, row]));
            let coordinateErrors = 0, identityErrors = 0, maxXError = 0, maxRatioError = 0;
            const seen = new Set();
            marks.forEach(mark => {
                const row = regionRows[mark.dataset.region];
                const chrom = mark.dataset.chrom;
                const start = Number(mark.dataset.start), end = Number(mark.dataset.end);
                const midpoint = Number(mark.dataset.midpoint), length = Number(lengths[chrom]);
                if (!row || row.chrom !== chrom || Number(row.start) !== start || Number(row.end) !== end) identityErrors++;
                if (seen.has(mark.dataset.region)) identityErrors++;
                seen.add(mark.dataset.region);
                if (!(1 <= start && start <= end && end <= length && start <= midpoint && midpoint <= end)) coordinateErrors++;
                const track = tracksByChrom[chrom];
                const expectedX = Number(track.getAttribute('x')) + midpoint / length * Number(track.getAttribute('width'));
                maxXError = Math.max(maxXError, Math.abs(Number(mark.getAttribute('x1')) - expectedX));
            });
            const chr1Width = Number(tracksByChrom.chr1.getAttribute('width'));
            tracks.forEach(track => {
                const expectedRatio = Number(lengths[track.dataset.chrom]) / Number(lengths.chr1);
                const observedRatio = Number(track.getAttribute('width')) / chr1Width;
                maxRatioError = Math.max(maxRatioError, Math.abs(expectedRatio - observedRatio));
            });
            const jsonLinks = [...document.querySelectorAll('a[href*=".json"]')];
            const scroller = document.querySelector('.ideogram-scroll');
            const buttons = [...document.querySelectorAll('[data-ideogram-mode]')];
            const bodyWidth = Math.max(document.documentElement.scrollWidth, document.body.scrollWidth);
            return {
                sample: data.sample,
                W_tree: Number(data.canonical_sample.W_tree),
                W_primary: Number(data.canonical_sample.W_primary),
                retained_sSNV: Number(data.canonical_sample.retained_sSNV),
                panel_metrics: panelMetrics,
                weighted_site_sum: data.region_index.reduce((sum, row) => sum + Number(row.n_sSNV), 0),
                region_index_count: data.region_index.length,
                tracks: tracks.length,
                marks: marks.length,
                unique_marks: seen.size,
                coordinate_errors: coordinateErrors,
                identity_errors: identityErrors,
                maximum_x_error: maxXError,
                maximum_track_ratio_error: maxRatioError,
                mode_buttons: buttons.length,
                mode_min_height: Math.min(...buttons.map(button => button.getBoundingClientRect().height)),
                body_overflow_px: Math.max(0, bodyWidth - innerWidth),
                json_links: jsonLinks.length,
                json_links_outside_closed_drawer: jsonLinks.filter(link => {
                    const details = link.closest('details');
                    return !details || details.open;
                }).length,
                scroller: {
                    role: scroller.getAttribute('role'),
                    tabindex: scroller.getAttribute('tabindex'),
                    aria_label: scroller.getAttribute('aria-label'),
                    overflow_x: getComputedStyle(scroller).overflowX,
                    scroll_width: scroller.scrollWidth,
                    client_width: scroller.clientWidth
                },
                h1_count: document.querySelectorAll('h1').length,
                overview_labelled: document.getElementById('sample-overview').getAttribute('aria-labelledby'),
                genome_labelled: document.getElementById('genome-overview').getAttribute('aria-labelledby'),
                svg_role: document.getElementById('ideogram-svg').getAttribute('role'),
                svg_title: Boolean(document.querySelector('#ideogram-svg title')),
                svg_desc: Boolean(document.querySelector('#ideogram-svg desc')),
                ui_contract: document.querySelector('meta[name="intersubmod-ui-contract"]')?.content || '',
                renderer_sha: document.querySelector('meta[name="intersubmod-renderer-sha256"]')?.content || ''
            };
        }""",
        {"lengths": GRCH38, "panelIds": PANEL_IDS},
    )
    expected_denominators = {
        "topology-count": ("W_primary", metrics["W_primary"]),
        "candidate-count": ("W_primary", metrics["W_primary"]),
        "determinacy": ("W_primary", metrics["W_primary"]),
        "hp-h3": ("W_tree", metrics["W_tree"]),
        "region-size": ("W_tree", metrics["W_tree"]),
    }
    metrics["expected_renderer_sha"] = EXPECTED_RENDERER_SHA256
    panel_failures = []
    for panel in metrics["panel_metrics"]:
        denominator_key, denominator = expected_denominators[panel["id"]]
        actual_bins = [{"key": item["key"], "count": item["count"]} for item in panel["bins"]]
        if not (
            panel["denominator_key"] == denominator_key
            and panel["total"] == denominator
            and panel["sum"] == denominator
            and actual_bins == panel["expected_bins"]
            and panel["visible_denominator"]
        ):
            panel_failures.append(panel["id"])
    metrics["panel_failures"] = panel_failures

    mode_states = []
    for mode in MODES:
        page.locator(f'[data-ideogram-mode="{mode}"]').click()
        mode_states.append(
            page.evaluate(
                """mode => ({
                    requested: mode,
                    active: document.getElementById('genome-ideogram').dataset.mode,
                    pressed: document.querySelectorAll('[data-ideogram-mode][aria-pressed="true"]').length,
                    marks: document.querySelectorAll('.ideogram-mark').length,
                    unmapped: [...document.querySelectorAll('.ideogram-mark')].filter(mark => !mark.dataset.modeValue).length,
                    legend: document.querySelectorAll('[data-legend-key]').length
                })""",
                mode,
            )
        )
        if sample == "HCC1395" and viewport in {"desktop1440", "mobile390"}:
            capture_section(
                page,
                "#genome-overview",
                SCREENSHOTS / f"{viewport}__{sample}__genome__{mode}.png",
            )
    metrics["mode_states"] = mode_states

    page.locator('[data-ideogram-mode="determinacy"]').click()
    first_legend = page.locator("[data-legend-key]").first
    first_legend.click()
    metrics["legend_isolation"] = page.evaluate(
        """() => ({
            pressed: document.querySelectorAll('[data-legend-key][aria-pressed="true"]').length,
            dimmed: document.querySelectorAll('.ideogram-mark.dimmed').length
        })"""
    )
    first_legend.click()

    hit = page.locator('.ideogram-chrom-hit[data-chrom="chr8"]')
    hit.focus()
    hit.press("Enter")
    page.wait_for_timeout(100)
    metrics["chromosome_sync"] = page.evaluate(
        """() => ({
            filter: document.getElementById('fchr').value,
            grid: document.querySelector('.chrom-button[data-chrom="chr8"]').getAttribute('aria-pressed'),
            ideogram: document.querySelector('.ideogram-chrom-hit[data-chrom="chr8"]').getAttribute('aria-pressed'),
            hash: new URLSearchParams(location.hash.slice(1)).get('chr')
        })"""
    )
    page.locator("#all-genome").click()

    metrics["pass"] = all(
        (
            not panel_failures,
            metrics["weighted_site_sum"] == metrics["retained_sSNV"],
            metrics["region_index_count"] == metrics["W_tree"],
            metrics["tracks"] == 22,
            metrics["marks"] == metrics["W_tree"],
            metrics["unique_marks"] == metrics["W_tree"],
            metrics["coordinate_errors"] == 0,
            metrics["identity_errors"] == 0,
            metrics["maximum_x_error"] <= 0.01,
            metrics["maximum_track_ratio_error"] <= 0.00001,
            metrics["mode_buttons"] == 5,
            metrics["mode_min_height"] >= 44,
            metrics["body_overflow_px"] <= 1,
            metrics["json_links"] >= 3,
            metrics["json_links_outside_closed_drawer"] == 0,
            metrics["scroller"]["role"] == "region",
            metrics["scroller"]["tabindex"] == "0",
            bool(metrics["scroller"]["aria_label"]),
            metrics["scroller"]["overflow_x"] == "auto",
            metrics["h1_count"] == 1,
            metrics["overview_labelled"] == "sample-overview-title",
            metrics["genome_labelled"] == "genome-title",
            metrics["svg_role"] == "img",
            metrics["svg_title"],
            metrics["svg_desc"],
            metrics["ui_contract"] == "layered-workstation-v5-grch38-overview-1",
            metrics["renderer_sha"] == metrics["expected_renderer_sha"],
            all(
                state["requested"] == state["active"]
                and state["pressed"] == 1
                and state["marks"] == metrics["W_tree"]
                and state["unmapped"] == 0
                and state["legend"] >= 1
                for state in mode_states
            ),
            metrics["legend_isolation"]["pressed"] == 1,
            metrics["legend_isolation"]["dimmed"] > 0,
            metrics["chromosome_sync"] == {
                "filter": "chr8",
                "grid": "true",
                "ideogram": "true",
                "hash": "chr8",
            },
        )
    )
    return metrics


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    SCREENSHOTS.mkdir(parents=True, exist_ok=True)
    runs = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        for viewport_name, viewport in VIEWPORTS.items():
            for identity in ("index", *SAMPLES):
                path = WORKSTATION / ("index.html" if identity == "index" else f"{identity}.html")
                page = browser.new_page(viewport=viewport, color_scheme="light", locale="zh-TW")
                console_errors = []
                page_errors = []
                request_failures = []
                external_requests = []
                page.on(
                    "console",
                    lambda message, bucket=console_errors: bucket.append(message.text)
                    if message.type == "error"
                    else None,
                )
                page.on("pageerror", lambda error, bucket=page_errors: bucket.append(str(error)))
                page.on(
                    "requestfailed",
                    lambda request, bucket=request_failures: bucket.append(
                        {"url": request.url, "failure": request.failure}
                    ),
                )
                page.on(
                    "request",
                    lambda request, bucket=external_requests: bucket.append(request.url)
                    if request.url.startswith(("http://", "https://"))
                    else None,
                )
                page.goto(path.as_uri(), wait_until="load", timeout=120_000)
                page.wait_for_load_state("networkidle", timeout=120_000)
                body_overflow = page.evaluate(
                    "Math.max(document.documentElement.scrollWidth, document.body.scrollWidth) - innerWidth"
                )
                run = {
                    "identity": identity,
                    "viewport": viewport_name,
                    "viewport_size": viewport,
                    "input": str(path),
                    "body_overflow_px": max(0, body_overflow),
                    "console_errors": console_errors,
                    "page_errors": page_errors,
                    "request_failures": request_failures,
                    "external_requests": external_requests,
                }
                if identity == "index":
                    run["sample_links"] = page.locator("a.genome-link").count()
                    run["pass"] = (
                        run["sample_links"] == 7
                        and run["body_overflow_px"] <= 1
                        and not console_errors
                        and not page_errors
                        and not request_failures
                        and not external_requests
                    )
                    capture_section(
                        page,
                        "#launch-title",
                        SCREENSHOTS / f"{viewport_name}__index__launch.png",
                    )
                else:
                    capture_section(
                        page,
                        "#sample-overview",
                        SCREENSHOTS / f"{viewport_name}__{identity}__overview.png",
                    )
                    capture_section(
                        page,
                        "#genome-overview",
                        SCREENSHOTS / f"{viewport_name}__{identity}__genome.png",
                    )
                    run["sample_metrics"] = audit_sample(page, identity, viewport_name)
                    run["pass"] = (
                        run["sample_metrics"]["pass"]
                        and not console_errors
                        and not page_errors
                        and not request_failures
                        and not external_requests
                    )
                runs.append(run)
                print(f"[{viewport_name}] {identity}: {'PASS' if run['pass'] else 'FAIL'}")
                page.close()
        browser.close()

    sample_runs = [run for run in runs if run["identity"] != "index"]
    aggregate = {
        "W_tree": sum(run["sample_metrics"]["W_tree"] for run in sample_runs if run["viewport"] == "desktop1440"),
        "W_primary": sum(
            run["sample_metrics"]["W_primary"] for run in sample_runs if run["viewport"] == "desktop1440"
        ),
        "retained_sSNV": sum(
            run["sample_metrics"]["retained_sSNV"]
            for run in sample_runs
            if run["viewport"] == "desktop1440"
        ),
    }
    report = {
        "schema_version": "1.0",
        "task_type": "B_comprehensive_validation",
        "partial": False,
        "browser": "Chromium via Python Playwright",
        "viewports": VIEWPORTS,
        "samples": list(SAMPLES),
        "runs": runs,
        "summary": {
            "status": "pass" if all(run["pass"] for run in runs) else "fail",
            "page_runs": len(runs),
            "sample_page_runs": len(sample_runs),
            "passed": sum(run["pass"] for run in runs),
            "failed": sum(not run["pass"] for run in runs),
            "aggregate_once_per_dataset": aggregate,
        },
    }
    (OUT / "metrics.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report["summary"], ensure_ascii=False, indent=2))
    raise SystemExit(0 if report["summary"]["status"] == "pass" else 1)


if __name__ == "__main__":
    main()
