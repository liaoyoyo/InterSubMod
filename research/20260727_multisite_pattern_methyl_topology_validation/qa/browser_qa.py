#!/usr/bin/env python3
"""Browser QA for the standalone multi-sSNV pattern-methyl report."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from playwright.sync_api import Page, sync_playwright


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_report(page: Page, report: Path, errors: list[str]) -> None:
    page.emulate_media(reduced_motion="reduce")
    page.on("console", lambda event: errors.append(f"console:{event.type}:{event.text}") if event.type == "error" else None)
    page.on("pageerror", lambda error: errors.append(f"pageerror:{error}"))
    page.goto(report.as_uri())
    page.wait_for_load_state("networkidle")
    page.wait_for_selector("#caseRows tr")
    require(page.locator("#filterStatus").inner_text() == "顯示 1,045 / 1,045 個 formal units", "unexpected initial case count")
    require(page.locator("[data-testid='technical-summary']").inner_text().strip() != "", "technical summary is blank")
    direct_answer = page.locator("[data-testid='direct-answer']").inner_text()
    require("formal full-four RR/RA/AR/AA 單元為 0" in direct_answer, "direct answer omits the full-four result")
    require("HCC1395_DORADO chr2 AR|RA" in direct_answer, "direct answer omits the strongest secondary case")


def assert_no_global_overflow(page: Page, label: str) -> None:
    dimensions = page.evaluate(
        """() => ({
          clientWidth: document.documentElement.clientWidth,
          scrollWidth: document.documentElement.scrollWidth
        })"""
    )
    require(
        dimensions["scrollWidth"] <= dimensions["clientWidth"] + 1,
        f"global horizontal overflow at {label}: {dimensions}",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--data-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--result-json", type=Path, required=True)
    args = parser.parse_args()

    report = args.report.resolve(strict=True)
    data_json = args.data_json.resolve(strict=True)
    external_data = json.loads(data_json.read_text(encoding="utf-8"))
    output_dir = args.output_dir.resolve()
    result_json = args.result_json.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    for path in (
        output_dir / "01_desktop_overview.png",
        output_dir / "02_desktop_robust_detail.png",
        output_dir / "03_secondary_exploratory_detail.png",
        output_dir / "04_mobile_robust_detail.png",
        output_dir / "05_print_overview.png",
        output_dir / "06_desktop_read_cluster_hover.png",
        result_json,
    ):
        require(not path.exists(), f"refusing to overwrite QA output: {path}")

    errors: list[str] = []
    checks: dict[str, object] = {}
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)

        desktop = browser.new_page(viewport={"width": 1440, "height": 1000}, device_scale_factor=1)
        load_report(desktop, report, errors)
        require(
            desktop.evaluate("() => DATA") == external_data,
            "embedded report data does not equal external data JSON",
        )
        read_cluster_partition = {
            "available": external_data["aggregate"]["read_cluster_records"],
            "non_evaluable": external_data["aggregate"][
                "read_cluster_non_evaluable_records"
            ],
            "no_detail": external_data["aggregate"][
                "read_cluster_evaluable_no_detail_records"
            ],
            "n_gt_160": external_data["aggregate"][
                "read_cluster_source_null_records"
            ],
        }
        require(
            read_cluster_partition
            == {
                "available": 557,
                "non_evaluable": 234,
                "no_detail": 253,
                "n_gt_160": 1,
            },
            f"unexpected read-cluster availability partition: {read_cluster_partition}",
        )
        require(
            "234 不可評估 · 253 未達 detail gate · 1 N>160"
            in desktop.locator("#metricGrid").inner_text(),
            "aggregate UI omits the read-cluster unavailability partition",
        )
        checks["read_cluster_availability_partition"] = read_cluster_partition
        assert_no_global_overflow(desktop, "desktop initial")
        require(
            desktop.locator("#filterDataset option[value='COLO829']").count() == 1,
            "COLO829 zero-formal dataset is absent from the scope filter",
        )
        desktop.locator("#filterDataset").select_option("COLO829")
        require(
            desktop.locator("#filterStatus").inner_text()
            == "顯示 0 / 1,045 個 formal units",
            "COLO829 does not expose its zero-formal result",
        )
        filtered_metrics = desktop.locator("#metricGrid").inner_text()
        require(
            "目前篩選 formal units\n0" in filtered_metrics,
            "COLO829 filtered formal-unit metric is not zero",
        )
        for global_label in (
            "全域 Formal full-four",
            "全域 k≥3 formal units",
            "全域 Artifact markers",
            "全域 Partial read projections",
        ):
            require(
                global_label in filtered_metrics,
                f"fixed aggregate metric lacks a global label: {global_label}",
            )
        checks["global_metric_labels_under_zero_filter"] = True
        desktop.locator("#resetFilters").click()
        desktop.screenshot(path=output_dir / "01_desktop_overview.png")
        checks["desktop_initial_rows"] = desktop.locator("#caseRows tr").count()
        checks["sentinel_chr22_visible"] = "AARR:18 · RRAA:3" in desktop.locator("#sentinelGrid").inner_text()
        checks["sentinel_chr5_visible"] = "AAAA:1 · AAAR:4 · AARA:19 · AARR:2" in desktop.locator("#sentinelGrid").inner_text()
        checks["sentinels_marked_nonformal"] = desktop.locator(
            "#sentinelGrid [data-assessment='NOT_IN_FORMAL_EVIDENCE']"
        ).count() == 4

        desktop.locator("#filterAssessment").select_option("ROBUST_ASSOCIATION")
        desktop.wait_for_timeout(100)
        require(desktop.locator("#filterStatus").inner_text() == "顯示 3 / 1,045 個 formal units", "robust filter did not yield three cases")
        require(desktop.locator("#caseRows tr").count() == 3, "robust table row count mismatch")
        desktop.locator("#caseRows .table-action").first.click()
        desktop.wait_for_timeout(100)
        desktop_bounds = desktop.evaluate(
            """() => ({
              detailTop: document.querySelector("#caseDetail").getBoundingClientRect().top,
              filterBottom: document.querySelector(".filter-shell").getBoundingClientRect().bottom
            })"""
        )
        require(
            desktop_bounds["detailTop"] >= desktop_bounds["filterBottom"] + 4,
            f"sticky filter obscures case detail heading: {desktop_bounds}",
        )
        desktop.screenshot(path=output_dir / "02_desktop_robust_detail.png")
        robust_text = desktop.locator("#caseDetail").inner_text()
        require("H1437 · chr22" in robust_text, "first robust case is not H1437 chr22")
        require("RRAR / RRRA" in robust_text, "H1437 robust state pair is absent")
        cluster_contract = desktop.evaluate(
            """() => {
              const row=DATA.cases.find(item =>
                item.dataset === "H1437" &&
                item.chrom === "chr22" &&
                item.assessment === "ROBUST_ASSOCIATION"
              );
              const cluster=row?.detail?.read_cluster;
              const canvas=document.querySelector("#readClusterCanvas");
              const context=canvas?.getContext("2d");
              const pixels=context?.getImageData(0,0,canvas.width,canvas.height).data;
              const colors=new Set();
              if(pixels){
                const pixelCount=pixels.length/4;
                const stride=Math.max(1,Math.floor(pixelCount/5000));
                for(let pixel=0;pixel<pixelCount;pixel+=stride){
                  const offset=pixel*4;
                  if(pixels[offset+3] > 0){
                    colors.add(`${pixels[offset]},${pixels[offset+1]},${pixels[offset+2]}`);
                  }
                }
              }
              return {
                available:Boolean(cluster?.available),
                n:cluster?.n_reads,
                matrixBytes:cluster ? atob(cluster.matrix_u8_b64).length : 0,
                expectedMatrixBytes:cluster ? cluster.n_reads ** 2 : 0,
                canvasWidth:canvas?.width ?? 0,
                canvasHeight:canvas?.height ?? 0,
                sampledColorCount:colors.size,
                orderUsesPatternLabels:cluster?.order_uses_pattern_labels,
                orderMethod:cluster?.order_method,
                orderTieBreaker:cluster?.order_tie_breaker,
                patternLevels:cluster?.pattern_levels ?? [],
                readGroupLevels:cluster?.read_group_levels ?? [],
                readGroupLevelsPseudonymized:
                  cluster?.read_group_levels_pseudonymized
              };
            }"""
        )
        require(cluster_contract["available"], f"robust read cluster is unavailable: {cluster_contract}")
        require(
            cluster_contract["matrixBytes"] == cluster_contract["expectedMatrixBytes"],
            f"read-cluster matrix payload length mismatch: {cluster_contract}",
        )
        require(
            cluster_contract["canvasWidth"] > 0
            and cluster_contract["canvasHeight"] > 0
            and cluster_contract["sampledColorCount"] >= 8,
            f"read-cluster canvas was not materially drawn: {cluster_contract}",
        )
        require(
            cluster_contract["orderUsesPatternLabels"] is False
            and "average" in cluster_contract["orderMethod"].lower(),
            f"read order is not documented as label-independent average linkage: {cluster_contract}",
        )
        require(
            cluster_contract["orderTieBreaker"]
            == "SOURCE_LEAF_ORDINAL_FOR_EXACT_DISTANCE_TIES",
            f"read-cluster tie breaker is not explicit: {cluster_contract}",
        )
        require(
            {"RRAR", "RRRA"}.issubset(set(cluster_contract["patternLevels"])),
            f"robust pattern annotations are incomplete: {cluster_contract}",
        )
        require(
            cluster_contract["readGroupLevelsPseudonymized"] is True
            and all(
                level == "." or level.isdigit()
                for level in cluster_contract["readGroupLevels"]
            ),
            f"read-group legend exposes a non-pseudonymized level: {cluster_contract}",
        )
        cluster_canvas = desktop.locator("#readClusterCanvas")
        cluster_canvas.hover(position={"x": 300, "y": 300})
        tooltip_text = desktop.locator("#readClusterTooltip").inner_text()
        require(
            "Bernoulli distance" in tooltip_text and "r0" in tooltip_text,
            f"read-cluster hover tooltip is missing anonymous pair detail: {tooltip_text}",
        )
        cluster_canvas.focus()
        cluster_canvas.press("ArrowRight")
        cluster_live_text = desktop.locator("#readClusterLive").inner_text()
        require(
            "Bernoulli distance" in cluster_live_text,
            f"keyboard read-cluster navigation is not announced: {cluster_live_text}",
        )
        desktop.locator("[data-testid='read-cluster-heatmap']").screenshot(
            path=output_dir / "06_desktop_read_cluster_hover.png"
        )
        checks["read_cluster_contract"] = cluster_contract
        checks["read_cluster_hover_anonymous"] = True
        checks["read_cluster_keyboard_accessible"] = True
        marker_alignment = desktop.evaluate(
            """() => {
              const row=DATA.cases.find(item =>
                item.dataset === "H1437" && item.assessment === "ROBUST_ASSOCIATION"
              );
              const positions=row.detail.cpg_positions_displayed,cell=18,left=112;
              const centers=positions.map((_,index)=>left+index*cell+(cell-1)/2);
              const offsets=row.active_positions.map(position=>{
                const x=ordinalMarkerX(position,positions,left,cell);
                return Math.min(...centers.map(center=>Math.abs(center-x)));
              });
              return {maxNearestOffset:Math.max(...offsets),offsets};
            }"""
        )
        require(
            marker_alignment["maxNearestOffset"] <= 18,
            f"active marker is not aligned to the CpG ordinal axis: {marker_alignment}",
        )
        checks["marker_ordinal_max_nearest_px"] = marker_alignment[
            "maxNearestOffset"
        ]
        checks["robust_filtered_rows"] = 3

        secondary_id = desktop.evaluate(
            """() => DATA.cases.find(row =>
              row.dataset === "HCC1395_DORADO" &&
              row.chrom === "chr2" &&
              row.phase_set === "4101751" &&
              row.hp_raw === "2-1"
            )?.id"""
        )
        require(bool(secondary_id), "strongest secondary case was not found")
        desktop.locator("#resetFilters").click()
        require(
            desktop.locator("#filterStatus").inner_text()
            == "顯示 1,045 / 1,045 個 formal units",
            "filter reset failed before secondary screenshot",
        )
        desktop.evaluate("(id) => selectCase(id, true)", secondary_id)
        desktop.wait_for_timeout(100)
        desktop.screenshot(path=output_dir / "03_secondary_exploratory_detail.png")
        secondary_text = desktop.locator("#caseDetail").inner_text()
        require("HCC1395_DORADO · chr2" in secondary_text, "secondary case identity mismatch")
        require("0.791" in secondary_text, "secondary R-squared is absent")
        require("Evaluable · no robust association" in secondary_text, "secondary case is overstated")
        secondary_cluster = desktop.evaluate(
            """id => {
              const row=DATA.cases.find(item => item.id === id);
              return {
                available:Boolean(row?.detail?.read_cluster?.available),
                n:row?.detail?.read_cluster?.n_reads,
                patterns:row?.detail?.read_cluster?.pattern_levels ?? [],
                canvasPresent:Boolean(document.querySelector("#readClusterCanvas"))
              };
            }""",
            secondary_id,
        )
        require(
            secondary_cluster["available"]
            and secondary_cluster["canvasPresent"]
            and {"AR", "RA"}.issubset(set(secondary_cluster["patterns"])),
            f"secondary AR/RA read-cluster view is incomplete: {secondary_cluster}",
        )
        require(
            "HCC1395_DORADO" in desktop.locator(
                "[data-testid='strongest-secondary']"
            ).inner_text(),
            "strongest secondary fixed card is absent",
        )
        checks["secondary_case_id"] = secondary_id
        checks["secondary_read_cluster"] = secondary_cluster

        mobile = browser.new_page(viewport={"width": 390, "height": 844}, device_scale_factor=1)
        load_report(mobile, report, errors)
        assert_no_global_overflow(mobile, "mobile initial")
        mobile_robust_id = mobile.evaluate(
            """() => DATA.cases.find(row =>
              row.dataset === "H2009" && row.assessment === "ROBUST_ASSOCIATION"
            )?.id"""
        )
        require(bool(mobile_robust_id), "H2009 robust case was not found")
        mobile.evaluate("(id) => selectCase(id, false)", mobile_robust_id)
        mobile.locator("#caseDetail").scroll_into_view_if_needed()
        mobile.wait_for_timeout(100)
        assert_no_global_overflow(mobile, "mobile robust detail")
        mobile_cluster = mobile.evaluate(
            """() => {
              const scroll=document.querySelector(".read-cluster-scroll");
              const canvas=document.querySelector("#readClusterCanvas");
              return {
                canvasPresent:Boolean(canvas),
                canvasWidth:canvas?.getBoundingClientRect().width ?? 0,
                scrollClientWidth:scroll?.clientWidth ?? 0,
                scrollWidth:scroll?.scrollWidth ?? 0,
                overflowX:scroll ? getComputedStyle(scroll).overflowX : ""
              };
            }"""
        )
        require(
            mobile_cluster["canvasPresent"]
            and mobile_cluster["canvasWidth"] > mobile_cluster["scrollClientWidth"]
            and mobile_cluster["scrollWidth"] > mobile_cluster["scrollClientWidth"]
            and mobile_cluster["overflowX"] in {"auto", "scroll"},
            f"mobile read-cluster view does not use local horizontal scrolling: {mobile_cluster}",
        )
        mobile.screenshot(path=output_dir / "04_mobile_robust_detail.png")
        checks["mobile_width"] = mobile.evaluate("document.documentElement.clientWidth")
        checks["mobile_read_cluster_scroll"] = mobile_cluster

        print_page = browser.new_page(viewport={"width": 1440, "height": 1000}, device_scale_factor=1)
        load_report(print_page, report, errors)
        print_page.emulate_media(media="print")
        require(print_page.locator(".filter-shell").evaluate("node => getComputedStyle(node).display") == "none", "print CSS does not hide filters")
        print_robust_id = print_page.evaluate(
            """() => DATA.cases.find(row =>
              row.dataset === "H1437" && row.assessment === "ROBUST_ASSOCIATION"
            )?.id"""
        )
        print_page.evaluate("(id) => selectCase(id, false)", print_robust_id)
        require(
            print_page.locator("#readClusterCanvas").evaluate(
                "node => getComputedStyle(node).display"
            )
            != "none",
            "print CSS hides the read-cluster canvas",
        )
        print_page.screenshot(path=output_dir / "05_print_overview.png")
        checks["print_filter_hidden"] = True
        checks["print_read_cluster_visible"] = True

        browser.close()

    require(not errors, f"browser console/page errors: {errors}")
    require(checks["sentinel_chr22_visible"], "chr22 sentinel state counts are absent")
    require(checks["sentinel_chr5_visible"], "chr5 sentinel state counts are absent")
    require(checks["sentinels_marked_nonformal"], "sentinel non-formal labels are incomplete")
    payload = {
        "status": "PASS",
        "report": str(report),
        "data_json": str(data_json),
        "report_sha256": sha256_file(report),
        "data_json_sha256": sha256_file(data_json),
        "embedded_external_data_equal": True,
        "checks": checks,
        "browser_errors": errors,
        "screenshots": [str(path) for path in sorted(output_dir.glob("*.png"))],
    }
    result_json.parent.mkdir(parents=True, exist_ok=True)
    result_json.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
