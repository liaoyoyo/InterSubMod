#!/usr/bin/env python3
"""Chromium desktop/mobile visual and interaction audit for exact-PS pages."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import struct
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from playwright.sync_api import Browser, Page, sync_playwright


HERE = Path(__file__).resolve().parent
TOPIC_ROOT = HERE.parent
REPO_ROOT = HERE.parents[2]
WORKSTATION = REPO_ROOT / "docs" / "methodology" / "_assets" / "layered_workstation"
BUILDER_SCRIPT = TOPIC_ROOT / "scripts" / "build_exact_ps_layered_workstation.py"
GENE_DRUG_BUILDER_SCRIPT = (
    TOPIC_ROOT / "scripts" / "build_exact_ps_gene_drug_annotation.py"
)
GENE_DRUG_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/"
    "observation_workspaces/20260725_exact_ps_gene_drug_annotation/all7_v2"
)
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
UI_CONTRACT = "layered-workstation-exact-ps-v3"
DISPLAY_LABELS = (
    "HCC1395_HKU",
    "HCC1395_NYGC",
    "HCC1937",
    "HCC1954",
    "H1437",
    "H2009",
    "COLO829",
)
ANNOTATION_REGION_COUNTS = {
    "HCC1395": (480, 174),
    "HCC1395_DORADO": (240, 91),
    "HCC1937": (189, 66),
    "HCC1954": (273, 90),
    "H1437": (537, 167),
    "H2009": (1379, 503),
    "COLO829": (456, 161),
}
PROFILE_MODE_LEGEND_COUNTS = {
    "resolution": 3,
    "morphology": 5,
    "active_k": 13,
}
SIMILARITY_MODE_LABELS = {
    "composite": "Composite",
    "resolution": "Determinacy",
    "morphology": "Topology / coarse",
    "active_k": "Complexity / active k",
    "status": "State",
    "hp_symmetric": "HP balance",
}


class AuditFailure(RuntimeError):
    """Raised on a browser contract violation."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditFailure(message)


def launch_browser(playwright: Any) -> tuple[Browser, Path]:
    try:
        return (
            playwright.chromium.launch(headless=True),
            Path(playwright.chromium.executable_path).resolve(),
        )
    except Exception:
        for candidate in (
            Path("/usr/bin/chromium"),
            Path("/usr/bin/chromium-browser"),
            Path("/usr/bin/google-chrome"),
        ):
            if candidate.is_file():
                return (
                    playwright.chromium.launch(
                        headless=True,
                        executable_path=str(candidate),
                        args=["--no-sandbox"],
                    ),
                    candidate.resolve(),
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


def displayed_region_count(page: Page) -> int:
    text = page.locator("#filterStatus").inner_text()
    match = re.search(r"顯示 ([0-9,]+) / ([0-9,]+) regions", text)
    require(match is not None, f"cannot parse filter status: {text}")
    return int(match.group(1).replace(",", ""))


def scroll_reachability(page: Page, selector: str) -> list[dict[str, Any]]:
    return page.locator(selector).evaluate_all(
        """elements => elements.map(element => {
          const original = element.scrollLeft;
          const maxScrollLeft = Math.max(0, element.scrollWidth - element.clientWidth);
          element.scrollLeft = maxScrollLeft;
          const reachedEnd = Math.abs(element.scrollLeft - maxScrollLeft) <= 1;
          element.scrollLeft = original;
          return {
            className: element.className,
            clientWidth: element.clientWidth,
            scrollWidth: element.scrollWidth,
            maxScrollLeft,
            reachedEnd,
            overflowX: getComputedStyle(element).overflowX
          };
        })"""
    )


def exercise_legend_mode(
    page: Page,
    *,
    sample: str,
    viewport: str,
    mode_key: str,
    total_groups: int,
    keep_second_selected: bool = False,
) -> dict[str, Any]:
    """Exercise all/none, multi-select union, and second-click semantics."""

    page.locator(f'[data-mode="{mode_key}"]').click()
    legend_buttons = page.locator("#legend .legend-btn")
    require(
        legend_buttons.count() >= 2,
        f"{sample} {viewport}: {mode_key} legend has fewer than two classes",
    )
    all_keys = legend_buttons.evaluate_all(
        "elements => elements.map(element => element.dataset.legendKey)"
    )
    for key in all_keys:
        page.locator(f'[data-legend-key="{key}"]').click()
    require(
        displayed_region_count(page) == total_groups,
        f"{sample} {viewport}: {mode_key} all-selected is not equivalent "
        "to none-selected",
    )
    for key in all_keys:
        page.locator(f'[data-legend-key="{key}"]').click()
    require(
        displayed_region_count(page) == total_groups
        and page.locator(
            '#legend .legend-btn[aria-pressed="true"]'
        ).count()
        == 0,
        f"{sample} {viewport}: {mode_key} none-selected reset failed",
    )

    selected_keys = all_keys[:2]
    for key in selected_keys:
        page.locator(f'[data-legend-key="{key}"]').click()
    expected_union = page.evaluate(
        """({modeKey, keys}) => {
          const rows=JSON.parse(document.getElementById('pageData').textContent).rows;
          const value=row => {
            if(modeKey==="status")return String(row.readStatus);
            if(modeKey==="k")return String(Math.min(6,Number(row.activeK)));
            if(modeKey==="annotation")return String(row.annotationClass);
            return String(row[modeKey]);
          };
          return rows.filter(row=>keys.includes(value(row))).length;
        }""",
        {"modeKey": mode_key, "keys": selected_keys},
    )
    require(
        all(
            page.locator(f'[data-legend-key="{key}"]').get_attribute(
                "aria-pressed"
            )
            == "true"
            for key in selected_keys
        )
        and displayed_region_count(page) == expected_union,
        f"{sample} {viewport}: {mode_key} multiselect union failed",
    )
    page.locator(f'[data-legend-key="{selected_keys[0]}"]').click()
    expected_second_only = page.evaluate(
        """({modeKey, key}) => {
          const rows=JSON.parse(document.getElementById('pageData').textContent).rows;
          const value=row => {
            if(modeKey==="status")return String(row.readStatus);
            if(modeKey==="k")return String(Math.min(6,Number(row.activeK)));
            if(modeKey==="annotation")return String(row.annotationClass);
            return String(row[modeKey]);
          };
          return rows.filter(row=>value(row)===key).length;
        }""",
        {"modeKey": mode_key, "key": selected_keys[1]},
    )
    require(
        page.locator(
            f'[data-legend-key="{selected_keys[0]}"]'
        ).get_attribute("aria-pressed")
        == "false"
        and displayed_region_count(page) == expected_second_only,
        f"{sample} {viewport}: {mode_key} second-click deselect failed",
    )
    if not keep_second_selected:
        page.locator(f'[data-legend-key="{selected_keys[1]}"]').click()
        require(
            displayed_region_count(page) == total_groups,
            f"{sample} {viewport}: {mode_key} final reset failed",
        )
    return {
        "keys": all_keys,
        "selectedKeys": selected_keys,
        "unionCount": expected_union,
        "secondOnlyCount": expected_second_only,
    }


def screenshot_record(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    require(
        payload[:8] == b"\x89PNG\r\n\x1a\n" and len(payload) >= 24,
        f"invalid PNG screenshot: {path}",
    )
    width, height = struct.unpack(">II", payload[16:24])
    return {
        "path": str(path),
        "bytes": len(payload),
        "width": width,
        "height": height,
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def file_record(path: Path) -> dict[str, Any]:
    path = path.resolve()
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            size += len(block)
            digest.update(block)
    return {
        "path": str(path),
        "bytes": size,
        "sha256": digest.hexdigest(),
    }


def audited_input_records() -> dict[str, Any]:
    return {
        "index": file_record(WORKSTATION / "index.html"),
        "samples": {
            sample: file_record(WORKSTATION / f"{sample}.html")
            for sample in SAMPLES
        },
        "scripts": {
            "audit": file_record(Path(__file__)),
            "builder": file_record(BUILDER_SCRIPT),
            "gene_drug_builder": file_record(GENE_DRUG_BUILDER_SCRIPT),
        },
        "gene_drug_sidecar": {
            "annotation": file_record(GENE_DRUG_ROOT / "annotation.v1.json"),
            "receipt": file_record(GENE_DRUG_ROOT / "receipt.json"),
        },
    }


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
    require(
        "CGC ∩ drug" in page.locator("#samples").inner_text()
        and bool(
            page.locator('meta[name="intersubmod-gene-drug-sha256"]').get_attribute(
                "content"
            )
        ),
        f"index {mode}: cancer-gene/drug intersection summary is missing",
    )
    require(
        page.locator("#profileChart .profile-row").count() == 7,
        f"index {mode}: expected seven sample profile rows",
    )
    require(
        page.locator("#similarityHeatmap tbody td").count() == 49,
        f"index {mode}: expected 7x7 similarity heatmap",
    )
    matrix_columns = page.locator(
        "#similarityHeatmap thead th:not(:first-child)"
    ).all_inner_texts()
    matrix_rows = page.locator(
        "#similarityHeatmap tbody th"
    ).all_inner_texts()
    require(
        tuple(matrix_columns) == DISPLAY_LABELS
        and tuple(matrix_rows) == DISPLAY_LABELS,
        f"index {mode}: comparison order/labels mismatch "
        f"columns={matrix_columns} rows={matrix_rows}",
    )
    profile_labels = page.locator("#profileChart .profile-sample").evaluate_all(
        """elements => elements.map(element =>
          Array.from(element.childNodes)
            .filter(node => node.nodeType === Node.TEXT_NODE)
            .map(node => node.textContent)
            .join("")
            .trim()
        )"""
    )
    require(
        tuple(profile_labels) == DISPLAY_LABELS,
        f"index {mode}: profile order/labels mismatch {profile_labels}",
    )
    require(
        page.locator("#topPairs li").count() == 5,
        f"index {mode}: expected top-five pair ranking",
    )
    for profile_mode, expected_legend_count in (
        PROFILE_MODE_LEGEND_COUNTS.items()
    ):
        button = page.locator(f'[data-profile-mode="{profile_mode}"]')
        button.click()
        require(
            button.get_attribute("aria-pressed") == "true"
            and page.locator("#profileChart .profile-row").count() == 7
            and page.locator("#profileLegend > span").count()
            == expected_legend_count,
            f"index {mode}: {profile_mode} profile switch failed",
        )
    for similarity_mode, expected_label in SIMILARITY_MODE_LABELS.items():
        button = page.locator(
            f'[data-similarity-mode="{similarity_mode}"]'
        )
        button.click()
        require(
            button.get_attribute("aria-pressed") == "true"
            and expected_label.lower()
            in page.locator("#pairRankingTitle").inner_text().lower()
            and page.locator("#similarityHeatmap tbody td").count() == 49
            and page.locator("#topPairs li").count() == 5,
            f"index {mode}: {similarity_mode} similarity switch failed",
        )
    page.locator('[data-profile-mode="active_k"]').click()
    page.locator('[data-similarity-mode="active_k"]').click()
    require(
        "0.926381" in page.locator("#pair").inner_text()
        and "443/564" in page.locator("#pair").inner_text(),
        f"index {mode}: HCC technical evidence ladder missing",
    )
    post_info = geometry(page)
    require(
        post_info["html"]["scrollWidth"]
        <= post_info["html"]["clientWidth"] + 1,
        f"index {mode}: post-interaction horizontal overflow {post_info['html']}",
    )
    page.screenshot(path=str(output / f"index_{mode}.png"), full_page=False)
    require(not console_errors, f"index {mode}: console errors {console_errors}")
    require(not external_requests, f"index {mode}: external requests {external_requests}")
    return {
        **post_info,
        "initial_html": info["html"],
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
    require(
        page.locator('[data-testid="annotation-overview"]').count() == 1
        and page.locator('[data-testid="annotation-filters"] button').count() == 3,
        f"{sample} {mode}: annotation overview/filter controls missing",
    )

    page.locator("#genome").scroll_into_view_if_needed()
    page.wait_for_timeout(250)
    total_groups = page.evaluate(
        "() => JSON.parse(document.getElementById('pageData').textContent).rows.length"
    )
    require(
        displayed_region_count(page) == total_groups,
        f"{sample} {mode}: initial unfiltered count mismatch",
    )
    legend_multiselect = {}
    for mode_key in ("coarse", "hp", "status", "k", "annotation", "resolution"):
        legend_multiselect[mode_key] = exercise_legend_mode(
            page,
            sample=sample,
            viewport=mode,
            mode_key=mode_key,
            total_groups=total_groups,
            keep_second_selected=mode_key == "resolution",
        )
    selected_keys = legend_multiselect["resolution"]["selectedKeys"]
    legend_buttons = page.locator("#legend .legend-btn")

    page.locator("#clearSearch").click()
    first_region = page.locator("#regionList .region-button").first
    require(first_region.count() == 1, f"{sample} {mode}: no region list item")
    first_region.click()
    page.wait_for_selector("#regionDetail .detail-id", timeout=10_000)
    require(
        page.locator("#regionDetail .detail-id").inner_text().startswith("chr"),
        f"{sample} {mode}: region detail did not render",
    )
    visible_tree_node_labels = page.locator(
        "#regionDetail .tree g.hidden text, "
        "#regionDetail .tree g.observed text"
    ).evaluate_all(
        "elements => elements.map(element => element.textContent || '')"
    )
    require(
        all(
            not str(label or "").startswith("H_")
            for label in visible_tree_node_labels
        ),
        f"{sample} {mode}: raw long state leaked into visible tree nodes "
        f"{visible_tree_node_labels}",
    )
    tied_overlay: dict[str, int] | None = None
    if sample == "H2009":
        tied_overlay = {
            "global_best_union_edges": page.locator(
                "#candidateExplorer .candidate-edge.edge-global-best"
            ).count(),
            "selected_exemplar_edges": page.locator(
                "#candidateExplorer .candidate-edge.edge-selected"
            ).count(),
            "table_global_best_edges": page.locator(
                ".candidate-edge-table tbody tr.edge-global-best"
            ).count(),
            "table_selected_edges": page.locator(
                ".candidate-edge-table tbody tr.edge-selected"
            ).count(),
        }
        require(
            tied_overlay["selected_exemplar_edges"] == 3
            and tied_overlay["table_selected_edges"] == 3
            and tied_overlay["global_best_union_edges"]
            == tied_overlay["table_global_best_edges"]
            and tied_overlay["global_best_union_edges"]
            > tied_overlay["selected_exemplar_edges"],
            f"{sample} {mode}: tied global-best union was confused with "
            f"selected deterministic exemplar {tied_overlay}",
        )

    annotation_summary = page.locator(
        '[data-testid="annotation-overview"]'
    ).inner_text()
    annotation_counts = page.evaluate(
        """() => {
          const summary=JSON.parse(document.getElementById('pageData').textContent).summary;
          return {
            cgcRegions:Number(summary.cgcRegions),
            cgcDrugRegions:Number(summary.cgcDrugRegions),
            metaCgc:Number(document.querySelector('meta[name="intersubmod-cgc-locus-regions"]').content),
            metaDrug:Number(document.querySelector('meta[name="intersubmod-cgc-drug-locus-regions"]').content)
          };
        }"""
    )
    expected_cgc, expected_cgc_drug = ANNOTATION_REGION_COUNTS[sample]
    require(
        annotation_counts
        == {
            "cgcRegions": expected_cgc,
            "cgcDrugRegions": expected_cgc_drug,
            "metaCgc": expected_cgc,
            "metaDrug": expected_cgc_drug,
        }
        and f"{expected_cgc:,}" in annotation_summary
        and f"{expected_cgc_drug:,}" in annotation_summary,
        f"{sample} {mode}: annotation summary/meta mismatch "
        f"{annotation_counts}",
    )

    page.locator(f'[data-legend-key="{selected_keys[0]}"]').click()
    page.locator('[data-annotation-filter="cgc_drug"]').click()
    page.wait_for_timeout(150)
    expected_and_count = page.evaluate(
        """keys => {
          const rows=JSON.parse(document.getElementById('pageData').textContent).rows;
          return rows.filter(row=>row.hasCgcDrug && keys.includes(String(row.resolution))).length;
        }""",
        selected_keys,
    )
    require(
        "癌症基因 ∩ approved 抗腫瘤藥"
        in page.locator("#filterStatus").inner_text()
        and displayed_region_count(page) == expected_and_count
        and 0 < expected_and_count <= expected_cgc_drug
        and page.locator("#regionList .region-button").count() >= 1,
        f"{sample} {mode}: topology-union ∩ CGC/drug AND mismatch "
        f"expected={expected_and_count}",
    )
    page.locator("#regionList .region-button").first.click()
    require(
        page.locator(
            '#regionDetail [data-testid="gene-drug-panel"][open] '
            ".gene-drug-table tbody tr.has-drug"
        ).count()
        >= 1,
        f"{sample} {mode}: filtered region lacks same-gene CGC/drug display",
    )
    annotation_scroll_geometry = scroll_reachability(
        page,
        "#regionDetail .locus-strip, "
        "#regionDetail .gene-drug-panel .table-wrap",
    )
    require(
        annotation_scroll_geometry
        and all(
            row["reachedEnd"]
            and (
                row["maxScrollLeft"] <= 1
                or row["overflowX"] in {"auto", "scroll"}
            )
            for row in annotation_scroll_geometry
        ),
        f"{sample} {mode}: annotation/locus horizontal content is not fully "
        f"reachable {annotation_scroll_geometry}",
    )

    page.locator('[data-annotation-filter="all"]').click()
    page.locator('[data-mode="annotation"]').click()
    require(
        page.locator("#legend .legend-btn").count() == 3
        and "癌症基因"
        in page.locator("#legend").inner_text(),
        f"{sample} {mode}: annotation genome legend failed",
    )
    page.locator('[data-legend-key="NO_CGC"]').click()
    page.locator('[data-annotation-filter="cgc"]').click()
    require(
        displayed_region_count(page) == expected_cgc
        and page.locator(
            '#legend .legend-btn[aria-pressed="true"]'
        ).count()
        == 0,
        f"{sample} {mode}: incompatible legend selection was orphaned "
        "after annotation filter switch",
    )
    page.locator('[data-annotation-filter="all"]').click()
    page.locator('[data-mode="resolution"]').click()

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
        require(
            page.locator("#regionDetail .locus-card").count() == 2
            and page.locator("#regionDetail .state-matrix tbody tr").count() == 5,
            f"{sample} {mode}: target locus/read matrix mismatch",
        )
        require(
            page.locator("#candidateExplorer").count() == 1
            and page.locator(".candidate-edge-table tbody tr").count() == 4
            and page.locator(".candidate-edge-table tbody tr.edge-selected").count()
            == 2,
            f"{sample} {mode}: complete candidate union/overlay mismatch",
        )
        detail_text = page.locator("#regionDetail").inner_text().replace(",", "")
        require(
            "87888228" in detail_text and "87928739" in detail_text,
            f"{sample} {mode}: singleton boundary explanation missing",
        )
        require(
            "PTEN" not in detail_text
            and "NO_HIT_EVALUATED"
            in page.locator(
                '#regionDetail [data-testid="gene-drug-panel"]'
            ).inner_text()
            and page.locator(
                "#regionDetail .gene-drug-table"
            ).count()
            == 0,
            f"{sample} {mode}: historical PTEN envelope leaked into "
            "exact-locus gene/drug display",
        )
        candidate_input = page.locator("#candidateIndex")
        require(
            candidate_input.count() == 1
            and candidate_input.get_attribute("max") == "3",
            f"{sample} {mode}: target candidate carousel count mismatch",
        )
        page.locator('[data-candidate-step="1"]').click()
        require(
            page.locator("#candidateIndex").input_value() == "3",
            f"{sample} {mode}: candidate next control failed",
        )
        page.locator('[data-candidate-step="-1"]').click()
        require(
            page.locator("#candidateIndex").input_value() == "2",
            f"{sample} {mode}: candidate previous control failed",
        )

    post_info = geometry(page)
    detail_geometry = page.locator("#regionDetail").evaluate(
        """element => ({
          clientWidth: element.clientWidth,
          scrollWidth: element.scrollWidth,
          left: element.getBoundingClientRect().left,
          right: element.getBoundingClientRect().right
        })"""
    )
    graph_geometry = page.locator(
        "#regionDetail .candidate-network, #regionDetail .tree"
    ).evaluate_all(
        """elements => elements.map(element => {
          const original = element.scrollLeft;
          const maxScrollLeft = Math.max(0, element.scrollWidth - element.clientWidth);
          element.scrollLeft = maxScrollLeft;
          const reachedEnd = Math.abs(element.scrollLeft - maxScrollLeft) <= 1;
          element.scrollLeft = original;
          return {
            className: element.className,
            clientWidth: element.clientWidth,
            scrollWidth: element.scrollWidth,
            scrollLeft: original,
            maxScrollLeft,
            reachedEnd,
            overflowX: getComputedStyle(element).overflowX
          };
        })"""
    )
    require(
        post_info["html"]["scrollWidth"]
        <= post_info["html"]["clientWidth"] + 1,
        f"{sample} {mode}: post-detail horizontal overflow {post_info['html']}",
    )
    require(
        detail_geometry["right"] <= post_info["viewport"]["w"] + 1
        or detail_geometry["clientWidth"] == 0,
        f"{sample} {mode}: detail panel exceeds viewport {detail_geometry}",
    )
    if mode == "desktop":
        require(
            all(
                row["scrollWidth"] <= row["clientWidth"] + 1
                for row in graph_geometry
            ),
            f"{sample} {mode}: graph clipped at default scroll position {graph_geometry}",
        )
    else:
        require(
            all(
                row["reachedEnd"]
                and (
                    row["maxScrollLeft"] <= 1
                    or row["overflowX"] in {"auto", "scroll"}
                )
                for row in graph_geometry
            ),
            f"{sample} {mode}: horizontally overflowing graph is not fully "
            f"reachable {graph_geometry}",
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
        **post_info,
        "initial_html": info["html"],
        "detail_geometry": detail_geometry,
        "graph_geometry": graph_geometry,
        "legend_buttons": legend_buttons.count(),
        "filter_status": page.locator("#filterStatus").inner_text(),
        "detail_id": page.locator("#regionDetail .detail-id").inner_text(),
        "visible_tree_node_labels": visible_tree_node_labels,
        "tied_overlay": tied_overlay,
        "annotation_summary": annotation_summary,
        "annotation_counts": annotation_counts,
        "annotation_and_count": expected_and_count,
        "annotation_scroll_geometry": annotation_scroll_geometry,
        "legend_multiselect": legend_multiselect,
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
        "schema_version": "1.1.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "authority": AUTHORITY,
        "ui_contract": UI_CONTRACT,
        "samples": list(samples),
        "viewports": VIEWPORTS,
        "checks": {},
        "pages": {},
    }
    try:
        audited_inputs_pre = audited_input_records()
        with sync_playwright() as playwright:
            browser, browser_executable = launch_browser(playwright)
            receipt["browser"] = {
                "version": browser.version,
                "executable": file_record(browser_executable),
            }
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
        audited_inputs_post = audited_input_records()
        require(
            audited_inputs_pre == audited_inputs_post,
            "audited HTML/scripts changed while Chromium audit was running",
        )
        receipt["audited_inputs"] = audited_inputs_pre
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
            "audited_inputs_stable": True,
        }
        expected_screenshots = {
            f"index_{mode}.png"
            for mode in VIEWPORTS
        }
        expected_screenshots.update(
            f"{sample}_{mode}_{panel}.png"
            for sample in samples
            for mode in VIEWPORTS
            for panel in ("genome", "detail", "top")
        )
        receipt["artifacts"] = {
            name: screenshot_record(output / name)
            for name in sorted(expected_screenshots)
        }
        receipt["checks"]["screenshots_recorded"] = (
            len(receipt["artifacts"]) == len(expected_screenshots)
            and all(
                row["bytes"] > 0 and row["width"] > 0 and row["height"] > 0
                for row in receipt["artifacts"].values()
            )
        )
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
