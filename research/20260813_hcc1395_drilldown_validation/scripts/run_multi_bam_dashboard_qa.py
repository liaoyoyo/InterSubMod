#!/usr/bin/env python3
"""Browser QA for the portable multi-BAM dashboard."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from playwright.sync_api import Page, sync_playwright


DEFAULT_BROWSER = Path(
    "/bip7_disk/liaoyoyo2001/.cache/ms-playwright/"
    "chromium_headless_shell-1223/chrome-headless-shell-linux64/chrome-headless-shell"
)

AUDIT_BLOCK_BY_TITLE = {
    "BAM／paired-normal／reference 輸入綁定與 bounded payload 身分": "b_bam_input_table",
    "LongPhase-S alignment tag 分母與守恆明細": "b_tag_table",
    "7 組資料的 topology 明細": "b_topology_table",
    "HCC1395 甲基化軸明細": "b_axis_table",
    "HCC1395 lineage/LCA 受控比較 gate": "b_lca_table",
    "Legacy／目前 browser QA 收據": "b_visual_table",
}


def normalized_text(locator: Any) -> str:
    return " ".join(locator.inner_text().split())


def metric_texts(page: Page) -> list[str]:
    cards = page.locator(".report-metric-card")
    return ["".join(cards.nth(index).inner_text().split()) for index in range(cards.count())]


def expected_metric_texts(
    scope: dict[str, Any], bam_scope: dict[str, Any]
) -> list[str]:
    def percent(value: Any) -> str:
        if value is None:
            return "n/a"
        return f"{round(float(value) * 100, 1):g}%"

    macro_tree = (
        "生物樣本macrotreecoverage"
        f"{percent(scope['biological_macro_tree_coverage_median'])}"
    )
    if scope["biological_macro_tree_coverage_q1"] is not None:
        macro_tree += (
            f"P25{percent(scope['biological_macro_tree_coverage_q1'])}"
            f"P75{percent(scope['biological_macro_tree_coverage_q3'])}"
        )
    macro_tree += f"生物樣本n{scope['biological_macro_n']}"

    macro_unique = (
        "生物樣本macrounique/tree"
        f"{percent(scope['biological_macro_unique_among_tree_median'])}"
    )
    if scope["biological_macro_unique_among_tree_q1"] is not None:
        macro_unique += (
            f"P25{percent(scope['biological_macro_unique_among_tree_q1'])}"
            f"P75{percent(scope['biological_macro_unique_among_tree_q3'])}"
        )
    macro_unique += f"生物樣本n{scope['biological_macro_n']}"

    return [
        f"目前資料集{scope['dataset_count']}技術replicate{scope['technical_replicate_n']}",
        f"生物樣本{scope['biological_count']}",
        f"HashPASS{scope['hash_match_n']}"
        f"IdentityPASS{scope['sample_identity_match_n']}"
        f"ReceiptPASS{scope['receipt_technical_pass_n']}"
        f"目前資料集{scope['technical_gate_total_n']}",
        f"全family完整{scope['all_families_complete_n']}目前資料集{scope['all_families_complete_total_n']}",
        f"全單位objective-certified{scope['all_objective_certified_n']}目前資料集{scope['all_objective_certified_total_n']}",
        macro_tree,
        macro_unique,
        f"已稽核下游資料集{scope['downstream_available_count']}目前資料集{scope['downstream_total_count']}",
        f"TumorpayloadMATCH{bam_scope['tumor_bounded_payload_match_n']}"
        f"NormalpayloadMATCH{bam_scope['normal_bounded_payload_match_n']}"
        f"StrictidentityMATCH{bam_scope['strict_storage_identity_match_n']}"
        f"Mountmetadatadrift{bam_scope['mount_device_drift_n']}",
        f"Tumorquickcheck{bam_scope['tumor_quickcheck_pass_n']}"
        f"Normalquickcheck{bam_scope['normal_quickcheck_pass_n']}"
        f"Referencecompatible{bam_scope['reference_compatible_n']}"
        f"TagreceiptPASS{bam_scope['tag_receipt_pass_n']}",
    ]


def filter_button(page: Page) -> Any:
    buttons = page.locator("button")
    for index in range(buttons.count()):
        candidate = buttons.nth(index)
        text = normalized_text(candidate)
        if text.startswith("選擇資料集") or text.startswith("Sample filter"):
            return candidate
    raise AssertionError("Sample filter button was not found")


def choose_sample(page: Page, sample: str) -> None:
    filter_button(page).click(force=True, timeout=5000)
    candidates = page.locator("button").filter(has_text=sample)
    candidates.first.wait_for(state="visible", timeout=3000)
    for index in range(candidates.count()):
        candidate = candidates.nth(index)
        if normalized_text(candidate) == sample:
            candidate.click(force=True, timeout=5000)
            page.wait_for_timeout(450)
            return
    raise AssertionError(f"Sample option was not found: {sample}")


def choose_sample_keyboard(page: Page, sample: str) -> None:
    control = filter_button(page)
    control.focus()
    page.keyboard.press("Enter")
    candidates = page.locator("button").filter(has_text=sample)
    candidates.first.wait_for(state="visible", timeout=3000)
    for index in range(candidates.count()):
        candidate = candidates.nth(index)
        if normalized_text(candidate) == sample:
            candidate.focus()
            page.keyboard.press("Enter")
            page.wait_for_timeout(450)
            return
    raise AssertionError(f"Keyboard sample option was not found: {sample}")


def viewport_receipt(page: Page) -> dict[str, int]:
    return page.evaluate(
        """
        () => ({
          innerWidth,
          clientWidth: document.documentElement.clientWidth,
          scrollWidth: document.documentElement.scrollWidth,
        })
        """
    )


def assert_no_overflow(receipt: dict[str, int]) -> None:
    assert receipt["scrollWidth"] <= receipt["clientWidth"] + 1, receipt


def table_text(page: Page, title: str) -> str:
    panel = table_panel(page, title)
    return normalized_text(panel)


def table_panel(page: Page, title: str) -> Any:
    block_id = AUDIT_BLOCK_BY_TITLE.get(title)
    if block_id is not None:
        toggle = page.locator(
            f'.intersubmod-audit-toggle[data-audit-block-id="{block_id}"]'
        )
        assert toggle.count() == 1, block_id
        if toggle.get_attribute("aria-expanded") != "true":
            toggle.click()
            page.wait_for_timeout(120)
    panel = page.locator("section.table-panel").filter(has_text=title).first
    assert panel.count() == 1, title
    return panel


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--html", type=Path, required=True)
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--figures-dir", type=Path, required=True)
    parser.add_argument("--browser", type=Path, default=DEFAULT_BROWSER)
    args = parser.parse_args()

    html_path = args.html.resolve()
    artifact_path = args.artifact.resolve()
    receipt_path = args.receipt.resolve()
    figures_dir = args.figures_dir.resolve()
    figures_dir.mkdir(parents=True, exist_ok=True)
    receipt_path.parent.mkdir(parents=True, exist_ok=True)

    all_screenshot = figures_dir / "19_multi_bam_dashboard_all_desktop.png"
    hcc_screenshot = figures_dir / "20_multi_bam_dashboard_hcc1395_selected.png"
    mobile_screenshot = figures_dir / "21_multi_bam_dashboard_mobile.png"
    fold_screenshot = figures_dir / "22_multi_bam_dashboard_folded_details.png"
    topology_screenshot = figures_dir / "23_multi_bam_dashboard_topology_chart.png"
    dorado_screenshot = figures_dir / "24_multi_bam_dashboard_dorado_selected.png"
    bam_input_screenshot = figures_dir / "25_multi_bam_dashboard_bam_input_table.png"
    tag_chart_screenshot = figures_dir / "26_multi_bam_dashboard_read_tag_chart.png"
    mobile_audit_screenshot = figures_dir / "27_multi_bam_dashboard_mobile_audit_detail.png"
    opportunity_screenshot = figures_dir / "28_multi_bam_dashboard_opportunity_fixed_labels.png"
    denominator_screenshot = figures_dir / "29_multi_bam_dashboard_denominator_rails.png"

    receipt: dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "html": str(html_path),
        "html_sha256": sha256(html_path),
        "artifact": str(artifact_path),
        "artifact_sha256": sha256(artifact_path),
        "browser": str(args.browser.resolve()),
        "assertions": [],
        "console_errors": [],
        "page_errors": [],
        "external_requests": [],
        "screenshots": {},
        "responsive": {},
    }
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    scope_rows = artifact["snapshot"]["datasets"]["scope_summary"]
    scopes = {row["sample_filter"]: row for row in scope_rows}
    bam_scope_rows = artifact["snapshot"]["datasets"]["bam_scope_summary"]
    bam_scopes = {row["sample_filter"]: row for row in bam_scope_rows}
    availability_rows = artifact["snapshot"]["datasets"]["availability"]
    availability_by_scope_component = {
        (row["sample_filter"], row["component"]): row for row in availability_rows
    }
    bam_identity_rows = artifact["snapshot"]["datasets"]["bam_input_identity"]
    sample_order = [row["sample_filter"] for row in scope_rows if row["sample_filter"] != "All"]
    assert len(sample_order) == 7 and len(scopes) == 8 and len(bam_scopes) == 8

    def passed(label: str) -> None:
        receipt["assertions"].append({"label": label, "status": "PASS"})

    def observe(page: Page) -> None:
        page.on(
            "console",
            lambda message: receipt["console_errors"].append(message.text)
            if message.type == "error"
            else None,
        )
        page.on("pageerror", lambda error: receipt["page_errors"].append(str(error)))
        page.on(
            "request",
            lambda request: receipt["external_requests"].append(request.url)
            if request.url.startswith(("http://", "https://"))
            else None,
        )

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            headless=True, executable_path=str(args.browser.resolve())
        )

        desktop = browser.new_page(viewport={"width": 1440, "height": 1000})
        observe(desktop)
        desktop.goto(html_path.as_uri())
        desktop.wait_for_selector("#data-analytics-portable-reader-root", timeout=8000)
        desktop.wait_for_timeout(650)

        title = normalized_text(desktop.locator("#data-analytics-portable-reader-root h1").first)
        assert title == "InterSubMod 多樣本分析總覽", title
        assert desktop.locator("html").get_attribute("lang") == "zh-Hant"
        passed("desktop title")
        partial_pill = desktop.locator(".intersubmod-partial-pill")
        assert partial_pill.count() == 1
        assert normalized_text(partial_pill) == "PARTIAL"
        assert desktop.locator("#intersubmod-filter-live").get_attribute("aria-live") == "polite"
        passed("persistent PARTIAL authority and polite filter live region are present")
        top_bar_contrast = desktop.locator(".analytics-top-bar h1").evaluate(
            "node => ({color: getComputedStyle(node).color, "
            "background: getComputedStyle(node.closest('.analytics-top-bar')).backgroundColor})"
        )
        assert top_bar_contrast == {
            "color": "rgb(246, 251, 248)",
            "background": "rgba(10, 55, 52, 0.96)",
        }, top_bar_contrast
        passed("top-bar title uses high-contrast light text on the dark evidence header")
        body_text = normalized_text(desktop.locator("body"))
        assert "下游證據缺口仍在；維持 PARTIAL" in body_text
        assert "source query could not complete" not in body_text
        passed("known data gaps use accurate PARTIAL copy rather than a query-failure message")
        blocker_details = desktop.locator(".access-issue-strip details")
        assert blocker_details.count() == 2
        assert desktop.locator(".access-issue-strip details[open]").count() == 0
        blocker_details.nth(0).locator("summary").click()
        assert "v1 source 缺失" in normalized_text(blocker_details.nth(0))
        blocker_details.nth(0).locator("summary").click()
        assert desktop.locator(".access-issue-strip details[open]").count() == 0
        passed("two blocker names stay visible while detailed reasons are collapsed on demand")

        all_viewport = viewport_receipt(desktop)
        assert_no_overflow(all_viewport)
        passed("desktop All has no page-level horizontal overflow")
        all_metrics = metric_texts(desktop)
        expected_all = expected_metric_texts(scopes["All"], bam_scopes["All"])
        assert all_metrics == expected_all, {"actual": all_metrics, "expected": expected_all}
        assert all_metrics[:8] == [
            "目前資料集7技術replicate1",
            "生物樣本6",
            "HashPASS7IdentityPASS7ReceiptPASS7目前資料集7",
            "全family完整0目前資料集7",
            "全單位objective-certified0目前資料集7",
            "生物樣本macrotreecoverage77.7%P2574%P7578.6%生物樣本n6",
            "生物樣本macrounique/tree62.4%P2546.7%P7584.2%生物樣本n6",
            "已稽核下游資料集1目前資料集7",
        ]
        passed("All-scope hero and bounded BAM readiness metrics match reviewed snapshot")
        topology_box = desktop.locator("#ch_topology_rates").bounding_box()
        assert topology_box and topology_box["y"] < 1000 and topology_box["y"] + topology_box["height"] > 0
        passed("primary topology chart begins inside the 1440x1000 decision viewport")
        desktop.screenshot(path=str(all_screenshot), full_page=False)
        desktop.locator("#ch_topology_rates").screenshot(path=str(topology_screenshot))
        assert desktop.locator(".portable-chart-data-has-vector").count() == 6
        passed("all six charts retain semantic tabular fallbacks in the portable HTML")

        audit_toggles = desktop.locator(".intersubmod-audit-toggle")
        assert audit_toggles.count() == 6
        assert all(
            audit_toggles.nth(index).get_attribute("aria-expanded") == "false"
            for index in range(audit_toggles.count())
        )
        assert all(
            desktop.locator(
                '[data-analytics-layout-item][data-artifact-block-id="'
                + audit_toggles.nth(index).get_attribute("data-audit-block-id")
                + '"]'
            ).is_hidden()
            for index in range(audit_toggles.count())
        )
        assert all(
            audit_toggles.nth(index).get_attribute("aria-controls")
            == "intersubmod-audit-panel-"
            + audit_toggles.nth(index).get_attribute("data-audit-block-id")
            for index in range(audit_toggles.count())
        )
        passed("six secondary audit tables are collapsed by default with semantic buttons")
        passed("all audit disclosure buttons bind to labeled regions with aria-controls")

        coverage_rail = table_panel(desktop, "Active-k exact N/D")
        assert coverage_rail.locator("tbody tr").count() == 8
        coverage_text = normalized_text(coverage_rail)
        assert "1 3,202 2,735 85.4% 2,753 86%" in coverage_text
        assert "8+ 989 406 41.1% 408 41.3%" in coverage_text
        axis_rail = table_panel(desktop, "甲基化軸 exact N/D 與 validity rail")
        assert axis_rail.locator("tbody tr").count() == 5
        axis_rail_text = normalized_text(axis_rail)
        assert "cluster INVALID · DOUBLE-DIPPING" in axis_rail_text
        assert "11,545 / 11,546 · 100.0%" in axis_rail_text
        assert "11,779 / 11,779 · 100.0%" in axis_rail_text
        axis_chart_rows = [
            row
            for row in artifact["snapshot"]["datasets"]["axis_rates"]
            if row["sample_filter"] == "All"
        ]
        assert axis_chart_rows and all(row["axis"] != "cluster" for row in axis_chart_rows)
        passed("active-k chart has an adjacent eight-row exact N/D rail")
        passed("invalid cluster rows retain exact N/D in the validity rail and are excluded from the valid-axis chart")

        read_tag_chart_text = normalized_text(desktop.locator("#ch_read_tag_rates"))
        assert "Duplicate identities / all alignments" in read_tag_chart_text
        all_tag_rates = [
            row
            for row in artifact["snapshot"]["datasets"]["read_tag_rates"]
            if row["sample_filter"] == "All"
            and row["metric"] == "Duplicate identities / all alignments"
        ]
        assert len(all_tag_rates) == 7
        assert any(
            row["sample"] == "HCC1395"
            and abs(row["rate"] - 0.16139643) < 1e-8
            for row in all_tag_rates
        )
        assert any(
            row["sample"] == "HCC1937"
            and abs(row["rate"] - 0.25942924) < 1e-8
            for row in all_tag_rates
        )
        passed("duplicate-alignment denominator guard is visible and preserves HCC1395/HCC1937 rates")

        all_bam_panel = table_panel(desktop, "BAM／paired-normal／reference 輸入綁定與 bounded payload 身分")
        all_bam_rows = all_bam_panel.locator("tbody tr")
        assert all_bam_rows.count() == 7
        all_bam_text = normalized_text(all_bam_panel)
        for row in (
            item for item in bam_identity_rows if item["sample_filter"] == "All"
        ):
            row_text = next(
                normalized_text(all_bam_rows.nth(index))
                for index in range(all_bam_rows.count())
                if row["sample"] in normalized_text(all_bam_rows.nth(index))
            )
            for field in (
                "tumor_strict_identity",
                "normal_strict_identity",
                "reference_strict_identity",
            ):
                assert row[field] in row_text
        assert all_bam_text.count("UNAVAILABLE_NO_RG") >= 7
        assert all_bam_text.count("CROSS_DIRECTORY_REVIEW_REQUIRED") == 1
        assert all_bam_text.count("NOT_COLLECTED") >= 7
        all_tag_panel = table_panel(desktop, "LongPhase-S alignment tag 分母與守恆明細")
        assert all_tag_panel.locator("tbody tr").count() == 7
        all_tag_text = normalized_text(all_tag_panel)
        assert "all_captured_mapped_alignment_records_including_primary_secondary_supplementary" in all_tag_text
        passed("All-scope BAM table preserves mount drift, no-RG, one pairing review, and no full SHA")
        passed("All-scope tag table preserves seven exact alignment-record denominator rows")

        per_scope: dict[str, Any] = {"All": {"metrics": all_metrics, "viewport": all_viewport}}
        for sample in sample_order:
            choose_sample(desktop, sample)
            assert normalized_text(filter_button(desktop)).endswith(sample)
            sample_viewport = viewport_receipt(desktop)
            assert_no_overflow(sample_viewport)
            actual_metrics = metric_texts(desktop)
            expected_metrics = expected_metric_texts(scopes[sample], bam_scopes[sample])
            assert actual_metrics == expected_metrics, {
                "sample": sample,
                "actual": actual_metrics,
                "expected": expected_metrics,
            }
            topology_text = table_text(desktop, "7 組資料的 topology 明細")
            assert sample in topology_text
            topology_rows = table_panel(
                desktop, "7 組資料的 topology 明細"
            ).locator("tbody tr")
            topology_datasets = [
                normalized_text(topology_rows.nth(index).locator("td").nth(0))
                for index in range(topology_rows.count())
                if topology_rows.nth(index).locator("td").count()
            ]
            assert topology_datasets == [sample], topology_datasets
            availability_text = table_text(desktop, "跨樣本資料可用性")
            expected_identity_availability = availability_by_scope_component[
                (sample, "Tumor + normal bounded payload identity")
            ]["status"]
            assert expected_identity_availability in availability_text
            assert "AVAILABLE_EXISTING_PRODUCER_RECEIPT" in availability_text
            bam_panel = table_panel(desktop, "BAM／paired-normal／reference 輸入綁定與 bounded payload 身分")
            bam_rows = bam_panel.locator("tbody tr")
            assert bam_rows.count() == 1
            bam_text = normalized_text(bam_panel)
            assert sample in bam_text
            expected_bam_row = next(
                row
                for row in bam_identity_rows
                if row["sample_filter"] == sample and row["sample"] == sample
            )
            for field in (
                "tumor_strict_identity",
                "normal_strict_identity",
                "reference_strict_identity",
            ):
                assert expected_bam_row[field] in bam_text
            assert "UNAVAILABLE_NO_RG" in bam_text
            assert "NOT_COLLECTED" in bam_text
            tag_panel = table_panel(desktop, "LongPhase-S alignment tag 分母與守恆明細")
            tag_rows = tag_panel.locator("tbody tr")
            assert tag_rows.count() == 1
            tag_text = normalized_text(tag_panel)
            assert sample in tag_text and "PASS" in tag_text
            all_scope_block = desktop.locator(
                '[data-artifact-block-id="b_four_layer"]'
            )
            assert all_scope_block.count() == 2
            assert all(
                all_scope_block.nth(index).is_hidden()
                for index in range(all_scope_block.count())
            )
            opportunity_block = desktop.locator(
                '[data-artifact-block-id="b_opportunity_fixed_labels"]'
            )
            assert opportunity_block.count() == 2
            assert all(
                opportunity_block.nth(index).is_hidden()
                for index in range(opportunity_block.count())
            )
            bundle_text = table_text(desktop, "HCC1395 v1/v3 bundle 診斷")
            axis_text = table_text(desktop, "HCC1395 甲基化軸明細")
            lca_text = table_text(desktop, "HCC1395 lineage/LCA 受控比較 gate")
            if sample == "HCC1395":
                assert "AVAILABLE_BLOCKED_OBSERVATION" in availability_text
                assert "BLOCKED" in bundle_text
                assert "DOUBLE_DIPPING_INVALID_AS_EVIDENCE" in axis_text
                assert "FAIL" in lca_text
                hcc_viewport = sample_viewport
                hcc_metrics = actual_metrics
                assert "CROSS_DIRECTORY_REVIEW_REQUIRED" in bam_text
                desktop.screenshot(path=str(hcc_screenshot), full_page=False)
            else:
                assert "CROSS_DIRECTORY_REVIEW_REQUIRED" not in bam_text
                assert availability_text.count("ABSENT_NO_EQUIVALENT_BUNDLE") >= 2
                assert "No rows match the selected filters" in bundle_text
                assert "No rows match the selected filters" in axis_text
                assert "No rows match the selected filters" in lca_text
                assert actual_metrics[7] == "已稽核下游資料集0目前資料集1"
            if sample == "HCC1395_DORADO":
                assert "HCC1395" in topology_text and "yes" in topology_text
                assert "macrotreecoveragen/a生物樣本n0" in actual_metrics[5]
                assert "macrounique/treen/a生物樣本n0" in actual_metrics[6]
                desktop.screenshot(path=str(dorado_screenshot), full_page=False)
            per_scope[sample] = {
                "metrics": actual_metrics,
                "viewport": sample_viewport,
                "availability_absent_n": availability_text.count(
                    "ABSENT_NO_EQUIVALENT_BUNDLE"
                ),
            }
        passed("all seven dataset selector states match artifact scope rows")
        passed("all seven selector states isolate exactly one BAM and one tag-detail row")
        passed("fixed All-scope evidence anchor is hidden outside the All selector state")
        passed("fixed-label All-scope opportunity panel is hidden outside All")
        passed("non-HCC states retain null/ABSENT and never borrow HCC1395 downstream rows")
        passed("HCC1395 exposes BLOCKED bundle, invalid cluster axis, and failed LCA A/B gate")
        passed("HCC1395_DORADO remains a technical replicate and is excluded from biological macro")
        passed("all eight selector states have no page-level overflow")

        choose_sample_keyboard(desktop, "H1437")
        assert normalized_text(filter_button(desktop)).endswith("H1437")
        assert metric_texts(desktop) == expected_metric_texts(
            scopes["H1437"], bam_scopes["H1437"]
        )
        assert normalized_text(desktop.locator(":focus")) == normalized_text(
            filter_button(desktop)
        )
        assert "H1437" in normalized_text(
            desktop.locator("#intersubmod-filter-live")
        )
        passed("sample selector supports Enter-key keyboard interaction")
        passed("sample selection returns focus to the control and announces the new state")

        choose_sample(desktop, "All")
        all_scope_blocks = desktop.locator('[data-artifact-block-id="b_four_layer"]')
        assert all_scope_blocks.count() == 2
        assert all(
            all_scope_blocks.nth(index).get_attribute("hidden") is None
            for index in range(all_scope_blocks.count())
        )
        assert desktop.locator(
            '[data-artifact-block-id="b_four_layer"]:visible'
        ).count() == 1
        opportunity_blocks = desktop.locator(
            '[data-artifact-block-id="b_opportunity_fixed_labels"]'
        )
        assert opportunity_blocks.count() == 2
        assert desktop.locator(
            '[data-artifact-block-id="b_opportunity_fixed_labels"]:visible'
        ).count() == 1
        passed("fixed four-layer evidence anchor returns in the All selector state")
        passed("fixed-label opportunity panel returns in the All selector state")

        html_frames = desktop.locator("iframe.report-html-frame")
        assert html_frames.count() == 3
        four_layer = desktop.frame_locator("iframe.report-html-frame").nth(0)
        assert four_layer.locator("script").count() == 0
        four_layer_text = normalized_text(four_layer.locator("body"))
        for required in (
            "Aggregate · n=6",
            "Canonical · HCC1395",
            "Extreme observed · H2009",
            "Technical control · HCC1395_DORADO",
            "PARTIAL · 描述性",
        ):
            assert required in four_layer_text
        passed("four-layer evidence strip is semantic, no-script, and uses reviewed anchors")

        opportunity_frame = desktop.frame_locator("iframe.report-html-frame").nth(1)
        opportunity_text = normalized_text(opportunity_frame.locator("body"))
        for required in (
            "COLO829 13,919",
            "H1437 17,598",
            "H2009 36,042",
            "HCC1395 11,590",
            "HCC1395_DORADO · technical 6,865",
            "HCC1937 5,195",
            "HCC1954 7,746",
            "Distinct sSNV",
        ):
            assert required in opportunity_text
        assert opportunity_frame.locator("script").count() == 0
        opportunity_frame.locator("body").screenshot(path=str(opportunity_screenshot))
        passed("fixed-label opportunity panel shows all seven datasets, exact counts, and no script")

        iframe = html_frames.nth(2)
        iframe.scroll_into_view_if_needed()
        desktop.wait_for_timeout(500)
        frame = desktop.frame_locator("iframe.report-html-frame").nth(2)
        details = frame.locator("details")
        assert details.count() == 4
        assert frame.locator("script").count() == 0
        assert frame.locator("details[open]").count() == 0
        details.nth(0).locator("summary").click()
        assert frame.locator("details[open]").count() == 1
        expanded_text = normalized_text(details.nth(0))
        assert "7 datasets = 6 biological + 1 technical replicate" in expanded_text
        assert "三個技術 gate 都是 7/7 PASS" in expanded_text
        details.nth(3).locator("summary").click()
        identity_detail = normalized_text(details.nth(3))
        all_bam_scope = bam_scopes["All"]
        assert (
            f"{all_bam_scope['strict_storage_identity_match_n']}/"
            f"{all_bam_scope['dataset_count']}"
        ) in identity_detail
        if all_bam_scope["mount_device_drift_n"]:
            assert (
                f"{all_bam_scope['mount_device_drift_n']}/"
                f"{all_bam_scope['dataset_count']}"
            ) in identity_detail
            assert "MOUNT_DEVICE_DRIFT_ONLY" in identity_detail
        else:
            assert "strict storage identity 也是" in identity_detail
        assert "沒有 @RG" in identity_detail
        assert "full BAM SHA-256 尚未執行" in identity_detail
        passed("folded detail block has four closed sections and opens on demand")
        passed("folded BAM identity detail matches dynamic strict/drift counts and preserves no-RG/full-hash caveats")
        passed("folded detail block is semantic no-script HTML")
        iframe.screenshot(path=str(fold_screenshot))

        desktop.locator("#ch_read_tag_rates").screenshot(path=str(tag_chart_screenshot))
        axis_rail.screenshot(path=str(denominator_screenshot))
        table_panel(desktop, "BAM／paired-normal／reference 輸入綁定與 bounded payload 身分").screenshot(
            path=str(bam_input_screenshot)
        )

        responsive_specs = [
            (1024, 768, "tablet_1024"),
            (512, 768, "zoom_200pct_equivalent"),
            (390, 844, "mobile_390"),
            (320, 568, "mobile_320"),
        ]
        mobile_viewport: dict[str, int] | None = None
        for width, height, label in responsive_specs:
            responsive = browser.new_page(viewport={"width": width, "height": height})
            observe(responsive)
            responsive.goto(html_path.as_uri())
            responsive.wait_for_selector(
                "#data-analytics-portable-reader-root", timeout=8000
            )
            responsive.wait_for_timeout(700)
            responsive_viewport = viewport_receipt(responsive)
            assert_no_overflow(responsive_viewport)
            assert len(metric_texts(responsive)) == 10
            assert normalized_text(filter_button(responsive)).endswith("All")
            if label == "mobile_390":
                responsive.evaluate("window.scrollTo(0, 0)")
                responsive.screenshot(path=str(mobile_screenshot), full_page=False)

            topology_table = table_panel(
                responsive, "7 組資料的 topology 明細"
            ).locator(".table-wrap")
            assert topology_table.count() == 1
            table_scroll = topology_table.evaluate(
                "node => ({scrollWidth: node.scrollWidth, clientWidth: node.clientWidth, "
                "overflowX: getComputedStyle(node).overflowX})"
            )
            assert table_scroll["scrollWidth"] >= table_scroll["clientWidth"]
            assert table_scroll["overflowX"] in {"auto", "scroll"}
            receipt["responsive"][label] = {
                "css_viewport": {"width": width, "height": height},
                "document": responsive_viewport,
                "topology_table": table_scroll,
            }
            if label == "zoom_200pct_equivalent":
                receipt["responsive"][label]["interpretation"] = (
                    "1024 physical pixels represented as a 512 CSS-pixel layout viewport"
                )
            if label == "mobile_390":
                mobile_viewport = responsive_viewport
                responsive.screenshot(
                    path=str(mobile_audit_screenshot), full_page=False
                )
            responsive.close()
        assert mobile_viewport is not None
        passed("1024px tablet view has no page-level horizontal overflow")
        passed("200% zoom-equivalent 512 CSS-pixel view has no page-level overflow")
        passed("390px mobile view retains filter, ten priority/input cards, and no page-level overflow")
        passed("320px narrow mobile view retains filter, ten priority/input cards, and no page-level overflow")
        passed("wide tables remain locally scrollable at every responsive viewport")

        browser.close()

    assert not receipt["console_errors"], receipt["console_errors"]
    assert not receipt["page_errors"], receipt["page_errors"]
    assert not receipt["external_requests"], receipt["external_requests"]
    passed("no browser console or page errors")
    passed("standalone dashboard makes no external HTTP(S) requests")

    receipt["desktop"] = {
        "all": {"viewport": all_viewport, "metrics": all_metrics},
        "hcc1395": {"viewport": hcc_viewport, "metrics": hcc_metrics},
        "per_scope": per_scope,
    }
    receipt["mobile"] = {"viewport": mobile_viewport}
    for label, path in {
        "all_desktop": all_screenshot,
        "hcc1395_selected": hcc_screenshot,
        "mobile": mobile_screenshot,
        "folded_details": fold_screenshot,
        "topology_chart": topology_screenshot,
        "dorado_selected": dorado_screenshot,
        "bam_input_table": bam_input_screenshot,
        "read_tag_chart": tag_chart_screenshot,
        "mobile_audit_detail": mobile_audit_screenshot,
        "opportunity_fixed_labels": opportunity_screenshot,
        "denominator_rails": denominator_screenshot,
    }.items():
        receipt["screenshots"][label] = {
            "path": str(path),
            "sha256": sha256(path),
        }
    receipt["status"] = "PASS"
    receipt["assertion_pass_n"] = len(receipt["assertions"])
    receipt_path.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": receipt["status"],
                "assertion_pass_n": receipt["assertion_pass_n"],
                "console_error_n": len(receipt["console_errors"]),
                "page_error_n": len(receipt["page_errors"]),
                "receipt": str(receipt_path),
                "screenshots": receipt["screenshots"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
