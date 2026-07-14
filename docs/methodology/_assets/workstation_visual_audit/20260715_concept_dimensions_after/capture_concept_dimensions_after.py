#!/usr/bin/env python3
"""Chromium audit for the canonical three-dimension workstation redesign."""

from __future__ import annotations

import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from playwright.sync_api import Page, sync_playwright


ROOT = Path(__file__).resolve().parents[5]
OUT = Path(__file__).resolve().parent
SHOTS = OUT / "screenshots"
SOURCE = ROOT / "docs/methodology/_assets/layered_workstation"
PAGES = [
    "index.html",
    "HCC1395.html",
    "COLO829.html",
    "H1437.html",
    "H2009.html",
    "HCC1395_DORADO.html",
    "HCC1937.html",
    "HCC1954.html",
]
VIEWPORTS = {
    "desktop": {"width": 1440, "height": 1000},
    "tablet": {"width": 768, "height": 1024},
    "mobile": {"width": 390, "height": 844},
    "narrow": {"width": 320, "height": 720},
}
CLASS_EXPECTATIONS = {
    "exact_and_topology_unique": {"value": "Exact 唯一", "c": "1", "topo": "1"},
    "topology_unique_exact_multiple": {"value": "Shape 唯一", "c": ">1", "topo": "1"},
    "topology_multiple_exact_multiple": {"value": "Multi-shape", "c": ">1", "topo": ">1"},
    "incomplete": {"value": "尚未評估", "c": "—", "topo": "—"},
    "no_primary_lineage": {"value": "不適用", "c": "—", "topo": "—"},
}


def add_check(checks: list[dict[str, Any]], name: str, passed: bool, actual: Any, expected: Any) -> None:
    checks.append(
        {
            "name": name,
            "status": "pass" if passed else "fail",
            "expected": expected,
            "actual": actual,
        }
    )


def number_from_count(text: str) -> int:
    match = re.search(r"[\d,]+", text)
    return int(match.group(0).replace(",", "")) if match else -1


def readout_matches(actual: str, expected: str) -> bool:
    if expected == ">1":
        return number_from_count(actual) > 1
    return actual == expected


def runtime_handlers(page: Page, record: dict[str, Any]) -> None:
    page.on(
        "console",
        lambda message: record["console_errors"].append(message.text) if message.type == "error" else None,
    )
    page.on("pageerror", lambda error: record["page_errors"].append(str(error)))
    page.on("requestfailed", lambda request: record["request_failures"].append(request.url))


def screenshot_first(page: Page, page_name: str, viewport_name: str) -> dict[str, Any]:
    selector = 'section[aria-labelledby="dimensions-title"]' if page_name == "index.html" else "#dimension-guide"
    locator = page.locator(selector)
    locator.wait_for(state="visible", timeout=30_000)
    locator.scroll_into_view_if_needed()
    page.wait_for_timeout(120)
    target = SHOTS / f"{viewport_name}_{page_name.removesuffix('.html')}_dimensions.png"
    locator.screenshot(path=str(target), animations="disabled", caret="hide", timeout=120_000)
    return {
        "kind": "dimension_guide",
        "path": str(target.relative_to(ROOT)),
        "bytes": target.stat().st_size,
        "box": locator.bounding_box(),
    }


def audit_index(page: Page, viewport_name: str, checks: list[dict[str, Any]]) -> None:
    state = page.evaluate(
        """() => ({
            cards: document.querySelectorAll('.dimension-grid > .dimension-card').length,
            labels: Array.from(document.querySelectorAll('.dimension-code')).map(el => el.innerText.trim()),
            text: document.querySelector('section[aria-labelledby="dimensions-title"]').innerText,
            beforeTopology: !!(document.querySelector('section[aria-labelledby="dimensions-title"]')
                .compareDocumentPosition(document.querySelector('section[aria-labelledby="topology-title"]'))
                & Node.DOCUMENT_POSITION_FOLLOWING),
            bodyOverflow: Math.max(document.body.scrollWidth, document.documentElement.scrollWidth) - innerWidth,
            cardColumns: getComputedStyle(document.querySelector('.dimension-grid')).gridTemplateColumns.split(' ').length,
            positionLeaders: Array.from(document.querySelectorAll('.position-leaders a')).map(el=>({
                text:el.innerText.trim(),href:el.getAttribute('href')
            }))
        })"""
    )
    add_check(checks, f"index/{viewport_name}/three_cards", state["cards"] == 3, state["cards"], 3)
    add_check(
        checks,
        f"index/{viewport_name}/three_labels",
        all(term in " ".join(state["labels"]) for term in ("拓樸型態", "可辨識度", "基因體位置")),
        state["labels"],
        ["拓樸型態", "可辨識度", "基因體位置"],
    )
    add_check(checks, f"index/{viewport_name}/before_topology", state["beforeTopology"], state["beforeTopology"], True)
    add_check(checks, f"index/{viewport_name}/no_overflow", state["bodyOverflow"] <= 1, state["bodyOverflow"], "<=1px")
    expected_columns = 3 if viewport_name == "desktop" else 1
    add_check(
        checks,
        f"index/{viewport_name}/responsive_columns",
        state["cardColumns"] == expected_columns,
        state["cardColumns"],
        expected_columns,
    )
    add_check(
        checks,
        f"index/{viewport_name}/seven_position_leaders",
        len(state["positionLeaders"]) == 7
        and all("#chr=chr" in item["href"] for item in state["positionLeaders"]),
        state["positionLeaders"],
        "7 hash-linked W_primary chromosome leaders",
    )
    for phrase in ("22,319", "11,582", "10,737", "19,921", "7,975", "7 × chr1–22"):
        add_check(checks, f"index/{viewport_name}/text/{phrase}", phrase in state["text"], state["text"], phrase)


def sample_data_integrity(page: Page) -> dict[str, Any]:
    return page.evaluate(
        """() => {
            const data=JSON.parse(document.getElementById('workstation-data').textContent);
            const rows=data.region_index;
            const counts={};rows.forEach(row=>counts[row.topology_class]=(counts[row.topology_class]||0)+1);
            const invalid=rows.filter(row=>{
                if(row.region!==row.chrom+':'+row.start+'-'+row.end||row.start>row.end)return true;
                if(row.endpoint_distance_bp!==row.end-row.start)return true;
                if(row.inclusive_span_bp!==row.end-row.start+1)return true;
                if(row.inclusive_span_bp!==row.endpoint_distance_bp+1)return true;
                if(row.topology_class==='no_primary_lineage')return !(row.n_primary===0&&row.candidate_complete===null&&row.C===null&&row.Topo===null);
                if(row.topology_class==='incomplete')return !(row.n_primary>0&&row.candidate_complete===false&&row.C===null&&row.Topo===null);
                if(row.topology_class==='exact_and_topology_unique')return !(row.candidate_complete===true&&row.C===1&&row.Topo===1);
                if(row.topology_class==='topology_unique_exact_multiple')return !(row.candidate_complete===true&&row.C>1&&row.Topo===1);
                if(row.topology_class==='topology_multiple_exact_multiple')return !(row.candidate_complete===true&&row.C>1&&row.Topo>1);
                return true;
            }).map(row=>row.region);
            const chroms={};for(let n=1;n<=22;n++)chroms['chr'+n]={chrom:'chr'+n,W_tree:0,W_primary:0,incomplete:0};
            rows.forEach(row=>{const c=chroms[row.chrom];c.W_tree++;if(row.topology_class!=='no_primary_lineage'){c.W_primary++;if(row.topology_class==='incomplete')c.incomplete++;}});
            const rank=chrom=>Number(chrom.replace('chr',''));
            const countLeader=Object.values(chroms).sort((a,b)=>b.W_primary-a.W_primary||rank(a.chrom)-rank(b.chrom))[0];
            const incompleteLeader=Object.values(chroms).filter(row=>row.W_primary).sort((a,b)=>
                b.incomplete/b.W_primary-a.incomplete/a.W_primary||rank(a.chrom)-rank(b.chrom))[0];
            return {sample:data.sample,rows:rows.length,invalid,counts,countLeader,incompleteLeader,canonical:data.canonical_sample};
        }"""
    )


def audit_sample(page: Page, page_name: str, viewport_name: str, checks: list[dict[str, Any]]) -> None:
    integrity = sample_data_integrity(page)
    state = page.evaluate(
        """() => ({
            cards: document.querySelectorAll('#dimension-cards > .dimension-card').length,
            labels: Array.from(document.querySelectorAll('.dimension-label')).map(el=>el.innerText.trim()),
            guideBeforeGenome: !!(document.getElementById('dimension-guide').compareDocumentPosition(
                document.getElementById('genome-overview')) & Node.DOCUMENT_POSITION_FOLLOWING),
            glossaryItems: document.querySelectorAll('.glossary-item').length,
            glossaryOpen: document.querySelector('.dimension-glossary').open,
            glossaryColumns: getComputedStyle(document.querySelector('.glossary-grid')).gridTemplateColumns.split(' ').length,
            actionSizes: Array.from(document.querySelectorAll('.dimension-action')).map(el=>{
                const r=el.getBoundingClientRect();return {width:r.width,height:r.height};
            }),
            firstRowDimensions: document.querySelector('.result-row .result-dimensions')?.children.length||0,
            bodyOverflow: Math.max(document.body.scrollWidth,document.documentElement.scrollWidth)-innerWidth,
            visibleJsonLinks: Array.from(document.querySelectorAll('a[href*=".json"]')).filter(el=>{
                if(el.closest('details:not([open])'))return false;const r=el.getBoundingClientRect();return r.width>0&&r.height>0;
            }).length,
            text: document.getElementById('dimension-guide').innerText
        })"""
    )
    prefix = f"{page_name}/{viewport_name}"
    add_check(checks, f"{prefix}/integrity", not integrity["invalid"], integrity["invalid"], [])
    add_check(checks, f"{prefix}/three_cards", state["cards"] == 3, state["cards"], 3)
    add_check(
        checks,
        f"{prefix}/labels",
        all(term in " ".join(state["labels"]) for term in ("拓樸型態", "可辨識度", "基因體位置")),
        state["labels"],
        ["拓樸型態", "可辨識度", "基因體位置"],
    )
    add_check(checks, f"{prefix}/guide_before_genome", state["guideBeforeGenome"], state["guideBeforeGenome"], True)
    add_check(checks, f"{prefix}/glossary", state["glossaryItems"] == 4 and not state["glossaryOpen"], state, "4 items; closed")
    expected_glossary_columns = 4 if viewport_name == "desktop" else 2 if viewport_name == "tablet" else 1
    add_check(
        checks,
        f"{prefix}/glossary_columns",
        state["glossaryColumns"] == expected_glossary_columns,
        state["glossaryColumns"],
        expected_glossary_columns,
    )
    add_check(
        checks,
        f"{prefix}/hit_targets",
        all(size["width"] >= 44 and size["height"] >= 44 for size in state["actionSizes"]),
        state["actionSizes"],
        ">=44x44",
    )
    add_check(checks, f"{prefix}/row_three_dimensions", state["firstRowDimensions"] == 3, state["firstRowDimensions"], 3)
    add_check(checks, f"{prefix}/no_overflow", state["bodyOverflow"] <= 1, state["bodyOverflow"], "<=1px")
    add_check(checks, f"{prefix}/json_hidden", state["visibleJsonLinks"] == 0, state["visibleJsonLinks"], 0)
    counts = integrity["counts"]
    expected_text = [
        f"{counts['exact_and_topology_unique']:,}",
        f"{counts['topology_unique_exact_multiple']:,}",
        f"{counts['topology_multiple_exact_multiple']:,}",
        f"{counts['incomplete']:,}",
        integrity["countLeader"]["chrom"],
        integrity["incompleteLeader"]["chrom"],
    ]
    add_check(
        checks,
        f"{prefix}/data_driven_text",
        all(text in state["text"] for text in expected_text),
        state["text"],
        expected_text,
    )

    for topology_class, expectation in CLASS_EXPECTATIONS.items():
        row = page.evaluate(
            """topologyClass => {
                const data=JSON.parse(document.getElementById('workstation-data').textContent);
                return data.region_index.find(row=>row.topology_class===topologyClass);
            }""",
            topology_class,
        )
        page.evaluate("region => selectRegion(region,null,false,'none')", row["region"])
        page.wait_for_function(
            "region => document.querySelector('#detail-title')?.textContent.trim()===region",
            arg=row["region"],
            timeout=30_000,
        )
        detail = page.evaluate(
            """() => ({
                cards: document.querySelectorAll('#detail .region-dimension').length,
                text: document.querySelector('#detail .region-dimensions').innerText,
                c: document.querySelectorAll('.ct-box strong')[0].innerText.trim(),
                topo: document.querySelectorAll('.ct-box strong')[1].innerText.trim(),
                beforeAssertions: !!(document.querySelector('.region-dimensions').compareDocumentPosition(
                    document.querySelector('.assertion-grid')) & Node.DOCUMENT_POSITION_FOLLOWING)
            })"""
        )
        passed = (
            detail["cards"] == 3
            and expectation["value"] in detail["text"]
            and readout_matches(detail["c"], expectation["c"])
            and readout_matches(detail["topo"], expectation["topo"])
            and detail["beforeAssertions"]
            and row["region"] in detail["text"]
        )
        add_check(checks, f"{prefix}/detail/{topology_class}", passed, detail, expectation)

    if viewport_name == "desktop":
        page.locator("#fevidence").select_option("partial_only")
        page.locator("#fsignal").select_option("recurrence")
        page.locator("#fq").fill("chr22:999999999")
        page.locator('[data-dimension-topology="topology_multiple_exact_multiple"]').click()
        page.wait_for_timeout(100)
        action = page.evaluate(
            """() => {
                const data=JSON.parse(document.getElementById('workstation-data').textContent);
                const rows=Array.from(document.querySelectorAll('.result-row')).map(el=>el.dataset.region);
                const map=new Map(data.region_index.map(row=>[row.region,row]));
                return {
                    topo:document.getElementById('ftopo').value,chrom:document.getElementById('fchr').value,
                    evidence:document.getElementById('fevidence').value,signal:document.getElementById('fsignal').value,
                    query:document.getElementById('fq').value,count:document.getElementById('result-count').innerText,
                    hash:new URLSearchParams(location.hash.slice(1)).get('topo'),focus:document.activeElement?.id,
                    allMatch:rows.every(region=>map.get(region).topology_class==='topology_multiple_exact_multiple')
                };
            }"""
        )
        expected_multi = counts["topology_multiple_exact_multiple"]
        passed = (
            action["topo"] == "topology_multiple_exact_multiple"
            and action["chrom"] == ""
            and action["evidence"] == action["signal"] == action["query"] == ""
            and number_from_count(action["count"]) == expected_multi
            and action["hash"] == "topology_multiple_exact_multiple"
            and action["focus"] == "ftopo"
            and action["allMatch"]
        )
        add_check(checks, f"{prefix}/topology_action_resets_filters", passed, action, expected_multi)

        page.locator("#fevidence").select_option("partial_only")
        page.locator("#fsignal").select_option("hidden")
        page.locator("#fq").fill("chr22:999999999")
        page.locator("[data-dimension-chrom]").click()
        page.wait_for_timeout(100)
        position = page.evaluate(
            """() => ({
                chrom:document.getElementById('fchr').value,topo:document.getElementById('ftopo').value,
                evidence:document.getElementById('fevidence').value,signal:document.getElementById('fsignal').value,
                query:document.getElementById('fq').value,count:document.getElementById('result-count').innerText,
                hash:new URLSearchParams(location.hash.slice(1)).get('chr'),focus:document.activeElement?.id,
                status:document.getElementById('genome-status').innerText
            })"""
        )
        leader = integrity["countLeader"]
        passed = (
            position["chrom"] == leader["chrom"]
            and position["topo"] == ""
            and position["evidence"] == position["signal"] == position["query"] == ""
            and number_from_count(position["count"]) == leader["W_tree"]
            and position["hash"] == leader["chrom"]
            and position["focus"] == "fchr"
            and "W_primary" in position["status"]
            and "W_tree" in position["status"]
        )
        add_check(checks, f"{prefix}/position_action_denominator", passed, position, leader)

    if page_name == "HCC1395.html" and viewport_name in {"desktop", "mobile"}:
        row = page.evaluate(
            """() => JSON.parse(document.getElementById('workstation-data').textContent).region_index
                .find(row=>row.topology_class==='topology_multiple_exact_multiple')"""
        )
        page.evaluate("region => selectRegion(region,null,false,'none')", row["region"])
        page.wait_for_timeout(100)
        target = SHOTS / f"{viewport_name}_hcc1395_region_dimensions.png"
        page.locator("#detail .region-dimensions").screenshot(
            path=str(target), animations="disabled", caret="hide", timeout=120_000
        )


def main() -> int:
    SHOTS.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, Any]] = []
    runs: list[dict[str, Any]] = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True, args=["--allow-file-access-from-files"])
        version = browser.version
        for page_name in PAGES:
            for viewport_name, viewport in VIEWPORTS.items():
                context = browser.new_context(viewport=viewport, device_scale_factor=1)
                page = context.new_page()
                record = {
                    "page": page_name,
                    "viewport": viewport_name,
                    "console_errors": [],
                    "page_errors": [],
                    "request_failures": [],
                }
                runtime_handlers(page, record)
                try:
                    page.goto((SOURCE / page_name).as_uri(), wait_until="load", timeout=120_000)
                    if page_name == "index.html":
                        page.wait_for_selector("#cohort-table tbody tr", timeout=30_000)
                    else:
                        page.wait_for_function(
                            "document.querySelectorAll('.chrom-button').length===22 && document.querySelectorAll('.result-row').length>0",
                            timeout=30_000,
                        )
                    page.wait_for_timeout(180)
                    record["screenshot"] = screenshot_first(page, page_name, viewport_name)
                    if page_name == "index.html":
                        audit_index(page, viewport_name, checks)
                    else:
                        audit_sample(page, page_name, viewport_name, checks)
                except Exception as exc:  # noqa: BLE001 - audit must record all page failures
                    record["fatal_error"] = str(exc)
                    add_check(checks, f"{page_name}/{viewport_name}/fatal", False, str(exc), "no exception")
                add_check(
                    checks,
                    f"{page_name}/{viewport_name}/runtime",
                    not record["console_errors"] and not record["page_errors"] and not record["request_failures"],
                    {
                        "console": record["console_errors"],
                        "page": record["page_errors"],
                        "request": record["request_failures"],
                    },
                    "no errors",
                )
                runs.append(record)
                context.close()
        browser.close()

    failed = [check for check in checks if check["status"] == "fail"]
    payload = {
        "schema_version": "1.0",
        "task_type": "B_comprehensive_validation",
        "partial": False,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input": str(SOURCE / "*.html"),
        "output": str(OUT),
        "browser": {"engine": "chromium", "version": version},
        "scope": {"pages": len(PAGES), "viewports": VIEWPORTS, "contexts": len(runs)},
        "summary": {
            "status": "pass" if not failed else "fail",
            "checks_total": len(checks),
            "checks_passed": len(checks) - len(failed),
            "checks_failed": len(failed),
        },
        "failed_checks": failed,
        "checks": checks,
        "runs": runs,
    }
    (OUT / "metrics.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": payload["summary"], "metrics": str(OUT / "metrics.json")}, ensure_ascii=False, indent=2))
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
