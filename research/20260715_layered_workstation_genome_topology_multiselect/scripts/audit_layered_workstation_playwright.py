#!/usr/bin/env python3
"""Audit the canonical layered workstation with local Chromium/Playwright.

This is a deterministic, read-only UI audit for the cohort index and seven
canonical sample pages.  It reconciles the rendered GRCh38 ideogram and index
distributions with current-v5 read-AF topology sidecars, exercises category
multi-selection, validates the sidecar 1.1 exhaustive edge census, opens real
ranking and zero-coverage region details, and captures desktop/mobile visual
evidence.  A fixed HCC1395 fixture guards full-union, selected-tree, and stored
preview consistency across panels.

Default output:
    ../qa/full/validation_receipt.json

Exit codes:
    0: every check passed
    1: one or more product checks failed
    2: audit configuration or Chromium startup failed
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

try:
    from playwright.sync_api import Browser, Page, sync_playwright
except ImportError as exc:  # Reported in the JSON receipt by main().
    Browser = Page = Any  # type: ignore[misc,assignment]
    sync_playwright = None  # type: ignore[assignment]
    PLAYWRIGHT_IMPORT_ERROR: Optional[str] = str(exc)
else:
    PLAYWRIGHT_IMPORT_ERROR = None


SCRIPT_PATH = Path(__file__).resolve()
TOPIC_ROOT = SCRIPT_PATH.parents[1]
REPO_ROOT = SCRIPT_PATH.parents[3]
DEFAULT_WORKSTATION_DIR = REPO_ROOT / "docs/methodology/_assets/layered_workstation"
DEFAULT_SIDECAR_DIR = TOPIC_ROOT / "data/current_v5_read_af_topology"
DEFAULT_OUTPUT_DIR = TOPIC_ROOT / "qa/full"

SAMPLES: Tuple[str, ...] = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)
VIEWPORTS: Dict[str, Dict[str, int]] = {
    "desktop": {"width": 1440, "height": 1000},
    "mobile": {"width": 390, "height": 844},
    "narrow": {"width": 320, "height": 800},
}
IDEOGRAM_MODES: Tuple[str, ...] = (
    "determinacy",
    "read-af-selection",
    "morphology",
    "evidence",
    "primary-hp",
    "n-ssnv",
    "cn-sidecar",
)
EXPECTED_OVERVIEW_IDS: Tuple[str, ...] = (
    "topology-count",
    "candidate-count",
    "determinacy",
    "read-af-selection",
    "morphology",
    "hp-h3",
    "region-size",
)
EXPECTED_UI_CONTRACT = "layered-workstation-v5-grch38-topology-multiselect-4"
EXPECTED_SIDECAR_SCHEMA = (
    "intersubmod.current_v5_read_af_topology_sample",
    "1.1.0",
)
EXPECTED_AGGREGATE_W_PRIMARY = 50_215
HCC1395_EDGE_CENSUS_FIXTURE: Dict[str, Any] = {
    "region": "chr10:87818272-87928739",
    "family": "1",
    "candidate_total": 74,
    "top_candidate_total": 1,
    "node_total": 16,
    "edge_total": 32,
    "stored_candidate_total": 32,
    "edge": ["H_RRAA", "H_ARAA"],
    "candidate_count": 2,
    "top_candidate_count": 1,
}
EDGE_CENSUS_ELIGIBLE_STATUSES = {
    "ranked",
    "missing_read_af",
    "recurrence_not_evaluable",
}
INDEX_SELECTION_BINS: Tuple[Tuple[str, str], ...] = (
    ("structural_exact_unique", "結構已 exact 唯一"),
    ("read_af_unique_first", "read-AF 唯一第一順位"),
    ("read_af_tied_same_topology", "並列第一 · 同一 Topo"),
    ("read_af_tied_different_topology", "並列第一 · 多 Topo"),
    ("read_af_unavailable", "read-AF 不可排序"),
    ("incomplete", "候選未完整"),
)
INDEX_MORPHOLOGY_BINS: Tuple[Tuple[str, str], ...] = (
    ("single_no_within_hp_relation", "單支／無 HP 內關係"),
    ("direct_chain", "直系鏈"),
    ("sister_branch", "旁系／分支"),
    ("direct_and_sister", "直系＋旁系"),
    ("unresolved", "形態未解"),
)


def utc_now() -> str:
    """Return a timezone-aware ISO timestamp."""

    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    """Hash a file without loading it all into memory."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_safe(value: Any) -> Any:
    """Convert audit values to JSON-safe primitives."""

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return str(value)


def add_check(
    host: Dict[str, Any],
    name: str,
    passed: bool,
    *,
    expected: Any,
    actual: Any,
) -> None:
    """Append one explicit assertion to a run or index receipt."""

    host.setdefault("checks", []).append(
        {
            "name": name,
            "status": "pass" if passed else "fail",
            "expected": json_safe(expected),
            "actual": json_safe(actual),
        }
    )


def load_json(path: Path) -> Dict[str, Any]:
    """Read one UTF-8 JSON object."""

    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return payload


def is_plain_int(value: Any) -> bool:
    """Return whether a JSON value is an integer rather than bool/float."""

    return isinstance(value, int) and not isinstance(value, bool)


def audit_full_edge_census_contract(
    sample: str,
    payload: Dict[str, Any],
) -> Dict[str, Any]:
    """Validate the sidecar 1.1 exhaustive edge-census contract.

    A census is required only for complete, successfully re-enumerated units
    with more than one exact candidate.  Incomplete or enumeration-mismatch
    units must not pretend to expose a complete edge universe.  The audit also
    proves that every durable read-AF top representative edge belongs to the
    exhaustive census and has positive top-candidate support.
    """

    violation_count = 0
    violations: List[Dict[str, Any]] = []

    def record(code: str, location: str, actual: Any) -> None:
        nonlocal violation_count
        violation_count += 1
        if len(violations) < 100:
            violations.append(
                {"code": code, "location": location, "actual": json_safe(actual)}
            )

    schema_actual = (payload.get("schema_name"), payload.get("schema_version"))
    if schema_actual != EXPECTED_SIDECAR_SCHEMA:
        record(
            "sidecar_schema_not_1_1",
            sample,
            {"expected": EXPECTED_SIDECAR_SCHEMA, "actual": schema_actual},
        )

    required_units = 0
    census_units = 0
    census_edge_rows = 0
    top_representative_edges = 0
    top_representative_edges_in_census = 0
    top_representative_edges_with_positive_top_count = 0

    for region in payload.get("regions") or []:
        region_id = str(region.get("region") or "<missing-region>")
        units = region.get("units") or {}
        if not isinstance(units, dict):
            record("units_not_object", region_id, units)
            continue
        for family, unit in units.items():
            location = f"{sample}/{region_id}/HP{family}"
            if not isinstance(unit, dict):
                record("unit_not_object", location, unit)
                continue
            status = unit.get("status")
            expected_total = unit.get("n_trees_expected")
            enumerated_total = unit.get("n_trees_enumerated")
            claims_complete_multitree = (
                status in EDGE_CENSUS_ELIGIBLE_STATUSES
                and is_plain_int(expected_total)
                and expected_total > 1
            )
            enumeration_matches = (
                is_plain_int(enumerated_total)
                and enumerated_total == expected_total
            )
            if claims_complete_multitree and not enumeration_matches:
                record(
                    "complete_status_enumeration_mismatch",
                    location,
                    {
                        "status": status,
                        "n_trees_expected": expected_total,
                        "n_trees_enumerated": enumerated_total,
                    },
                )
            eligible = claims_complete_multitree and enumeration_matches
            census = unit.get("full_edge_census")
            if eligible:
                required_units += 1
            if eligible and not isinstance(census, dict):
                record("eligible_unit_missing_full_edge_census", location, census)
            if census is None:
                representatives = unit.get("top_shape_representatives") or []
                for representative in representatives:
                    top_representative_edges += len(
                        {tuple(edge) for edge in representative.get("edges") or []}
                    )
                continue
            if not isinstance(census, dict):
                record("full_edge_census_not_object", location, census)
                continue
            census_units += 1
            expected_census_keys = {
                "complete",
                "candidate_total",
                "top_candidate_total",
                "node_total",
                "edge_total",
                "edge_rows",
            }
            if set(census) != expected_census_keys:
                record(
                    "full_edge_census_field_set_mismatch",
                    location,
                    {
                        "expected": sorted(expected_census_keys),
                        "actual": sorted(census),
                    },
                )
            if not eligible:
                record(
                    "full_edge_census_present_for_ineligible_unit",
                    location,
                    {
                        "status": status,
                        "n_trees_expected": expected_total,
                        "n_trees_enumerated": enumerated_total,
                    },
                )

            complete = census.get("complete")
            candidate_total = census.get("candidate_total")
            top_candidate_total = census.get("top_candidate_total")
            node_total = census.get("node_total")
            edge_total = census.get("edge_total")
            edge_rows = census.get("edge_rows")
            scalar_contract = {
                "complete": complete is True,
                "candidate_total": is_plain_int(candidate_total)
                and candidate_total > 1,
                "top_candidate_total": is_plain_int(top_candidate_total)
                and top_candidate_total >= 0
                and is_plain_int(candidate_total)
                and top_candidate_total <= candidate_total,
                "node_total": is_plain_int(node_total) and node_total >= 1,
                "edge_total": is_plain_int(edge_total) and edge_total >= 1,
                "edge_rows": isinstance(edge_rows, list),
            }
            if not all(scalar_contract.values()):
                record("full_edge_census_scalar_contract", location, scalar_contract)
                continue
            if candidate_total != expected_total or candidate_total != enumerated_total:
                record(
                    "candidate_total_mismatch",
                    location,
                    {
                        "census": candidate_total,
                        "expected": expected_total,
                        "enumerated": enumerated_total,
                    },
                )
            expected_top_total = (
                unit.get("n_top_exact") if status == "ranked" else 0
            )
            if top_candidate_total != expected_top_total:
                record(
                    "top_candidate_total_mismatch",
                    location,
                    {
                        "census": top_candidate_total,
                        "expected": expected_top_total,
                        "status": status,
                    },
                )
            if edge_total != len(edge_rows):
                record(
                    "edge_total_mismatch",
                    location,
                    {"edge_total": edge_total, "edge_rows": len(edge_rows)},
                )
            census_edge_rows += len(edge_rows)

            row_map: Dict[Tuple[str, str], Tuple[int, int]] = {}
            normalized_keys: List[Tuple[str, str]] = []
            nodes = set()
            for row_index, row in enumerate(edge_rows):
                row_location = f"{location}/edge_rows/{row_index}"
                if not isinstance(row, list) or len(row) != 4:
                    record("edge_row_not_four_item_array", row_location, row)
                    continue
                parent, child, count, top_count = row
                if not (
                    isinstance(parent, str)
                    and parent
                    and isinstance(child, str)
                    and child
                    and parent != child
                    and is_plain_int(count)
                    and 1 <= count <= candidate_total
                    and is_plain_int(top_count)
                    and 0 <= top_count <= top_candidate_total
                    and top_count <= count
                ):
                    record("edge_row_value_contract", row_location, row)
                    continue
                key = (parent, child)
                if key in row_map:
                    record("duplicate_edge_row", row_location, key)
                row_map[key] = (count, top_count)
                normalized_keys.append(key)
                nodes.update(key)
            if normalized_keys != sorted(normalized_keys):
                record("edge_rows_not_lexicographically_sorted", location, normalized_keys)
            if len(row_map) != edge_total:
                record(
                    "edge_rows_not_unique_or_complete",
                    location,
                    {"unique_rows": len(row_map), "edge_total": edge_total},
                )
            if len(nodes) != node_total:
                record(
                    "node_total_mismatch",
                    location,
                    {"nodes_from_edges": len(nodes), "node_total": node_total},
                )

            representative_edges = {
                tuple(edge)
                for representative in unit.get("top_shape_representatives") or []
                for edge in representative.get("edges") or []
            }
            top_representative_edges += len(representative_edges)
            for edge in sorted(representative_edges):
                census_counts = row_map.get(edge)
                if census_counts is None:
                    record("top_representative_edge_missing_from_census", location, edge)
                    continue
                top_representative_edges_in_census += 1
                if census_counts[1] < 1:
                    record(
                        "top_representative_edge_has_zero_top_count",
                        location,
                        {"edge": edge, "counts": census_counts},
                    )
                else:
                    top_representative_edges_with_positive_top_count += 1

    return {
        "status": "pass" if violation_count == 0 else "fail",
        "schema_expected": list(EXPECTED_SIDECAR_SCHEMA),
        "schema_actual": list(schema_actual),
        "required_units": required_units,
        "census_units": census_units,
        "census_edge_rows": census_edge_rows,
        "top_representative_edges": top_representative_edges,
        "top_representative_edges_in_census": top_representative_edges_in_census,
        "top_representative_edges_with_positive_top_count": (
            top_representative_edges_with_positive_top_count
        ),
        "violation_count": violation_count,
        "violations_sample": violations,
    }


def extract_hcc1395_edge_census_fixture(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Extract the fixed PTEN-region HP1 census used by the UI regression."""

    expected = HCC1395_EDGE_CENSUS_FIXTURE
    region = next(
        (
            row
            for row in payload.get("regions") or []
            if row.get("region") == expected["region"]
        ),
        None,
    )
    unit = (region.get("units") or {}).get(expected["family"]) if region else None
    census = unit.get("full_edge_census") if isinstance(unit, dict) else None
    edge_row = None
    if isinstance(census, dict):
        edge_row = next(
            (
                row
                for row in census.get("edge_rows") or []
                if isinstance(row, list)
                and len(row) == 4
                and row[:2] == expected["edge"]
            ),
            None,
        )
    return {
        "region_found": region is not None,
        "family_found": isinstance(unit, dict),
        "census": json_safe(census),
        "edge_row": json_safe(edge_row),
    }


def choose_ranking_fixture(regions: Iterable[Dict[str, Any]]) -> Dict[str, Any]:
    """Choose a complete ranked region with a table and drawable top tree.

    A representative whose original candidate index exceeds the historical
    display cap of 32 is preferred because it verifies that the independent
    current-v5 ranking view does not silently lose that result.
    """

    candidates: List[Dict[str, Any]] = []
    for region in regions:
        if region.get("selection_class") not in {
            "read_af_unique_first",
            "read_af_tied_same_topology",
            "read_af_tied_different_topology",
        }:
            continue
        for family, unit in (region.get("units") or {}).items():
            preview = unit.get("ranked_preview") or []
            representatives = unit.get("top_shape_representatives") or []
            if unit.get("status") != "ranked" or not preview or not representatives:
                continue
            top_indices = [int(item.get("candidate_index", 0)) for item in representatives]
            candidates.append(
                {
                    "region": region["region"],
                    "family": str(family),
                    "selection_class": region.get("selection_class"),
                    "exact_top_class": unit.get("exact_top_class"),
                    "n_trees_expected": int(unit.get("n_trees_expected", 0)),
                    "top_candidate_indices": top_indices,
                    "outside_original_cap_32": any(index > 32 for index in top_indices),
                }
            )
    if not candidates:
        raise ValueError("No ranked region has both ranked_preview and top_shape_representatives")
    candidates.sort(
        key=lambda item: (
            not item["outside_original_cap_32"],
            item["n_trees_expected"],
            item["region"],
            item["family"],
        )
    )
    return candidates[0]


def load_sidecar_expectations(
    sidecar_dir: Path,
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Any]]:
    """Load, hash, and validate sidecars plus their aggregate index contract."""

    index_path = sidecar_dir / "current_v5_read_af_topology.index.json"
    index = load_json(index_path)
    index_sha256 = sha256_file(index_path)
    index_by_sample = {row["sample"]: row for row in index.get("samples", [])}
    missing = sorted(set(SAMPLES) - set(index_by_sample))
    if missing:
        raise ValueError(f"Sidecar index is missing samples: {', '.join(missing)}")

    expectations: Dict[str, Dict[str, Any]] = {}
    for sample in SAMPLES:
        path = sidecar_dir / f"{sample}.current_v5_read_af_topology.json"
        if not path.is_file():
            raise FileNotFoundError(f"Missing current-v5 sidecar: {path}")
        digest = sha256_file(path)
        index_row = index_by_sample[sample]
        if digest != index_row.get("output_sha256"):
            raise ValueError(
                f"Sidecar SHA-256 mismatch for {sample}: {digest} != "
                f"{index_row.get('output_sha256')}"
            )
        payload = load_json(path)
        if payload.get("sample") != sample:
            raise ValueError(f"Sidecar sample mismatch: {path}")
        summary = payload.get("summary") or {}
        if not summary.get("all_checks_pass"):
            raise ValueError(f"Sidecar builder checks are not all pass: {path}")
        fixture = choose_ranking_fixture(payload.get("regions") or [])
        edge_census_contract = audit_full_edge_census_contract(sample, payload)
        expectations[sample] = {
            "path": str(path.resolve()),
            "sha256": digest,
            "schema_name": payload.get("schema_name"),
            "schema_version": payload.get("schema_version"),
            "W_tree": int(summary["W_tree"]),
            "W_primary": int(summary["W_primary"]),
            "no_primary": int(summary["no_primary"]),
            "selection_classes": {
                key: int(value) for key, value in summary["selection_classes"].items()
            },
            "morphology_classes": {
                key: int(value) for key, value in summary["morphology_classes"].items()
            },
            "ranking_fixture": fixture,
            "full_edge_census_contract": edge_census_contract,
            "hcc1395_edge_census_fixture": (
                extract_hcc1395_edge_census_fixture(payload)
                if sample == "HCC1395"
                else None
            ),
            "zero_coverage_fixture": None,
        }
    selection_aggregate = {
        key: sum(item["selection_classes"][key] for item in expectations.values())
        for key, _label in INDEX_SELECTION_BINS
    }
    morphology_aggregate = {
        key: sum(item["morphology_classes"][key] for item in expectations.values())
        for key, _label in INDEX_MORPHOLOGY_BINS
    }
    aggregate_w_primary = sum(item["W_primary"] for item in expectations.values())
    index_contract = {
        "path": str(index_path.resolve()),
        "sha256": index_sha256,
        "schema_name": index.get("schema_name"),
        "schema_version": index.get("schema_version"),
        "dataset_count": int(index.get("dataset_count", 0)),
        "W_primary": aggregate_w_primary,
        "expected_canonical_W_primary": EXPECTED_AGGREGATE_W_PRIMARY,
        "selection_classes": selection_aggregate,
        "morphology_classes": morphology_aggregate,
        "selection_sum": sum(selection_aggregate.values()),
        "morphology_sum": sum(morphology_aggregate.values()),
        "sidecar_schema_versions": sorted(
            {item["schema_version"] for item in expectations.values()}
        ),
        "full_edge_census_contracts_pass": all(
            item["full_edge_census_contract"]["status"] == "pass"
            for item in expectations.values()
        ),
    }
    return expectations, index_contract


def index_audit(
    browser: Browser,
    workstation_dir: Path,
    output_dir: Path,
    viewport_name: str,
    sidecar_index_contract: Dict[str, Any],
    timeout_ms: int,
) -> Dict[str, Any]:
    """Audit the canonical cohort index at one desktop/mobile viewport."""

    started = time.monotonic()
    path = (workstation_dir / "index.html").resolve()
    result: Dict[str, Any] = {
        "kind": "index",
        "sample": None,
        "viewport": viewport_name,
        "viewport_size": VIEWPORTS[viewport_name],
        "input_path": str(path),
        "url": path.as_uri(),
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
    page.on(
        "console",
        lambda message: result["console_errors"].append(message.text)
        if message.type == "error"
        else None,
    )
    page.on("pageerror", lambda error: result["page_errors"].append(str(error)))
    try:
        page.goto(path.as_uri(), wait_until="load", timeout=timeout_ms)
        page.wait_for_load_state("networkidle", timeout=timeout_ms)
        page.add_style_tag(content="html{scroll-behavior:auto!important}")
        add_check(
            result,
            "index_document_load",
            True,
            expected="index reaches load and networkidle",
            actual={"ready_state": page.evaluate("document.readyState")},
        )
        links = page.evaluate(
            """samples => Object.fromEntries(samples.map(sample => {
                const expected = sample + '.html';
                const link = [...document.querySelectorAll('a[href]')].find(anchor => {
                    try { return new URL(anchor.getAttribute('href'), location.href).pathname.endsWith('/' + expected); }
                    catch (_) { return false; }
                });
                return [sample, link ? link.getAttribute('href') : null];
            }))""",
            list(SAMPLES),
        )
        missing = [sample for sample, href in links.items() if not href]
        add_check(
            result,
            "index_links_all_seven_samples",
            not missing,
            expected={"samples": list(SAMPLES)},
            actual={"links": links, "missing": missing},
        )

        meta_sha256 = page.locator(
            'meta[name="intersubmod-read-af-topology-index-sha256"]'
        ).get_attribute("content")
        add_check(
            result,
            "index_sidecar_index_sha256_binding",
            meta_sha256 == sidecar_index_contract["sha256"],
            expected={"sha256": sidecar_index_contract["sha256"]},
            actual={"sha256": meta_sha256},
        )

        index_ui_contract = page.locator(
            'meta[name="intersubmod-ui-contract"]'
        ).get_attribute("content")
        add_check(
            result,
            "index_ui_contract",
            index_ui_contract == EXPECTED_UI_CONTRACT,
            expected={"ui_contract": EXPECTED_UI_CONTRACT},
            actual={"ui_contract": index_ui_contract},
        )

        distribution = page.evaluate(
            """() => {
                const visible = element => {
                    const rect = element.getBoundingClientRect();
                    const style = getComputedStyle(element);
                    return style.display !== 'none' && style.visibility !== 'hidden' &&
                        rect.width > 0 && rect.height > 0;
                };
                const parseCount = text => {
                    const match = String(text || '').match(/[0-9][0-9,]*/);
                    return match ? Number(match[0].replaceAll(',', '')) : null;
                };
                const parseWPrimary = text => {
                    const match = String(text || '').match(/W_primary\s+([0-9][0-9,]*)/i);
                    return match ? Number(match[1].replaceAll(',', '')) : null;
                };
                const title = document.getElementById('read-af-morphology-title');
                const section = title?.closest('section') || null;
                const cards = [...(section?.querySelectorAll('.distribution-card') || [])]
                    .map((card, index) => {
                        const total = card.querySelector('.distribution-head > b');
                        const rows = [...card.querySelectorAll('.distribution-legend li')]
                            .map(row => {
                                const label = row.querySelector('strong');
                                const value = row.querySelector('b');
                                return {
                                    label: (label?.innerText || '').trim(),
                                    count: parseCount(value?.innerText || ''),
                                    row_visible: visible(row),
                                    label_visible: !!label && visible(label),
                                    count_visible: !!value && visible(value)
                                };
                            });
                        return {
                            index,
                            visible: visible(card),
                            layer: (card.querySelector('.distribution-head span')?.innerText || '').trim(),
                            total_text: (total?.innerText || '').trim(),
                            total: parseWPrimary(total?.innerText || ''),
                            total_visible: !!total && visible(total),
                            rows,
                            counts: rows.map(row => row.count),
                            labels: rows.map(row => row.label),
                            sum: rows.reduce((sum, row) => sum + Number(row.count || 0), 0)
                        };
                    });
                return {
                    title_exists: !!title,
                    title_visible: !!title && visible(title),
                    title_text: (title?.innerText || '').trim(),
                    cards
                };
            }"""
        )
        add_check(
            result,
            "index_interpretation_section_and_two_cards_visible",
            distribution["title_exists"]
            and distribution["title_visible"]
            and len(distribution["cards"]) == 2
            and all(card["visible"] for card in distribution["cards"]),
            expected={
                "title": "#read-af-morphology-title visible",
                "distribution_cards": 2,
                "cards_visible": True,
            },
            actual=distribution,
        )

        expected_cards = [
            {
                "labels": [label for _key, label in INDEX_SELECTION_BINS],
                "counts": [
                    sidecar_index_contract["selection_classes"][key]
                    for key, _label in INDEX_SELECTION_BINS
                ],
            },
            {
                "labels": [label for _key, label in INDEX_MORPHOLOGY_BINS],
                "counts": [
                    sidecar_index_contract["morphology_classes"][key]
                    for key, _label in INDEX_MORPHOLOGY_BINS
                ],
            },
        ]
        card_counts_pass = len(distribution["cards"]) == 2
        if card_counts_pass:
            for card, expected_card in zip(distribution["cards"], expected_cards):
                card_counts_pass = card_counts_pass and (
                    card["labels"] == expected_card["labels"]
                    and card["counts"] == expected_card["counts"]
                    and card["total"] == sidecar_index_contract["W_primary"]
                    and card["sum"] == sidecar_index_contract["W_primary"]
                    and card["total_visible"]
                    and all(
                        row["row_visible"]
                        and row["label_visible"]
                        and row["count_visible"]
                        for row in card["rows"]
                    )
                )
        card_counts_pass = card_counts_pass and (
            sidecar_index_contract["W_primary"] == EXPECTED_AGGREGATE_W_PRIMARY
            and sidecar_index_contract["selection_sum"] == EXPECTED_AGGREGATE_W_PRIMARY
            and sidecar_index_contract["morphology_sum"] == EXPECTED_AGGREGATE_W_PRIMARY
        )
        add_check(
            result,
            "index_distribution_counts_match_seven_sidecars_and_conserve_W_primary",
            card_counts_pass,
            expected={
                "W_primary": EXPECTED_AGGREGATE_W_PRIMARY,
                "selection": expected_cards[0],
                "morphology": expected_cards[1],
                "all_counts_visible": True,
                "each_card_sum": EXPECTED_AGGREGATE_W_PRIMARY,
            },
            actual={
                "sidecar_index_contract": sidecar_index_contract,
                "cards": distribution["cards"],
            },
        )

        overflow = body_overflow_metrics(page)
        add_check(
            result,
            "index_body_horizontal_overflow",
            float(overflow["overflow_px"]) <= 1,
            expected={"overflow_px_max": 1},
            actual=overflow,
        )

        json_links = page.evaluate(
            """() => {
                const links = [...document.querySelectorAll('a[href]')]
                    .filter(link => (link.getAttribute('href') || '').toLowerCase().includes('.json'));
                return {
                    count: links.length,
                    inside_details: links.filter(link => !!link.closest('details')).length,
                    outside_details: links.filter(link => !link.closest('details')).map(link => ({
                        href: link.getAttribute('href'),
                        text: (link.innerText || '').trim()
                    }))
                };
            }"""
        )
        add_check(
            result,
            "index_all_json_links_inside_details",
            json_links["count"] > 0
            and json_links["inside_details"] == json_links["count"]
            and not json_links["outside_details"],
            expected={"json_links_minimum": 1, "all_inside_details": True},
            actual=json_links,
        )

        contrast_failures = text_contrast_failures(page)
        add_check(
            result,
            "index_visible_text_wcag_aa_contrast",
            not contrast_failures,
            expected={"failures": 0, "normal_text_ratio": 4.5, "large_text_ratio": 3.0},
            actual={"failures": len(contrast_failures), "examples": contrast_failures},
        )

        scroll_element_to_viewport_top(
            page,
            'section[aria-labelledby="read-af-morphology-title"]',
        )
        page.wait_for_timeout(180)
        screenshot_path = output_dir / f"index__{viewport_name}__interpretation.png"
        result["screenshots"]["interpretation"] = screenshot(page, screenshot_path)

        page.wait_for_timeout(100)
        add_check(
            result,
            "index_console_and_page_errors",
            not result["console_errors"] and not result["page_errors"],
            expected={"console_errors": 0, "page_errors": 0},
            actual={
                "console_errors": result["console_errors"],
                "page_errors": result["page_errors"],
            },
        )
    except Exception as exc:
        add_check(
            result,
            "index_document_load",
            False,
            expected="index reaches load and networkidle",
            actual={"error": str(exc)},
        )
    finally:
        result["duration_seconds"] = round(time.monotonic() - started, 3)
        context.close()
    return result


def page_overview_metrics(page: Page) -> Dict[str, Any]:
    """Read the seven overview panels and their denominator closure."""

    return page.evaluate(
        """() => {
            const data = JSON.parse(document.getElementById('workstation-data').textContent);
            const panels = [...document.querySelectorAll('[data-overview-panel]')].map(panel => {
                const bins = [...panel.querySelectorAll('[data-overview-bin]')].map(bin => {
                    const track = bin.querySelector('.overview-bar-track');
                    const fill = bin.querySelector('[data-overview-bar]');
                    const trackRect = track?.getBoundingClientRect();
                    const fillRect = fill?.getBoundingClientRect();
                    return {
                        key: bin.dataset.binKey,
                        count: Number(bin.dataset.count),
                        denominator: Number(bin.dataset.denominator),
                        style_width_pct: Number.parseFloat(fill?.style.width || '0'),
                        measured_width_pct: trackRect?.width ? 100 * (fillRect?.width || 0) / trackRect.width : 0,
                        printed: bin.querySelector('.overview-bin-value')?.textContent || ''
                    };
                });
                return {
                    id: panel.dataset.overviewPanel,
                    denominator_key: panel.dataset.denominatorKey,
                    total: Number(panel.dataset.total),
                    sum: bins.reduce((total, bin) => total + bin.count, 0),
                    bins
                };
            });
            return {
                panels,
                ids: panels.map(panel => panel.id),
                W_tree: Number(data.canonical_sample.W_tree),
                W_primary: Number(data.canonical_sample.W_primary),
                region_index_count: data.region_index.length
            };
        }"""
    )


def mode_metrics(page: Page, mode: str) -> Dict[str, Any]:
    """Switch one mode and return its rendered category histogram."""

    page.locator(f'[data-ideogram-mode="{mode}"]').click()
    return page.evaluate(
        """requested => {
            const marks = [...document.querySelectorAll('.ideogram-mark')];
            const histogram = {};
            marks.forEach(mark => {
                const key = mark.dataset.modeValue || '';
                histogram[key] = (histogram[key] || 0) + 1;
            });
            const legend = [...document.querySelectorAll('[data-legend-key]')].map(button => ({
                key: button.dataset.legendKey,
                aria_label: button.getAttribute('aria-label') || '',
                pressed: button.getAttribute('aria-pressed')
            }));
            const legendKeys = legend.map(item => item.key);
            const countLabelsExact = legend.every(item => {
                const count = histogram[item.key] || 0;
                return item.aria_label.includes(Number(count).toLocaleString());
            });
            const signatures = {};
            Object.keys(histogram).filter(Boolean).forEach(key => {
                const mark = marks.find(item => item.dataset.modeValue === key);
                const button = document.querySelector('[data-legend-key="' + CSS.escape(key) + '"]');
                const swatch = button?.querySelector('.ideogram-legend-swatch');
                const markStyle = mark ? getComputedStyle(mark) : null;
                const swatchStyle = swatch ? getComputedStyle(swatch) : null;
                signatures[key] = {
                    stroke: markStyle?.stroke || '',
                    stroke_width: markStyle?.strokeWidth || '',
                    stroke_dasharray: markStyle?.strokeDasharray || '',
                    swatch_background_image: swatchStyle?.backgroundImage || 'none'
                };
            });
            return {
                requested,
                active_mode: document.getElementById('genome-ideogram').dataset.mode,
                active_mode_buttons: document.querySelectorAll('[data-ideogram-mode][aria-pressed="true"]').length,
                marks: marks.length,
                histogram,
                histogram_sum: Object.values(histogram).reduce((total, count) => total + count, 0),
                unmapped: histogram[''] || 0,
                legend,
                legend_keys_match_histogram: JSON.stringify([...legendKeys].sort()) ===
                    JSON.stringify(Object.keys(histogram).filter(Boolean).sort()),
                legend_count_labels_exact: countLabelsExact
                ,signatures
            };
        }""",
        mode,
    )


def selection_state(page: Page) -> Dict[str, Any]:
    """Read current ideogram category selection and visibility."""

    return page.evaluate(
        """() => ({
            mode: document.getElementById('genome-ideogram').dataset.mode,
            selected: [...document.querySelectorAll('[data-legend-key][aria-pressed="true"]')]
                .map(button => button.dataset.legendKey).sort(),
            visible: document.querySelectorAll('.ideogram-mark:not(.dimmed)').length,
            dimmed: document.querySelectorAll('.ideogram-mark.dimmed').length,
            all_pressed: document.querySelector('[data-legend-clear]')?.getAttribute('aria-pressed') || null
        })"""
    )


def body_overflow_metrics(page: Page) -> Dict[str, Any]:
    """Measure page-level overflow; local scrollers remain allowed."""

    return page.evaluate(
        """() => {
            const root = document.documentElement;
            const body = document.body;
            const overflow = Math.max(root.scrollWidth, body ? body.scrollWidth : 0) - innerWidth;
            return {
                viewport_width: innerWidth,
                root_scroll_width: root.scrollWidth,
                body_scroll_width: body ? body.scrollWidth : null,
                overflow_px: Math.max(0, Math.round(overflow * 10) / 10)
            };
        }"""
    )


def touch_target_metrics(page: Page) -> Dict[str, Any]:
    """Measure the mobile genome controls that users repeatedly tap."""

    return page.evaluate(
        """() => {
            const selectors = [
                '[data-ideogram-mode]',
                '[data-legend-key]',
                '[data-legend-clear]',
                '[data-overview-panel="read-af-selection"] [data-overview-bin]',
                '.chrom-button'
            ];
            const controls = selectors.flatMap(selector =>
                [...document.querySelectorAll(selector)].map(element => {
                    const rect = element.getBoundingClientRect();
                    return {
                        selector,
                        tag: element.tagName.toLowerCase(),
                        key: element.dataset.ideogramMode || element.dataset.legendKey ||
                            element.dataset.binKey || element.dataset.chrom || element.id || null,
                        width: Math.round(rect.width * 10) / 10,
                        height: Math.round(rect.height * 10) / 10
                    };
                })
            );
            return {
                controls,
                minimum_height: controls.length ? Math.min(...controls.map(item => item.height)) : 0,
                below_44px: controls.filter(item => item.height < 44)
            };
        }"""
    )


def find_zero_coverage_region(page: Page) -> Optional[Dict[str, Any]]:
    """Find a primary HP site whose stored REF and ALT evidence are both zero."""

    return page.evaluate(
        """() => {
            for (const script of document.querySelectorAll('script[id^="chunk-"]')) {
                const text = script.textContent || '';
                if (!text.includes('[0,0]')) continue;
                const regions = JSON.parse(text);
                for (const region of regions) {
                    for (const line of (region.lineages || [])) {
                        if (!line.is_primary_lineage) continue;
                        for (const [position, counts] of Object.entries(line.obs_col_coverage || {})) {
                            if (Number(counts?.[0] || 0) === 0 && Number(counts?.[1] || 0) === 0) {
                                return {
                                    region: region.region,
                                    family: String(line.family),
                                    position: String(position)
                                };
                            }
                        }
                    }
                }
            }
            return null;
        }"""
    )


def select_region(page: Page, region: str, timeout_ms: int) -> None:
    """Open an exact region through the page's public interaction function."""

    available = page.evaluate("() => typeof selectRegion === 'function'")
    if not available:
        raise RuntimeError("Global selectRegion(region, ...) function is unavailable")
    page.evaluate("region => selectRegion(region, null, false, 'none')", region)
    page.wait_for_function(
        """region => {
            const heading = document.querySelector('#detail h3');
            return !!heading && (heading.textContent || '').includes(region);
        }""",
        arg=region,
        timeout=timeout_ms,
    )


def ranking_detail_metrics(page: Page, region: str) -> Dict[str, Any]:
    """Read visible ranking evidence for a selected real region."""

    return page.evaluate(
        """region => {
            const detail = document.getElementById('detail');
            const text = detail ? (detail.innerText || '') : '';
            const rankedTreeViews = [...(detail?.querySelectorAll('.read-af-top-tree') || [])]
                .map((tree, treeIndex) => {
                    const scroller = tree.querySelector('.network-scroll');
                    const svg = scroller?.querySelector('svg');
                    const scrollerRect = scroller?.getBoundingClientRect();
                    const nodes = [...(svg?.querySelectorAll('circle, rect') || [])].map(node => {
                        const rect = node.getBoundingClientRect();
                        const visible = !!scrollerRect && rect.width > 0 && rect.height > 0 &&
                            rect.right > scrollerRect.left && rect.left < scrollerRect.right &&
                            rect.bottom > scrollerRect.top && rect.top < scrollerRect.bottom;
                        return {
                            tag: node.tagName.toLowerCase(),
                            x: Math.round(rect.x * 10) / 10,
                            y: Math.round(rect.y * 10) / 10,
                            width: Math.round(rect.width * 10) / 10,
                            height: Math.round(rect.height * 10) / 10,
                            intersects_scroller_viewport: visible
                        };
                    });
                    return {
                        tree_index: treeIndex,
                        scroller: scrollerRect ? {
                            x: Math.round(scrollerRect.x * 10) / 10,
                            y: Math.round(scrollerRect.y * 10) / 10,
                            width: Math.round(scrollerRect.width * 10) / 10,
                            height: Math.round(scrollerRect.height * 10) / 10,
                            client_width: scroller.clientWidth,
                            scroll_width: scroller.scrollWidth,
                            initial_scroll_left: scroller.scrollLeft
                        } : null,
                        svg: svg ? {
                            rendered_width: Math.round(svg.getBoundingClientRect().width * 10) / 10,
                            view_box: svg.getAttribute('viewBox')
                        } : null,
                        nodes,
                        node_count: nodes.length,
                        visible_nodes: nodes.filter(node => node.intersects_scroller_viewport).length
                    };
                });
            return {
                region,
                heading: detail?.querySelector('h3')?.innerText || null,
                ranking_cards: detail?.querySelectorAll('.read-af-card').length || 0,
                ranking_tables: detail?.querySelectorAll('.ranking-table').length || 0,
                ranked_trees: detail?.querySelectorAll('.read-af-top-tree').length || 0,
                ranked_networks: detail?.querySelectorAll('.read-af-top-tree svg').length || 0,
                ranked_tree_views: rankedTreeViews,
                ranked_trees_with_visible_nodes: rankedTreeViews
                    .filter(view => view.visible_nodes >= 1).length,
                exact_fraction_text: text.includes('Exact Fraction'),
                probability_boundary_text: /非 posterior|not posterior/i.test(text)
            };
        }""",
        region,
    )


def hcc1395_edge_census_cross_panel_metrics(page: Page) -> Dict[str, Any]:
    """Inspect the fixed full-union/stored-preview/read-AF cross-panel fixture."""

    fixture = HCC1395_EDGE_CENSUS_FIXTURE
    page.locator("#detail details").evaluate_all(
        "elements => elements.forEach(element => { element.open = true; })"
    )
    page.wait_for_timeout(120)
    return page.evaluate(
        r"""fixture => {
            const detail = document.getElementById('detail');
            const edgeKey = fixture.edge.join('>');
            const visible = element => {
                if (!element) return false;
                const style = getComputedStyle(element);
                const rect = element.getBoundingClientRect();
                return style.display !== 'none' && style.visibility !== 'hidden' &&
                    Number(style.opacity || 1) > 0 &&
                    (rect.width > 0 || rect.height > 0);
            };
            const numberData = (element, key) => {
                const value = element?.dataset?.[key];
                return value == null || value === '' ? null : Number(value);
            };
            const fullNetworks = [...(detail?.querySelectorAll(
                '.full-edge-network[data-edge-scope="complete"]'
            ) || [])];
            const full = fullNetworks.find(network =>
                [...network.querySelectorAll('.network-edge[data-edge-key]')]
                    .some(edge => edge.dataset.edgeKey === edgeKey)
            ) || null;
            const edges = [...(full?.querySelectorAll(
                '.network-edge[data-edge-key]'
            ) || [])];
            const targetEdge = edges.find(edge => edge.dataset.edgeKey === edgeKey) || null;
            const overlay = [...(full?.querySelectorAll(
                '.network-edge-selected-overlay[data-edge-key]'
            ) || [])].find(edge => edge.dataset.edgeKey === edgeKey) || null;
            const edgeRecords = edges.map(edge => ({
                key: edge.dataset.edgeKey || null,
                candidate_count: numberData(edge, 'candidateCount'),
                candidate_total: numberData(edge, 'candidateTotal'),
                top_count: numberData(edge, 'topCount'),
                top_total: numberData(edge, 'topTotal'),
                status: edge.dataset.edgeStatus || null,
                selected: edge.classList.contains('edge-selected'),
                visible: visible(edge)
            }));
            const fullCard = full?.closest('.network-card') || full;
            const fullText = (fullCard?.innerText || '').trim();
            const tables = [...(detail?.querySelectorAll('.edge-census-table') || [])];
            const table = tables.find(element => {
                const text = element.innerText || '';
                return text.includes('RRAA') && text.includes('ARAA');
            }) || tables[0] || null;
            const scopes = [...(detail?.querySelectorAll(
                '.stored-preview-scope[data-stored-count][data-candidate-total]'
            ) || [])];
            const storedScope = scopes.find(element =>
                numberData(element, 'storedCount') === fixture.stored_candidate_total &&
                numberData(element, 'candidateTotal') === fixture.candidate_total
            ) || null;
            const topEdges = [...(detail?.querySelectorAll(
                '.read-af-top-tree .network-edge[data-edge-key]'
            ) || [])];
            const topEdge = topEdges.find(edge => edge.dataset.edgeKey === edgeKey) || null;
            return {
                fixture,
                detail_heading: detail?.querySelector('h3')?.innerText || null,
                full_network_count: fullNetworks.length,
                matching_full_network_found: !!full,
                full_edge_count: edges.length,
                full_unique_edge_count: new Set(edgeRecords.map(edge => edge.key)).size,
                full_edge_records: edgeRecords,
                all_full_edges_have_expected_totals: edgeRecords.length > 0 &&
                    edgeRecords.every(edge =>
                        edge.candidate_total === fixture.candidate_total &&
                        edge.top_total === fixture.top_candidate_total
                    ),
                all_full_edge_statuses_match_counts: edgeRecords.length > 0 &&
                    edgeRecords.every(edge =>
                        edge.status === (
                            edge.candidate_count === fixture.candidate_total
                                ? 'forced'
                                : 'variable'
                        )
                    ),
                all_selected_classes_match_positive_top_count: edgeRecords.length > 0 &&
                    edgeRecords.every(edge => edge.selected === (edge.top_count > 0)),
                target_edge: targetEdge ? {
                    key: targetEdge.dataset.edgeKey,
                    candidate_count: numberData(targetEdge, 'candidateCount'),
                    candidate_total: numberData(targetEdge, 'candidateTotal'),
                    top_count: numberData(targetEdge, 'topCount'),
                    top_total: numberData(targetEdge, 'topTotal'),
                    status: targetEdge.dataset.edgeStatus || null,
                    selected: targetEdge.classList.contains('edge-selected'),
                    visible: visible(targetEdge)
                } : null,
                selected_overlay: overlay ? {
                    key: overlay.dataset.edgeKey,
                    visible: visible(overlay)
                } : null,
                edge_census_table_count: tables.length,
                edge_census_table_visible: visible(table),
                stored_preview_scope: storedScope ? {
                    stored_count: numberData(storedScope, 'storedCount'),
                    candidate_total: numberData(storedScope, 'candidateTotal'),
                    text: (storedScope.innerText || '').trim(),
                    visible: visible(storedScope)
                } : null,
                read_af_top_edge: topEdge ? {
                    key: topEdge.dataset.edgeKey,
                    visible: visible(topEdge)
                } : null,
                full_union_boundary_text: /聯集|union/i.test(fullText),
                not_single_tree_boundary_text: /不是單棵|非單棵|not\s+(?:a\s+)?single/i.test(fullText),
                non_probability_boundary_text: /非[^\n]{0,20}(?:機率|概率|probab)|不是[^\n]{0,20}(?:機率|概率)|not[^\n]{0,20}probab/i.test(fullText),
                full_card_text: fullText
            };
        }""",
        fixture,
    )


def zero_coverage_detail_metrics(page: Page, fixture: Dict[str, Any]) -> Dict[str, Any]:
    """Confirm that primary HP 0/0 evidence renders as N/A, never as 0%."""

    return page.evaluate(
        """fixture => {
            const detail = document.getElementById('detail');
            const candidates = [...detail.querySelectorAll('.family-af-site.na, td')]
                .filter(element => (element.innerText || '').includes('0 A / 0 R'))
                .map(element => ({
                    text: (element.innerText || '').trim(),
                    has_na: /(^|\\s)N\\/?A(\\s|$)/i.test(element.innerText || ''),
                    has_zero_percent: (element.innerText || '').includes('0.0%')
                }));
            return {
                fixture,
                heading: detail?.querySelector('h3')?.innerText || null,
                zero_zero_elements: candidates,
                has_primary_hp_explanation: (detail?.innerText || '').includes('primary HP 無有效 read coverage')
            };
        }""",
        fixture,
    )


def screenshot(page: Page, path: Path) -> str:
    """Capture the current viewport at a stable absolute path."""

    last_error: Optional[Exception] = None
    for attempt in range(3):
        try:
            page.screenshot(path=str(path), full_page=False)
            return str(path.resolve())
        except Exception as exc:  # Chromium can transiently reject captureScreenshot.
            last_error = exc
            if attempt < 2:
                page.wait_for_timeout(180 * (attempt + 1))
    assert last_error is not None
    raise last_error


def text_contrast_failures(page: Page) -> List[Dict[str, Any]]:
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


def scroll_element_to_viewport_top(page: Page, selector: str, offset: int = 8) -> None:
    """Place a visual-audit target at a predictable viewport position."""

    payload = {"selector": selector, "offset": offset}
    last_position: Dict[str, Any] = {}
    for attempt in range(4):
        last_position = page.evaluate(
            """payload => {
                const element = document.querySelector(payload.selector);
                if (!element) throw new Error('Screenshot target not found: ' + payload.selector);
                document.documentElement.style.scrollBehavior = 'auto';
                const navHeight = document.querySelector('.local-nav')?.getBoundingClientRect().height || 0;
                const top = element.getBoundingClientRect().top + window.scrollY - navHeight - payload.offset;
                window.scrollTo({top: Math.max(0, top), behavior: 'auto'});
                return {target_top: element.getBoundingClientRect().top, nav_height: navHeight};
            }""",
            payload,
        )
        page.wait_for_timeout(80 * (attempt + 1))
        last_position = page.evaluate(
            """payload => {
                const element = document.querySelector(payload.selector);
                const navHeight = document.querySelector('.local-nav')?.getBoundingClientRect().height || 0;
                return {target_top: element?.getBoundingClientRect().top ?? null, nav_height: navHeight};
            }""",
            payload,
        )
        if (
            last_position["target_top"] is not None
            and abs(last_position["target_top"] - last_position["nav_height"] - offset) <= 2
        ):
            return
    raise RuntimeError(
        f"Could not align screenshot target {selector}: {last_position} offset={offset}"
    )


def audit_sample_viewport(
    browser: Browser,
    workstation_dir: Path,
    output_dir: Path,
    sample: str,
    viewport_name: str,
    expectation: Dict[str, Any],
    timeout_ms: int,
) -> Dict[str, Any]:
    """Run the complete product contract for one sample and viewport."""

    started = time.monotonic()
    html_path = (workstation_dir / f"{sample}.html").resolve()
    run: Dict[str, Any] = {
        "kind": "sample",
        "sample": sample,
        "viewport": viewport_name,
        "viewport_size": VIEWPORTS[viewport_name],
        "input_path": str(html_path),
        "url": html_path.as_uri(),
        "sidecar_path": expectation["path"],
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
    page.on(
        "console",
        lambda message: run["console_errors"].append(
            {"text": message.text, "location": json_safe(message.location)}
        )
        if message.type == "error"
        else None,
    )
    page.on("pageerror", lambda error: run["page_errors"].append(str(error)))

    try:
        if not html_path.is_file():
            raise FileNotFoundError(f"Missing sample HTML: {html_path}")
        page.goto(html_path.as_uri(), wait_until="load", timeout=timeout_ms)
        page.wait_for_load_state("networkidle", timeout=timeout_ms)
        page.add_style_tag(content="html{scroll-behavior:auto!important}")
        page.wait_for_function(
            """expected => document.querySelectorAll('.ideogram-mark').length === expected""",
            arg=expectation["W_tree"],
            timeout=timeout_ms,
        )
        add_check(
            run,
            "document_load",
            True,
            expected="file:// document reaches load/networkidle and renders W_tree marks",
            actual={"ready_state": page.evaluate("document.readyState")},
        )

        meta = page.evaluate(
            """() => ({
                contract: document.querySelector('meta[name="intersubmod-ui-contract"]')?.content || null,
                sample: document.querySelector('meta[name="intersubmod-canonical-sample"]')?.content || null,
                sidecar_sha256: document.querySelector('meta[name="intersubmod-read-af-topology-sha256"]')?.content || null
            })"""
        )
        add_check(
            run,
            "ui_contract_and_sidecar_hash",
            meta == {
                "contract": EXPECTED_UI_CONTRACT,
                "sample": sample,
                "sidecar_sha256": expectation["sha256"],
            },
            expected={
                "contract": EXPECTED_UI_CONTRACT,
                "sample": sample,
                "sidecar_sha256": expectation["sha256"],
            },
            actual=meta,
        )

        edge_census_contract = expectation["full_edge_census_contract"]
        add_check(
            run,
            "sidecar_1_1_full_edge_census_contract_and_top_edge_inclusion",
            edge_census_contract["status"] == "pass"
            and edge_census_contract["required_units"]
            == edge_census_contract["census_units"]
            and edge_census_contract["top_representative_edges"]
            == edge_census_contract["top_representative_edges_in_census"]
            == edge_census_contract[
                "top_representative_edges_with_positive_top_count"
            ],
            expected={
                "schema": list(EXPECTED_SIDECAR_SCHEMA),
                "every_eligible_complete_T_gt_1_unit_has_census": True,
                "complete": True,
                "candidate_and_top_counts_bounded": True,
                "edge_rows_unique_sorted_and_totals_close": True,
                "every_top_representative_edge_in_census": True,
                "every_top_representative_edge_has_positive_top_count": True,
            },
            actual=edge_census_contract,
        )

        if sample == "HCC1395":
            fixed_data = expectation["hcc1395_edge_census_fixture"]
            fixed_census = fixed_data.get("census") or {}
            expected_fixed = HCC1395_EDGE_CENSUS_FIXTURE
            expected_edge_row = [
                *expected_fixed["edge"],
                expected_fixed["candidate_count"],
                expected_fixed["top_candidate_count"],
            ]
            add_check(
                run,
                "hcc1395_fixed_full_edge_census_fixture",
                fixed_data.get("region_found") is True
                and fixed_data.get("family_found") is True
                and fixed_census.get("complete") is True
                and fixed_census.get("candidate_total")
                == expected_fixed["candidate_total"]
                and fixed_census.get("top_candidate_total")
                == expected_fixed["top_candidate_total"]
                and fixed_census.get("node_total") == expected_fixed["node_total"]
                and fixed_census.get("edge_total") == expected_fixed["edge_total"]
                and fixed_data.get("edge_row") == expected_edge_row,
                expected={
                    "region": expected_fixed["region"],
                    "family": expected_fixed["family"],
                    "candidate_total": expected_fixed["candidate_total"],
                    "top_candidate_total": expected_fixed["top_candidate_total"],
                    "node_total": expected_fixed["node_total"],
                    "edge_total": expected_fixed["edge_total"],
                    "edge_row": expected_edge_row,
                },
                actual=fixed_data,
            )

        overview = page_overview_metrics(page)
        expected_denominators = {
            "topology-count": ("W_primary", expectation["W_primary"]),
            "candidate-count": ("W_primary", expectation["W_primary"]),
            "determinacy": ("W_primary", expectation["W_primary"]),
            "read-af-selection": ("W_primary", expectation["W_primary"]),
            "morphology": ("W_primary", expectation["W_primary"]),
            "hp-h3": ("W_tree", expectation["W_tree"]),
            "region-size": ("W_tree", expectation["W_tree"]),
        }
        panel_failures = []
        bar_failures = []
        for panel in overview["panels"]:
            denominator = expected_denominators.get(panel["id"])
            if not denominator or not (
                panel["denominator_key"] == denominator[0]
                and panel["total"] == denominator[1]
                and panel["sum"] == denominator[1]
            ):
                panel_failures.append(panel)
            for bin_record in panel["bins"]:
                expected_pct = 100.0 * bin_record["count"] / panel["total"] if panel["total"] else 0.0
                printed_pct = f"{expected_pct:.1f}%"
                if not (
                    bin_record["denominator"] == panel["total"]
                    and abs(bin_record["style_width_pct"] - expected_pct) <= 0.01
                    and abs(bin_record["measured_width_pct"] - expected_pct) <= 0.75
                    and printed_pct in bin_record["printed"]
                ):
                    bar_failures.append(
                        {
                            "panel": panel["id"],
                            "bin": bin_record,
                            "expected_pct": expected_pct,
                            "printed_pct": printed_pct,
                        }
                    )
        add_check(
            run,
            "seven_overview_panels_and_denominator_closure",
            overview["ids"] == list(EXPECTED_OVERVIEW_IDS)
            and not panel_failures
            and overview["W_tree"] == expectation["W_tree"]
            and overview["W_primary"] == expectation["W_primary"],
            expected={
                "ids": list(EXPECTED_OVERVIEW_IDS),
                "W_tree": expectation["W_tree"],
                "W_primary": expectation["W_primary"],
                "all_panel_sums_close": True,
            },
            actual={"overview": overview, "panel_failures": panel_failures},
        )
        add_check(
            run,
            "overview_bars_use_canonical_denominator",
            not bar_failures,
            expected={"style_and_measured_width": "100 * count / panel denominator", "printed_pct": "same denominator"},
            actual={"failure_count": len(bar_failures), "failures": bar_failures[:20]},
        )

        modes: Dict[str, Dict[str, Any]] = {}
        for mode in IDEOGRAM_MODES:
            modes[mode] = mode_metrics(page, mode)
        generic_modes_pass = all(
            state["active_mode"] == mode
            and state["active_mode_buttons"] == 1
            and state["marks"] == expectation["W_tree"]
            and state["histogram_sum"] == expectation["W_tree"]
            and state["unmapped"] == 0
            and state["legend_keys_match_histogram"]
            and state["legend_count_labels_exact"]
            for mode, state in modes.items()
        )
        add_check(
            run,
            "seven_ideogram_modes_category_conservation",
            generic_modes_pass and len(modes) == 7,
            expected={
                "modes": list(IDEOGRAM_MODES),
                "marks_and_histogram_sum": expectation["W_tree"],
                "unmapped": 0,
                "legend_matches_histogram": True,
            },
            actual=modes,
        )

        expected_selection = dict(expectation["selection_classes"])
        expected_selection["no_primary"] = expectation["no_primary"]
        expected_morphology = dict(expectation["morphology_classes"])
        expected_morphology["not_applicable"] = expectation["no_primary"]
        add_check(
            run,
            "read_af_and_morphology_histograms_match_sidecar",
            modes["read-af-selection"]["histogram"] == expected_selection
            and modes["morphology"]["histogram"] == expected_morphology,
            expected={
                "read-af-selection": expected_selection,
                "morphology": expected_morphology,
            },
            actual={
                "read-af-selection": modes["read-af-selection"]["histogram"],
                "morphology": modes["morphology"]["histogram"],
            },
        )

        morphology_signatures = modes["morphology"]["signatures"]
        single_signature = morphology_signatures.get("single_no_within_hp_relation", {})
        na_signature = morphology_signatures.get("not_applicable", {})
        single_noncolor = (
            single_signature.get("stroke_width"),
            single_signature.get("stroke_dasharray"),
            single_signature.get("swatch_background_image"),
        )
        na_noncolor = (
            na_signature.get("stroke_width"),
            na_signature.get("stroke_dasharray"),
            na_signature.get("swatch_background_image"),
        )
        add_check(
            run,
            "morphology_single_and_na_have_distinct_noncolor_encoding",
            bool(single_signature)
            and bool(na_signature)
            and single_noncolor != na_noncolor
            and na_signature.get("stroke_dasharray") not in {"", "none"}
            and na_signature.get("swatch_background_image") not in {"", "none"},
            expected={"single": "solid", "not_applicable": "dashed mark + hatched legend"},
            actual={"single": single_signature, "not_applicable": na_signature},
        )

        all_mode_toggle_receipts: Dict[str, Any] = {}
        for mode in IDEOGRAM_MODES:
            state = mode_metrics(page, mode)
            first_key = state["legend"][0]["key"]
            expected_count = state["histogram"][first_key]
            button = page.locator(f'[data-legend-key="{first_key}"]')
            button.click()
            selected = selection_state(page)
            button.click()
            empty = selection_state(page)
            all_mode_toggle_receipts[mode] = {
                "key": first_key,
                "expected_count": expected_count,
                "selected": selected,
                "empty": empty,
            }
        add_check(
            run,
            "all_ideogram_modes_single_toggle_then_empty_means_all",
            all(
                receipt["selected"]["selected"] == [receipt["key"]]
                and receipt["selected"]["visible"] == receipt["expected_count"]
                and receipt["empty"]["selected"] == []
                and receipt["empty"]["visible"] == expectation["W_tree"]
                and receipt["empty"]["dimmed"] == 0
                and receipt["empty"]["all_pressed"] == "true"
                for receipt in all_mode_toggle_receipts.values()
            ),
            expected={"modes": list(IDEOGRAM_MODES), "sequence": "one category → second click → empty/all"},
            actual=all_mode_toggle_receipts,
        )

        page.locator('[data-ideogram-mode="read-af-selection"]').click()
        count_a = expected_selection["structural_exact_unique"]
        count_b = expected_selection["read_af_unique_first"]
        button_a = page.locator('[data-legend-key="structural_exact_unique"]')
        button_b = page.locator('[data-legend-key="read_af_unique_first"]')
        button_a.click()
        state_a = selection_state(page)
        button_b.click()
        state_ab = selection_state(page)
        button_a.click()
        state_b = selection_state(page)
        button_b.click()
        state_empty = selection_state(page)
        multiselect_pass = (
            state_a["selected"] == ["structural_exact_unique"]
            and state_a["visible"] == count_a
            and state_ab["selected"]
            == ["read_af_unique_first", "structural_exact_unique"]
            and state_ab["visible"] == count_a + count_b
            and state_b["selected"] == ["read_af_unique_first"]
            and state_b["visible"] == count_b
            and state_empty["selected"] == []
            and state_empty["visible"] == expectation["W_tree"]
            and state_empty["dimmed"] == 0
            and state_empty["all_pressed"] == "true"
        )
        add_check(
            run,
            "read_af_multiselect_A_union_B_toggle_and_empty_all",
            multiselect_pass,
            expected={
                "sequence": "A -> A+B -> B -> empty",
                "visible": [count_a, count_a + count_b, count_b, expectation["W_tree"]],
                "empty_means_all": True,
            },
            actual={"A": state_a, "A+B": state_ab, "B": state_b, "empty": state_empty},
        )

        button_a.click()
        page.locator('[data-ideogram-mode="read-af-selection"]').click()
        same_mode = selection_state(page)
        page.locator('[data-ideogram-mode="morphology"]').click()
        different_mode = selection_state(page)
        add_check(
            run,
            "selection_lifecycle_same_mode_retains_different_mode_clears",
            same_mode["selected"] == ["structural_exact_unique"]
            and same_mode["visible"] == count_a
            and different_mode["selected"] == []
            and different_mode["visible"] == expectation["W_tree"]
            and different_mode["all_pressed"] == "true",
            expected={"same_mode_retains": True, "different_mode_clears": True},
            actual={"same_mode": same_mode, "different_mode": different_mode},
        )

        overview_button = page.locator(
            '[data-overview-panel="read-af-selection"] '
            '[data-bin-key="read_af_unique_first"]'
        )
        overview_button.click()
        overview_sync = page.evaluate(
            """() => ({
                mode: document.getElementById('genome-ideogram').dataset.mode,
                mode_pressed: document.querySelector('[data-ideogram-mode="read-af-selection"]')
                    ?.getAttribute('aria-pressed') || null,
                legend_pressed: document.querySelector('[data-legend-key="read_af_unique_first"]')
                    ?.getAttribute('aria-pressed') || null,
                overview_pressed: document.querySelector(
                    '[data-overview-panel="read-af-selection"] [data-bin-key="read_af_unique_first"]'
                )?.getAttribute('aria-pressed') || null,
                visible: document.querySelectorAll('.ideogram-mark:not(.dimmed)').length
            })"""
        )
        add_check(
            run,
            "overview_read_af_bin_synchronizes_ideogram",
            overview_sync
            == {
                "mode": "read-af-selection",
                "mode_pressed": "true",
                "legend_pressed": "true",
                "overview_pressed": "true",
                "visible": count_b,
            },
            expected={
                "mode": "read-af-selection",
                "mode_pressed": "true",
                "legend_pressed": "true",
                "overview_pressed": "true",
                "visible": count_b,
            },
            actual=overview_sync,
        )
        page.locator("[data-legend-clear]").click()

        information_architecture = page.evaluate(
            """() => {
                const ids = ['sample-summary','genome-overview','region-browser','dimension-guide','sample-overview','method-evidence'];
                const nodes = ids.map(id => document.getElementById(id));
                const jsonLinks = [...document.querySelectorAll('a[href]')].filter(link => (link.getAttribute('href') || '').toLowerCase().includes('.json'));
                const rail = document.getElementById('ideogram-mobile-rail');
                return {
                    ids_present: nodes.every(Boolean),
                    dom_order: nodes.every((node,index) => index === nodes.length - 1 || Boolean(node.compareDocumentPosition(nodes[index + 1]) & Node.DOCUMENT_POSITION_FOLLOWING)),
                    nav_links: document.querySelectorAll('.local-nav a').length,
                    nav_height: document.querySelector('.local-nav')?.getBoundingClientRect().height || 0,
                    current_nav_links: document.querySelectorAll('.local-nav a[aria-current="location"]').length,
                    chrom_structural_legend: (document.getElementById('chrom-structure-key')?.innerText || '').includes('structural determinacy'),
                    rail_labels: rail?.querySelectorAll('span').length || 0,
                    rail_display: rail ? getComputedStyle(rail).display : null,
                    evidence_closed: document.getElementById('method-evidence')?.open === false,
                    json_links: jsonLinks.length,
                    json_links_inside_evidence: jsonLinks.filter(link => link.closest('#method-evidence')).length
                };
            }"""
        )
        expected_rail_display = "none" if viewport_name == "desktop" else "block"
        add_check(
            run,
            "answer_first_local_nav_structural_legend_and_mobile_rail",
            information_architecture["ids_present"]
            and information_architecture["dom_order"]
            and information_architecture["nav_links"] == 6
            and information_architecture["nav_height"] >= 44
            and information_architecture["current_nav_links"] == 1
            and information_architecture["chrom_structural_legend"]
            and information_architecture["rail_labels"] == 22
            and information_architecture["rail_display"] == expected_rail_display
            and information_architecture["evidence_closed"]
            and information_architecture["json_links"] > 0
            and information_architecture["json_links_inside_evidence"] == information_architecture["json_links"],
            expected={"order": ["sample-summary", "genome-overview", "region-browser", "dimension-guide", "sample-overview", "method-evidence"], "nav_links": 6, "rail_display": expected_rail_display, "json_hidden": True},
            actual=information_architecture,
        )

        base_overflow = body_overflow_metrics(page)
        add_check(
            run,
            "body_horizontal_overflow_at_genome",
            float(base_overflow["overflow_px"]) <= 1,
            expected={"overflow_px_max": 1},
            actual=base_overflow,
        )

        if viewport_name != "desktop":
            touch = touch_target_metrics(page)
            add_check(
                run,
                "mobile_genome_controls_touch_height",
                bool(touch["controls"]) and not touch["below_44px"],
                expected={"minimum_height_px": 44, "below_44px": []},
                actual=touch,
            )

        genome_target = "#genome-overview" if viewport_name == "desktop" else ".ideogram-head"
        scroll_element_to_viewport_top(page, genome_target)
        page.wait_for_timeout(180)
        genome_path = output_dir / f"{sample}__{viewport_name}__genome.png"
        run["screenshots"]["genome"] = screenshot(page, genome_path)
        genome_alignment = page.evaluate(
            """selector => {
                const target = document.querySelector(selector);
                const nav = document.querySelector('.local-nav');
                const heading = document.querySelector('.ideogram-head h3');
                return {
                    target_top: target?.getBoundingClientRect().top || null,
                    nav_bottom: nav?.getBoundingClientRect().bottom || null,
                    heading_top: heading?.getBoundingClientRect().top || null,
                    current_href: document.querySelector('.local-nav a[aria-current="location"]')?.getAttribute('href') || null
                };
            }""",
            genome_target,
        )
        add_check(
            run,
            "genome_screenshot_heading_unobscured_and_nav_current",
            genome_alignment["target_top"] is not None
            and genome_alignment["nav_bottom"] is not None
            and genome_alignment["heading_top"] is not None
            and genome_alignment["target_top"] >= genome_alignment["nav_bottom"] + 6
            and genome_alignment["heading_top"] >= genome_alignment["nav_bottom"] + 6
            and genome_alignment["current_href"] == "#genome-overview",
            expected={"target_gap_min_px": 6, "heading_gap_min_px": 6, "current_href": "#genome-overview"},
            actual=genome_alignment,
        )

        ranking_fixture = expectation["ranking_fixture"]
        select_region(page, ranking_fixture["region"], timeout_ms)
        ranking = ranking_detail_metrics(page, ranking_fixture["region"])
        add_check(
            run,
            "real_read_af_detail_has_card_table_tree_and_exact_fraction",
            ranking["heading"] is not None
            and ranking_fixture["region"] in ranking["heading"]
            and ranking["ranking_cards"] >= 1
            and ranking["ranking_tables"] >= 1
            and ranking["ranked_trees"] >= 1
            and ranking["ranked_networks"] >= 1
            and ranking["exact_fraction_text"]
            and ranking["probability_boundary_text"],
            expected={
                "fixture": ranking_fixture,
                "ranking_cards_min": 1,
                "ranking_tables_min": 1,
                "ranked_trees_min": 1,
                "Exact Fraction": True,
                "probability_claim_boundary": True,
            },
            actual=ranking,
        )
        add_check(
            run,
            "read_af_top_trees_have_initially_visible_nodes",
            ranking["ranked_trees"] >= 1
            and ranking["ranked_trees_with_visible_nodes"] == ranking["ranked_trees"]
            and all(
                view["scroller"] is not None
                and view["svg"] is not None
                and view["node_count"] >= 1
                and view["visible_nodes"] >= 1
                for view in ranking["ranked_tree_views"]
            ),
            expected={
                "scope": "every read-AF top-tree in every sample ranking fixture",
                "measurement": "circle/rect bounding-box intersection with local network scroller viewport",
                "node_count_minimum_each": 1,
                "visible_nodes_minimum_each": 1,
            },
            actual={
                "ranked_trees": ranking["ranked_trees"],
                "ranked_trees_with_visible_nodes": ranking["ranked_trees_with_visible_nodes"],
                "ranked_tree_views": ranking["ranked_tree_views"],
            },
        )

        if sample == "HCC1395":
            fixed = HCC1395_EDGE_CENSUS_FIXTURE
            select_region(page, fixed["region"], timeout_ms)
            census_ui = hcc1395_edge_census_cross_panel_metrics(page)
            expected_target = {
                "key": ">".join(fixed["edge"]),
                "candidate_count": fixed["candidate_count"],
                "candidate_total": fixed["candidate_total"],
                "top_count": fixed["top_candidate_count"],
                "top_total": fixed["top_candidate_total"],
                "status": "variable",
                "selected": True,
                "visible": True,
            }
            add_check(
                run,
                "hcc1395_full_edge_union_selected_overlay_and_stored_preview_regression",
                census_ui["detail_heading"] is not None
                and fixed["region"] in census_ui["detail_heading"]
                and census_ui["matching_full_network_found"]
                and census_ui["full_edge_count"] == fixed["edge_total"]
                and census_ui["full_unique_edge_count"] == fixed["edge_total"]
                and census_ui["all_full_edges_have_expected_totals"]
                and census_ui["all_full_edge_statuses_match_counts"]
                and census_ui["all_selected_classes_match_positive_top_count"]
                and census_ui["target_edge"] == expected_target
                and census_ui["selected_overlay"]
                == {"key": expected_target["key"], "visible": True}
                and census_ui["edge_census_table_count"] >= 1
                and census_ui["edge_census_table_visible"]
                and census_ui["stored_preview_scope"] is not None
                and census_ui["stored_preview_scope"]["stored_count"]
                == fixed["stored_candidate_total"]
                and census_ui["stored_preview_scope"]["candidate_total"]
                == fixed["candidate_total"]
                and census_ui["stored_preview_scope"]["visible"]
                and str(fixed["stored_candidate_total"])
                in census_ui["stored_preview_scope"]["text"]
                and str(fixed["candidate_total"])
                in census_ui["stored_preview_scope"]["text"]
                and census_ui["read_af_top_edge"]
                == {"key": expected_target["key"], "visible": True}
                and census_ui["full_union_boundary_text"]
                and census_ui["not_single_tree_boundary_text"]
                and census_ui["non_probability_boundary_text"],
                expected={
                    "fixture": fixed,
                    "full_union_selector": (
                        '.full-edge-network[data-edge-scope="complete"]'
                    ),
                    "full_edge_count": fixed["edge_total"],
                    "target_edge": expected_target,
                    "selected_overlay": {
                        "selector": ".network-edge-selected-overlay[data-edge-key]",
                        "visible": True,
                    },
                    "edge_census_table_visible": True,
                    "stored_preview": {
                        "stored_count": fixed["stored_candidate_total"],
                        "candidate_total": fixed["candidate_total"],
                        "counts_explicit_in_text": True,
                    },
                    "read_af_top_contains_same_edge": True,
                    "claim_boundary": (
                        "union is not one tree; n/N and top m/M are not probability"
                    ),
                },
                actual=census_ui,
            )

        detail_overflow = body_overflow_metrics(page)
        add_check(
            run,
            "body_horizontal_overflow_with_ranking_detail",
            float(detail_overflow["overflow_px"]) <= 1,
            expected={"overflow_px_max": 1},
            actual=detail_overflow,
        )

        contrast_failures = text_contrast_failures(page)
        add_check(
            run,
            "sample_visible_text_wcag_aa_contrast",
            not contrast_failures,
            expected={"failures": 0, "normal_text_ratio": 4.5, "large_text_ratio": 3.0},
            actual={"failures": len(contrast_failures), "examples": contrast_failures},
        )

        if sample == "HCC1395":
            # Put the evidence under audit—not only the region heading—inside
            # the visual receipt.  This makes the ranking table and first-ranked
            # network directly inspectable on both desktop and mobile.
            scroll_element_to_viewport_top(page, "#detail .read-af-card")
            page.wait_for_timeout(180)
            detail_path = output_dir / f"{sample}__{viewport_name}__detail.png"
            run["screenshots"]["detail"] = screenshot(page, detail_path)

        if expectation["zero_coverage_fixture"] is None:
            expectation["zero_coverage_fixture"] = find_zero_coverage_region(page)
        zero_fixture = expectation["zero_coverage_fixture"]
        if zero_fixture:
            select_region(page, zero_fixture["region"], timeout_ms)
            zero_metrics = zero_coverage_detail_metrics(page, zero_fixture)
            zero_passed = (
                zero_metrics["heading"] is not None
                and zero_fixture["region"] in zero_metrics["heading"]
                and bool(zero_metrics["zero_zero_elements"])
                and all(item["has_na"] for item in zero_metrics["zero_zero_elements"])
                and not any(item["has_zero_percent"] for item in zero_metrics["zero_zero_elements"])
            )
        else:
            zero_metrics = {"fixture": None, "reason": "no primary HP 0/0 fixture found"}
            zero_passed = False
        add_check(
            run,
            "primary_HP_zero_zero_renders_NA_not_zero_percent",
            zero_passed,
            expected={
                "fixture_found": True,
                "0_A_0_R_elements_min": 1,
                "all_show_NA": True,
                "none_show_0.0_percent": True,
            },
            actual=zero_metrics,
        )

        page.wait_for_timeout(120)
        add_check(
            run,
            "console_and_page_errors",
            not run["console_errors"] and not run["page_errors"],
            expected={"console_errors": 0, "page_errors": 0},
            actual={
                "console_errors": run["console_errors"],
                "page_errors": run["page_errors"],
            },
        )
    except Exception as exc:
        add_check(
            run,
            "audit_completed_without_exception",
            False,
            expected="all page assertions complete",
            actual={"error": str(exc)},
        )
    finally:
        run["duration_seconds"] = round(time.monotonic() - started, 3)
        context.close()
    return run


def summarize(report: Dict[str, Any]) -> None:
    """Attach aggregate result counts and failure coordinates."""

    index_runs = report.get("index_runs", [])
    sample_runs = report.get("runs", [])
    all_runs = list(index_runs) + list(sample_runs)
    checks: List[Dict[str, Any]] = []
    for run in all_runs:
        checks.extend(run.get("checks", []))
    failed = [check for check in checks if check.get("status") != "pass"]
    failing_runs = []
    for run in all_runs:
        names = [
            check["name"]
            for check in run.get("checks", [])
            if check.get("status") != "pass"
        ]
        if names:
            failing_runs.append(
                {
                    "kind": run.get("kind", "sample"),
                    "sample": run.get("sample"),
                    "viewport": run["viewport"],
                    "failed_checks": names,
                }
            )
    screenshots = [
        path
        for run in all_runs
        for path in run.get("screenshots", {}).values()
        if path
    ]
    report["summary"] = {
        "status": "pass" if not failed and not report.get("fatal_errors") else "fail",
        "samples_tested": len({run["sample"] for run in sample_runs}),
        "documents_tested": 1 + len({run["sample"] for run in sample_runs}),
        "sample_runs": len(sample_runs),
        "index_runs": len(index_runs),
        "page_runs": len(all_runs),
        "page_runs_total": len(all_runs),
        "checks_total": len(checks),
        "checks_passed": len(checks) - len(failed),
        "checks_failed": len(failed),
        "page_runs_failed": len(failing_runs),
        "failing_runs": failing_runs,
        "screenshots_written": len(screenshots),
        "screenshot_paths": screenshots,
    }


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line interface."""

    parser = argparse.ArgumentParser(
        description=(
            "Audit the canonical cohort index and all seven layered workstation "
            "sample pages in local headless Chromium at desktop/mobile viewports."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--workstation-dir",
        type=Path,
        default=DEFAULT_WORKSTATION_DIR,
        help="Directory containing index.html and the seven generated sample HTML files.",
    )
    parser.add_argument(
        "--sidecar-dir",
        type=Path,
        default=DEFAULT_SIDECAR_DIR,
        help="Directory containing current-v5 read-AF topology sidecars and index.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "Directory for 14 genome screenshots, two HCC1395 detail screenshots, "
            "two index interpretation screenshots, and receipt."
        ),
    )
    parser.add_argument(
        "--timeout-ms",
        type=int,
        default=180_000,
        help="Timeout for each navigation, selector, and interaction.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the audit and write the complete JSON receipt."""

    parser = build_parser()
    args = parser.parse_args(argv)
    started_at = utc_now()
    workstation_dir = args.workstation_dir.expanduser().resolve()
    sidecar_dir = args.sidecar_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    receipt_path = output_dir / "validation_receipt.json"
    report: Dict[str, Any] = {
        "schema_name": "intersubmod.layered_workstation_playwright_audit",
        "schema_version": "1.1.0",
        "tool": str(SCRIPT_PATH),
        "started_at": started_at,
        "finished_at": None,
        "scope": {
            "task_type": "B_comprehensive_validation",
            "partial": False,
            "samples": list(SAMPLES),
            "viewports": VIEWPORTS,
            "sample_page_runs_expected": len(SAMPLES) * len(VIEWPORTS),
            "index_page_runs_expected": len(VIEWPORTS),
            "page_runs_expected": (len(SAMPLES) + 1) * len(VIEWPORTS),
            "screenshots_expected": len(SAMPLES) * len(VIEWPORTS) + 2 * len(VIEWPORTS),
        },
        "inputs": {
            "workstation_index": str((workstation_dir / "index.html").resolve()),
            "workstation_dir": str(workstation_dir),
            "sidecar_dir": str(sidecar_dir),
        },
        "output": {"receipt": str(receipt_path)},
        "browser": None,
        "sidecars": {},
        "sidecar_index": {},
        "index_runs": [],
        "runs": [],
        "fatal_errors": [],
    }

    if args.timeout_ms <= 0:
        report["fatal_errors"].append("--timeout-ms must be greater than zero")
    if PLAYWRIGHT_IMPORT_ERROR:
        report["fatal_errors"].append(
            f"Python Playwright import failed: {PLAYWRIGHT_IMPORT_ERROR}"
        )

    expectations: Dict[str, Dict[str, Any]] = {}
    sidecar_index_contract: Dict[str, Any] = {}
    if not report["fatal_errors"]:
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            expectations, sidecar_index_contract = load_sidecar_expectations(sidecar_dir)
            report["sidecar_index"] = json_safe(sidecar_index_contract)
            report["sidecars"] = {
                sample: {
                    key: json_safe(value)
                    for key, value in expectation.items()
                    if key != "zero_coverage_fixture"
                }
                for sample, expectation in expectations.items()
            }
        except Exception as exc:
            report["fatal_errors"].append(str(exc))

    playwright = None
    browser: Optional[Browser] = None
    if not report["fatal_errors"]:
        try:
            playwright = sync_playwright().start()  # type: ignore[union-attr]
            browser = playwright.chromium.launch(
                headless=True,
                args=["--allow-file-access-from-files"],
            )
            report["browser"] = {
                "engine": "chromium",
                "version": browser.version,
                "headless": True,
                "file_access_from_files": True,
            }
            for viewport_name in VIEWPORTS:
                print(f"[RUN] index {viewport_name}", flush=True)
                index_run = index_audit(
                    browser,
                    workstation_dir,
                    output_dir,
                    viewport_name,
                    sidecar_index_contract,
                    args.timeout_ms,
                )
                report["index_runs"].append(index_run)
                index_failures = [
                    check["name"]
                    for check in index_run["checks"]
                    if check["status"] != "pass"
                ]
                index_status = "PASS" if not index_failures else "FAIL"
                print(
                    f"[{index_status}] index {viewport_name} "
                    f"{index_run['duration_seconds']}s"
                    + (
                        f" failures={','.join(index_failures)}"
                        if index_failures
                        else ""
                    ),
                    flush=True,
                )
            for sample in SAMPLES:
                for viewport_name in VIEWPORTS:
                    print(f"[RUN] {sample} {viewport_name}", flush=True)
                    run = audit_sample_viewport(
                        browser,
                        workstation_dir,
                        output_dir,
                        sample,
                        viewport_name,
                        expectations[sample],
                        args.timeout_ms,
                    )
                    report["runs"].append(run)
                    failures = [
                        check["name"]
                        for check in run["checks"]
                        if check["status"] != "pass"
                    ]
                    status = "PASS" if not failures else "FAIL"
                    print(
                        f"[{status}] {sample} {viewport_name} "
                        f"{run['duration_seconds']}s"
                        + (f" failures={','.join(failures)}" if failures else ""),
                        flush=True,
                    )
        except Exception as exc:
            report["fatal_errors"].append(str(exc))
        finally:
            if browser is not None:
                try:
                    browser.close()
                except Exception as exc:
                    report["fatal_errors"].append(f"Browser close failed: {exc}")
            if playwright is not None:
                try:
                    playwright.stop()
                except Exception as exc:
                    report["fatal_errors"].append(f"Playwright stop failed: {exc}")

    report["finished_at"] = utc_now()
    summarize(report)
    if report["fatal_errors"]:
        report["summary"]["status"] = "error"
        exit_code = 2
    elif report["summary"]["checks_failed"]:
        exit_code = 1
    else:
        exit_code = 0
    report["exit_code"] = exit_code
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        receipt_path.write_text(
            json.dumps(report, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
    except Exception as exc:
        print(f"[ERROR] Cannot write receipt: {exc}", file=sys.stderr)
        return 2
    print(
        f"[SUMMARY] status={report['summary']['status']} "
        f"runs={report['summary']['page_runs']} "
        f"checks={report['summary']['checks_passed']}/{report['summary']['checks_total']} "
        f"screenshots={report['summary']['screenshots_written']}",
        flush=True,
    )
    print(f"[OUTPUT] {receipt_path}", flush=True)
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
