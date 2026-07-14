#!/usr/bin/env python3
"""Build the canonical seven-dataset layered reconstruction workstation.

The repository machine summary is the only authority for sample paths and all
landing-page measurements.  The command fails closed before writing index.html
when that authority, any bound artifact, or any sample-page provenance marker
is stale or inconsistent.

Default mode delegates the seven large sample pages to
``build_layered_workstation_v5.py`` and then writes the cohort index.  ``--index-only``
never rebuilds sample pages; it succeeds only when all seven existing pages are
already bound to the current summary and region-view hashes.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import subprocess
import sys
from datetime import datetime
from html.parser import HTMLParser
from pathlib import Path
from string import Template
from typing import Any


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]
BUILD = HERE / "build_layered_workstation_v5.py"
OUTDIR = (HERE / ".." / ".." / "layered_workstation").resolve()
CURRENT_SUMMARY = (
    REPO_ROOT
    / "research"
    / "20260710_layered_reconstruction_v2"
    / "current_layered_topology_v3_raw_all_v1.json"
)
PYTHON = sys.executable or "python3"

EXPECTED_SAMPLES = {
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
}
EXPECTED_SCOPE = "7 datasets / 6 biological samples / chr1-22"
EXPECTED_BACKBONE = "longphase_s_recalibrated_FILTER_PASS"
EXPECTED_TREE_SOURCE = "longphase_s_recalibrated_filter_pass"
EXPECTED_REGION_SCOPE = "chr1-22 primary; chrX/chrY out-of-scope census only"
PAGE_META_NAMES = {
    "summary_sha256": "intersubmod-current-summary-sha256",
    "region_sha256": "intersubmod-region-view-sha256",
    "sample": "intersubmod-canonical-sample",
}
TOPOLOGY_CLASSES = (
    ("exact_and_topology_unique", "C=1 / Topo=1", "exact 與拓撲皆唯一", "unique"),
    ("topology_unique_exact_multiple", "C>1 / Topo=1", "exact 多解、拓撲唯一", "shape"),
    ("topology_multiple_exact_multiple", "C>1 / Topo>1", "exact 與拓撲皆多解", "multiple"),
    ("incomplete", "Incomplete", "候選集合未完整", "incomplete"),
)


class AuthorityError(RuntimeError):
    """Raised when canonical provenance cannot be proven exactly."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuthorityError(message)


def escaped(value: Any) -> str:
    return html.escape(str(value), quote=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AuthorityError(f"cannot read canonical JSON {path}: {exc}") from exc
    require(isinstance(payload, dict), f"canonical JSON is not an object: {path}")
    return payload


def parse_timestamp(value: str, label: str) -> datetime:
    try:
        timestamp = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except (TypeError, ValueError) as exc:
        raise AuthorityError(f"invalid {label} timestamp: {value!r}") from exc
    require(timestamp.tzinfo is not None, f"{label} timestamp has no timezone")
    return timestamp


def resolve_authority_path(value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (REPO_ROOT / path).resolve()


def require_under(path: Path, parent: Path, label: str) -> None:
    try:
        path.relative_to(parent)
    except ValueError as exc:
        raise AuthorityError(f"{label} escapes canonical run root: {path}") from exc


def verify_bound_file(path: Path, expected_sha256: str, label: str) -> None:
    require(path.is_file(), f"missing {label}: {path}")
    actual = sha256_file(path)
    require(actual == expected_sha256, f"{label} hash drift: expected={expected_sha256} actual={actual}")


def validate_aggregate(aggregate: dict[str, Any], samples: list[dict[str, Any]]) -> None:
    require(aggregate.get("dataset_count") == 7, "canonical aggregate dataset_count must be 7")
    require(aggregate.get("all_invariants_pass") is True, "canonical aggregate invariants are not all PASS")

    additive = (
        "tree_input_records",
        "autosomal_biallelic_sSNV",
        "retained_sSNV",
        "W_tree",
        "W_primary",
        "no_primary_lineage",
        "primary_units",
        "primary_units_full_and_partial",
        "primary_units_partial_only",
        "primary_units_full_only",
        "complete_regions",
        "incomplete_regions",
        "read_tag_exposures",
        "read_tag_exact_matches",
        "mixed_PS_regions",
    )
    for field in additive:
        observed = sum(int(sample[field]) for sample in samples)
        require(observed == int(aggregate[field]), f"canonical aggregate mismatch for {field}")

    topology = aggregate["topology_classes"]
    for key, _, _, _ in TOPOLOGY_CLASSES:
        observed = sum(int(sample["topology_classes"][key]) for sample in samples)
        require(observed == int(topology[key]), f"canonical topology aggregate mismatch for {key}")
    require(int(topology.get("impossible_exact_unique_topology_multiple", -1)) == 0,
            "impossible C=1 / Topo>1 state is non-zero")
    require(int(aggregate["W_tree"]) == int(aggregate["W_primary"]) + int(aggregate["no_primary_lineage"]),
            "W_tree does not conserve W_primary + no-primary")
    require(int(aggregate["W_primary"]) == int(aggregate["complete_regions"]) + int(aggregate["incomplete_regions"]),
            "W_primary does not conserve complete + incomplete")
    require(int(aggregate["W_primary"]) == sum(int(topology[key]) for key, _, _, _ in TOPOLOGY_CLASSES),
            "topology classes do not conserve W_primary")
    require(int(aggregate["read_tag_exposures"]) == int(aggregate["read_tag_exact_matches"]),
            "read-tag exposures do not equal exact matches")


def load_authority() -> dict[str, Any]:
    """Load and exhaustively bind the current canonical machine summary."""
    require(CURRENT_SUMMARY.is_file(), f"missing current summary: {CURRENT_SUMMARY}")
    summary_sha256 = sha256_file(CURRENT_SUMMARY)
    summary = load_json(CURRENT_SUMMARY)
    require(summary.get("schema_name") == "intersubmod.current_layered_topology_summary",
            "unexpected current-summary schema_name")
    require(summary.get("schema_version") == "1.0.0", "unexpected current-summary schema_version")
    require(summary.get("all_pass") is True, "current summary all_pass is not true")
    require(summary.get("task_type") == "B_comprehensive_validation", "current summary is not Task B")
    require(summary.get("scope") == EXPECTED_SCOPE, f"current summary scope is not {EXPECTED_SCOPE}")
    require(summary.get("claim_scope") == "regional mutation-state candidate trees; not confirmed cell clones",
            "current summary claim boundary drifted")
    definitions = summary.get("definitions", {})
    require(definitions.get("PS") == "phase-block QC context only; not a topology edge or lineage label",
            "current summary PS contract drifted")

    canonical = summary.get("canonical")
    require(isinstance(canonical, dict), "current summary has no canonical object")
    require(canonical.get("label") == EXPECTED_BACKBONE, "canonical backbone is not LongPhase-S recalibrated PASS")
    run_root = resolve_authority_path(canonical["run_root"])
    require(run_root.is_dir(), f"missing canonical run root: {run_root}")

    success_path = run_root / "_SUCCESS"
    verify_bound_file(success_path, canonical["success_sha256"], "canonical _SUCCESS")
    success = load_json(success_path)
    require(success.get("run_id") == run_root.name, "_SUCCESS run_id/root mismatch")
    require(success.get("extra", {}).get("dataset_count") == 7, "_SUCCESS dataset_count is not 7")
    require(success.get("extra", {}).get("biological_sample_count") == 6,
            "_SUCCESS biological_sample_count is not 6")
    require(success.get("extra", {}).get("scope") == "chr1-22", "_SUCCESS scope is not chr1-22")
    require(success.get("extra", {}).get("mode") == "full", "_SUCCESS mode is not full")

    verification_path = resolve_authority_path(canonical["verification_path"])
    require_under(verification_path, run_root, "verification summary")
    verify_bound_file(verification_path, canonical["verification_sha256"], "verification summary")
    verification = load_json(verification_path)
    require(verification.get("all_pass") is True, "verification all_pass is not true")
    require(verification.get("dataset_count") == 7, "verification dataset_count is not 7")
    require(verification.get("n_pass") == 7 and verification.get("n_fail") == 0,
            "verification is not 7/7 PASS")

    records = canonical.get("samples", [])
    require(isinstance(records, list) and len(records) == 7, "canonical sample list is not exactly seven records")
    names = [record.get("sample") for record in records]
    require(len(names) == len(set(names)), "canonical sample list contains duplicates")
    require(set(names) == EXPECTED_SAMPLES, f"canonical sample set drifted: {names}")
    verification_names = {record.get("sample") for record in verification.get("samples", [])}
    require(verification_names == EXPECTED_SAMPLES, "verification sample set drifted")

    summary_generated = parse_timestamp(summary["generated_at"], "current summary generated_at")
    max_source_mtime = 0.0
    enriched_records: list[dict[str, Any]] = []
    for record in records:
        sample = record["sample"]
        require(record.get("pass") is True, f"canonical sample is not PASS: {sample}")
        require(record.get("all_invariants_pass") is True, f"sample invariants are not PASS: {sample}")
        require(record.get("tree_backbone_source") == EXPECTED_TREE_SOURCE,
                f"tree source drifted for {sample}")
        require(int(record["W_tree"]) == int(record["W_primary"]) + int(record["no_primary_lineage"]),
                f"W conservation failed for {sample}")
        require(int(record["W_primary"]) == int(record["complete_regions"]) + int(record["incomplete_regions"]),
                f"complete/incomplete conservation failed for {sample}")
        sample_topology = record["topology_classes"]
        require(int(record["W_primary"]) == sum(int(sample_topology[key]) for key, _, _, _ in TOPOLOGY_CLASSES),
                f"topology conservation failed for {sample}")
        require(int(sample_topology.get("impossible_exact_unique_topology_multiple", -1)) == 0,
                f"impossible topology state present for {sample}")

        region_path = resolve_authority_path(record["paths"]["layered_region_view"])
        layered_path = resolve_authority_path(record["paths"]["layered_reconstruction"])
        require_under(region_path, run_root, f"{sample} region view")
        require_under(layered_path, run_root, f"{sample} layered reconstruction")
        verify_bound_file(region_path, record["sha256"]["layered_region_view"], f"{sample} region view")
        verify_bound_file(layered_path, record["sha256"]["layered_reconstruction"],
                          f"{sample} layered reconstruction")
        max_source_mtime = max(max_source_mtime, region_path.stat().st_mtime, layered_path.stat().st_mtime)

        region_payload = load_json(region_path)
        require(region_payload.get("sample") == sample, f"region-view sample mismatch: {sample}")
        require(region_payload.get("schema_version") == "2.0", f"region-view schema drifted: {sample}")
        census = region_payload.get("census", {})
        require(census.get("U1_backbone_source") == EXPECTED_TREE_SOURCE,
                f"region-view backbone drifted: {sample}")
        require(census.get("analysis_scope") == EXPECTED_REGION_SCOPE,
                f"region-view scope drifted: {sample}")
        require(int(census.get("n_regions", -1)) == int(record["W_tree"]),
                f"region-view W_tree drifted: {sample}")
        require(int(census.get("U1_sSNV_somatic_total", -1)) == int(record["tree_input_records"]),
                f"region-view tree-input count drifted: {sample}")
        l3 = region_payload.get("L3_methyl", {})
        require(l3.get("status") == "not_evaluated" and l3.get("role") == "bounded_auxiliary",
                f"L3 contract drifted: {sample}")
        ps_contract = region_payload.get("analysis_contract", {}).get("PS", "")
        require(ps_contract == "preserved for phase-block QC; not used as a topology edge or lineage label",
                f"PS contract drifted: {sample}")

        enriched = dict(record)
        enriched["region_path"] = region_path
        enriched["layered_path"] = layered_path
        enriched["region_sha256"] = record["sha256"]["layered_region_view"]
        enriched["l3_status"] = l3["status"]
        enriched["canonical_record"] = record
        enriched_records.append(enriched)

    newest_source = datetime.fromtimestamp(max_source_mtime, tz=summary_generated.tzinfo)
    require(summary_generated >= newest_source,
            f"current summary is stale: generated={summary_generated.isoformat()} newest_source={newest_source.isoformat()}")
    validate_aggregate(canonical["aggregate"], enriched_records)

    for role, provenance in summary.get("code_provenance", {}).items():
        source_path = resolve_authority_path(provenance["path"])
        verify_bound_file(source_path, provenance["sha256"], f"code provenance {role}")

    comparison = summary.get("backbone_comparison", {})
    comparison_path = resolve_authority_path(comparison["path"])
    verify_bound_file(comparison_path, comparison["sha256"], "backbone comparison")
    require(comparison.get("aggregate", {}).get("verdict") == "backbone_sensitive",
            "backbone comparison verdict drifted")

    return {
        "summary": summary,
        "summary_path": CURRENT_SUMMARY,
        "summary_sha256": summary_sha256,
        "summary_generated": summary_generated,
        "canonical": canonical,
        "aggregate": canonical["aggregate"],
        "samples": enriched_records,
        "run_root": run_root,
        "success": success,
        "verification": verification,
        "verification_path": verification_path,
    }


class _MetaParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.values: dict[str, str] = {}

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "meta":
            return
        values = {key.lower(): value for key, value in attrs if value is not None}
        name = values.get("name")
        if name in PAGE_META_NAMES.values():
            self.values[name] = values.get("content", "")


def read_page_meta(path: Path) -> dict[str, str]:
    parser = _MetaParser()
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            parser.feed(handle.read(256 * 1024))
    except OSError as exc:
        raise AuthorityError(f"cannot inspect sample page {path}: {exc}") from exc
    return parser.values


def page_freshness(output: Path, sample: dict[str, Any], authority: dict[str, Any]) -> tuple[bool, str]:
    if not output.is_file():
        return False, "missing page"
    if output.stat().st_size < 1024:
        return False, "page is unexpectedly small"
    markers = read_page_meta(output)
    expected = {
        PAGE_META_NAMES["summary_sha256"]: authority["summary_sha256"],
        PAGE_META_NAMES["region_sha256"]: sample["region_sha256"],
        PAGE_META_NAMES["sample"]: sample["sample"],
    }
    for name, value in expected.items():
        if markers.get(name) != value:
            return False, f"marker {name} mismatch"
    newest_input = max(authority["summary_path"].stat().st_mtime, sample["region_path"].stat().st_mtime)
    if output.stat().st_mtime < newest_input:
        return False, "page mtime predates its summary or region view"
    return True, "hash-bound"


def collect_rows(authority: dict[str, Any], build_samples: bool) -> list[dict[str, Any]]:
    """Build or verify seven hash-bound sample pages without fallback data."""
    rows: list[dict[str, Any]] = []
    for sample in authority["samples"]:
        name = sample["sample"]
        output = OUTDIR / f"{name}.html"
        if build_samples:
            env = dict(
                os.environ,
                SM_RV=str(sample["region_path"]),
                SM_OUT=str(output),
                SM_CURRENT_SUMMARY=str(authority["summary_path"]),
                SM_CURRENT_SUMMARY_SHA256=authority["summary_sha256"],
                SM_CANONICAL_SAMPLE=json.dumps(
                    sample["canonical_record"], ensure_ascii=False, sort_keys=True, separators=(",", ":")
                ),
                SM_REGION_VIEW_SHA256=sample["region_sha256"],
            )
            result = subprocess.run(
                [PYTHON, str(BUILD)],
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            if result.returncode != 0:
                output_excerpt = result.stdout.decode(errors="replace")[-1200:]
                raise AuthorityError(f"sample build failed for {name} (rc={result.returncode}):\n{output_excerpt}")

        ready, reason = page_freshness(output, sample, authority)
        require(ready, f"{'built' if build_samples else 'existing'} sample page is stale for {name}: {reason}")
        row = dict(sample)
        row.update({"output": output, "page_ready": True, "page_size": output.stat().st_size})
        rows.append(row)
        print(f"[{name}] {'BUILT' if build_samples else 'VERIFIED'} hash-bound page")
    require(len(rows) == 7, "did not collect exactly seven canonical sample pages")
    return rows


def percentage(value: int, total: int) -> float:
    return 100.0 * value / total if total else 0.0


def topology_segments(topology: dict[str, Any], total: int) -> str:
    segments = []
    for key, short_label, long_label, class_name in TOPOLOGY_CLASSES:
        count = int(topology[key])
        share = percentage(count, total)
        segments.append(
            f'<span class="bar-segment segment-{class_name}" style="width:{share:.4f}%" '
            f'title="{escaped(short_label)} · {count:,} · {share:.1f}%">'
            f'<span class="sr-only">{escaped(long_label)} {count:,}，{share:.1f}%</span></span>'
        )
    return "".join(segments)


def topology_legend(aggregate: dict[str, Any]) -> str:
    total = int(aggregate["W_primary"])
    topology = aggregate["topology_classes"]
    items = []
    for key, short_label, long_label, class_name in TOPOLOGY_CLASSES:
        count = int(topology[key])
        items.append(
            f'<li><span class="legend-swatch segment-{class_name}" aria-hidden="true"></span>'
            f'<span><strong>{escaped(short_label)}</strong><small>{escaped(long_label)}</small></span>'
            f'<b>{count:,}<small>{percentage(count, total):.1f}%</small></b></li>'
        )
    return "".join(items)


def sample_topology_rows(rows: list[dict[str, Any]]) -> str:
    rendered = []
    for row in rows:
        total = int(row["W_primary"])
        topology = row["topology_classes"]
        label = ", ".join(
            f"{short} {int(topology[key]):,}"
            for key, short, _, _ in TOPOLOGY_CLASSES
        )
        rendered.append(
            f'<article class="sample-topology"><div class="sample-topology-head">'
            f'<a href="{escaped(row["sample"])}.html"><strong>{escaped(row["sample"])}</strong></a>'
            f'<span>W_primary {total:,}</span></div>'
            f'<div class="topology-bar" role="img" aria-label="{escaped(row["sample"])}：{escaped(label)}">'
            f'{topology_segments(topology, total)}</div>'
            f'<div class="bar-scale"><span>0</span><span>候選拓撲狀態</span><span>100%</span></div></article>'
        )
    return "".join(rendered)


def genome_launchers(rows: list[dict[str, Any]]) -> str:
    links = []
    for index, row in enumerate(rows, start=1):
        links.append(
            f'<a class="genome-link" href="{escaped(row["sample"])}.html">'
            f'<span class="launcher-index">{index:02d}</span>'
            f'<span><strong>{escaped(row["sample"])}</strong>'
            f'<small>chr1–22 · {int(row["W_tree"]):,} regions</small></span>'
            f'<span class="launcher-action">進入巡覽</span></a>'
        )
    return "".join(links)


def sample_position_leaders(rows: list[dict[str, Any]]) -> str:
    """Render per-dataset W_primary chromosome leaders from bound region views."""
    rendered = []
    for row in rows:
        region_view = load_json(row["region_path"])
        counts = {f"chr{number}": 0 for number in range(1, 23)}
        regions = region_view.get("regions") or []
        for region in regions:
            chrom = region.get("chrom")
            if chrom not in counts:
                continue
            if any(line.get("is_primary_lineage") is True for line in (region.get("lineages") or [])):
                counts[chrom] += 1
        observed_primary = sum(counts.values())
        expected_primary = int(row["W_primary"])
        require(
            observed_primary == expected_primary,
            f'{row["sample"]}: position leader W_primary mismatch ({observed_primary} != {expected_primary})',
        )
        leader, count = min(counts.items(), key=lambda item: (-item[1], int(item[0][3:])))
        rendered.append(
            f'<a href="{escaped(row["sample"])}.html#chr={escaped(leader)}">'
            f'<span>{escaped(row["sample"])}</span><b>{escaped(leader)} · {count:,}</b></a>'
        )
    return "".join(rendered)


def cohort_rows(rows: list[dict[str, Any]]) -> str:
    rendered = []
    for row in rows:
        complete = int(row["complete_regions"])
        primary = int(row["W_primary"])
        rendered.append(
            '<tr>'
            f'<th scope="row" class="sample-cell"><a href="{escaped(row["sample"])}.html">'
            f'<strong>{escaped(row["sample"])}</strong></a><small>canonical PASS</small></th>'
            f'<td class="numeric" data-label="LPS PASS">{int(row["tree_input_records"]):,}</td>'
            f'<td class="numeric" data-label="Retained sSNV">{int(row["retained_sSNV"]):,}</td>'
            f'<td class="numeric" data-label="W tree">{int(row["W_tree"]):,}</td>'
            f'<td class="numeric" data-label="W primary">{primary:,}</td>'
            f'<td class="numeric" data-label="Complete"><strong>{complete:,}</strong>'
            f'<small>{percentage(complete, primary):.1f}% of W_primary</small></td>'
            f'<td class="numeric" data-label="Incomplete">{int(row["incomplete_regions"]):,}</td>'
            f'<td class="numeric" data-label="Primary units">{int(row["primary_units"]):,}</td>'
            f'<td><a class="table-action" href="{escaped(row["sample"])}.html">開啟全基因頁</a></td>'
            '</tr>'
        )
    return "".join(rendered)


def evidence_links(authority: dict[str, Any], rows: list[dict[str, Any]]) -> str:
    summary_href = os.path.relpath(authority["summary_path"], OUTDIR)
    comparison_meta = authority["summary"]["backbone_comparison"]
    comparison_path = resolve_authority_path(comparison_meta["path"])
    comparison_href = os.path.relpath(comparison_path, OUTDIR)
    links = [
        f'<li><a href="{escaped(summary_href)}">Current canonical machine summary JSON</a>'
        f'<code>{escaped(authority["summary_sha256"])}</code></li>',
        f'<li><a href="{escaped(comparison_href)}">Backbone sensitivity comparison JSON</a>'
        f'<code>{escaped(comparison_meta["sha256"])}</code></li>',
    ]
    for row in rows:
        links.append(
            f'<li><a href="{escaped(row["region_path"].as_uri())}">{escaped(row["sample"])} region-view JSON</a>'
            f'<code>{escaped(row["region_sha256"])}</code></li>'
        )
    return "".join(links)


def build_index(authority: dict[str, Any], rows: list[dict[str, Any]]) -> str:
    aggregate = authority["aggregate"]
    comparison = authority["summary"]["backbone_comparison"]["aggregate"]
    generated_at = datetime.now().astimezone().isoformat(timespec="minutes")
    verified_at = parse_timestamp(authority["verification"]["verified_at_utc"], "verification").astimezone()
    run_id = authority["run_root"].name
    topology = aggregate["topology_classes"]
    exact_unique = int(topology["exact_and_topology_unique"])
    shape_unique_exact_multiple = int(topology["topology_unique_exact_multiple"])
    multi_shape = int(topology["topology_multiple_exact_multiple"])
    incomplete = int(topology["incomplete"])
    shape_unique = exact_unique + shape_unique_exact_multiple
    complete_rate = percentage(int(aggregate["complete_regions"]), int(aggregate["W_primary"]))
    incomplete_rate = percentage(int(aggregate["incomplete_regions"]), int(aggregate["W_primary"]))
    no_primary_rate = percentage(int(aggregate["no_primary_lineage"]), int(aggregate["W_tree"]))

    template = Template("""<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light">
<link rel="icon" href="data:,">
<meta name="intersubmod-current-summary-sha256" content="$summary_sha256">
<meta name="intersubmod-backbone-comparison-sha256" content="$comparison_sha256">
<meta name="intersubmod-canonical-run" content="$run_id">
<title>Layered reconstruction · 全基因 cohort command center</title>
<style>
:root{--ink:#122b32;--ink-soft:#345058;--muted:#596a6e;--paper:#f2efe7;--panel:#fffdf7;--line:#c8c6bc;--rail:#0f3138;--teal:#168b8d;--cyan:#53c4bd;--amber:#d47a2b;--brick:#b84f3d;--sand:#dfd5bc;--blue:#426f9e;--focus:#f2a900;--shadow:0 15px 45px rgba(19,45,52,.09)}
*{box-sizing:border-box}
html{background:var(--paper);scroll-behavior:smooth}
body{margin:0;color:var(--ink);background:linear-gradient(90deg,rgba(15,49,56,.035) 1px,transparent 1px),linear-gradient(180deg,#e8e4d9 0,var(--paper) 360px);background-size:48px 100%,100% 100%;font-family:"IBM Plex Sans","Noto Sans TC","PingFang TC","Microsoft JhengHei",sans-serif;line-height:1.58}
a{color:var(--rail);text-decoration-color:rgba(15,49,56,.42);text-underline-offset:3px}
a:hover{text-decoration-thickness:2px}
a:focus-visible,summary:focus-visible,.table-scroll:focus-visible{outline:3px solid var(--focus);outline-offset:3px}
.sr-only{position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0}
.skip-link{position:fixed;z-index:100;left:14px;top:10px;transform:translateY(-180%);padding:9px 13px;border:2px solid var(--ink);background:#fff;color:#111;font-weight:800}
.skip-link:focus{transform:none}
.shell{width:min(100%,1360px);margin:0 auto;padding:22px clamp(14px,3vw,42px) 54px}
.hero{position:relative;display:grid;grid-template-columns:minmax(0,1.55fr) minmax(290px,.75fr);overflow:hidden;border:1px solid #789095;background:var(--rail);color:#f7f3e9;box-shadow:var(--shadow)}
.hero:after{content:"CHR 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22";position:absolute;left:0;right:0;bottom:0;padding:7px 24px;border-top:1px solid rgba(255,255,255,.18);color:rgba(255,255,255,.55);font:600 9px/1.2 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.18em;word-spacing:.55em;white-space:nowrap;overflow:hidden}
.hero-main{padding:clamp(28px,5vw,66px) clamp(22px,5vw,64px) 60px}
.eyebrow{display:flex;flex-wrap:wrap;gap:8px;margin:0 0 20px;font:700 10px/1.2 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.11em;text-transform:uppercase}
.eyebrow span{padding:5px 8px;border:1px solid rgba(255,255,255,.28)}
.eyebrow .canonical{border-color:var(--cyan);background:rgba(83,196,189,.12);color:#aaf0e8}
h1{max-width:830px;margin:0;font-family:"Iowan Old Style","Noto Serif TC","Songti TC",serif;font-size:clamp(34px,5.3vw,72px);font-weight:600;letter-spacing:-.045em;line-height:.98}
.lede{max-width:780px;margin:22px 0 0;color:#cbd9d9;font-size:clamp(15px,1.6vw,19px)}
.hero-aside{display:flex;flex-direction:column;justify-content:space-between;padding:30px 28px 54px;border-left:1px solid rgba(255,255,255,.18);background:rgba(0,0,0,.12)}
.authority-label{margin:0;color:#9bb1b4;font:700 10px/1.2 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.1em;text-transform:uppercase}
.authority-value{display:block;margin-top:7px;overflow-wrap:anywhere;font:700 13px/1.45 "IBM Plex Mono","Cascadia Code",monospace}
.authority-list{display:grid;gap:17px;margin:0}.authority-list div{padding-bottom:14px;border-bottom:1px solid rgba(255,255,255,.13)}.authority-list div:last-child{border-bottom:0}
.sensitivity-banner{display:grid;grid-template-columns:minmax(235px,.8fr) minmax(0,1.7fr);gap:20px;margin-top:12px;padding:17px 19px;border:1px solid #c88a4f;border-left:7px solid var(--amber);background:#fff7e9}.sensitivity-banner h2{font-family:"IBM Plex Sans","Noto Sans TC",sans-serif;font-size:18px;font-weight:800;letter-spacing:0}.sensitivity-banner p{margin:6px 0 0;color:#664522;font-size:12px}.sensitivity-metrics{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:7px}.sensitivity-metric{padding:9px 10px;border:1px solid #dec4a5;background:#fffdf8}.sensitivity-metric span{display:block;color:#765737;font-size:9.5px}.sensitivity-metric strong{display:block;margin-top:5px;font:800 15px/1.1 "IBM Plex Mono","Cascadia Code",monospace}
.section{margin-top:34px}
.section-head{display:flex;align-items:end;justify-content:space-between;gap:24px;margin-bottom:13px;padding-bottom:10px;border-bottom:1px solid var(--line)}
.section-kicker{margin:0 0 4px;color:var(--teal);font:800 10px/1.2 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.12em;text-transform:uppercase}
h2{margin:0;font-family:"Iowan Old Style","Noto Serif TC","Songti TC",serif;font-size:clamp(23px,3vw,37px);font-weight:600;letter-spacing:-.025em;line-height:1.1}
.section-note{max-width:610px;margin:0;color:var(--muted);font-size:12.5px}
.metric-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:10px}
.metric-card{position:relative;min-height:170px;padding:20px;border:1px solid var(--line);background:var(--panel);box-shadow:0 6px 20px rgba(19,45,52,.035)}
.metric-card:before{content:"";position:absolute;left:-1px;right:-1px;top:-1px;height:4px;background:var(--teal)}
.metric-card:nth-child(2):before{background:var(--blue)}.metric-card:nth-child(3):before{background:var(--cyan)}.metric-card:nth-child(4):before{background:var(--brick)}
.metric-code{color:var(--muted);font:800 10px/1.2 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.1em;text-transform:uppercase}
.metric-number{display:block;margin:18px 0 7px;font:700 clamp(29px,3vw,43px)/1 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:-.05em;font-variant-numeric:tabular-nums}
.metric-card p{margin:0;color:var(--ink-soft);font-size:12.5px}
.metric-sub{display:block;margin-top:8px;color:var(--muted);font:700 10.5px/1.4 "IBM Plex Mono","Cascadia Code",monospace}
.topology-board{display:grid;grid-template-columns:minmax(0,1.55fr) minmax(320px,.75fr);gap:12px}
.topology-main,.topology-legend,.sample-topology,.dimension-card,.claim-card,.layer-card,.table-shell,.evidence-drawer{border:1px solid var(--line);background:var(--panel)}
.topology-main{padding:24px}.topology-main h3{margin:0 0 5px;font-size:16px}.topology-main>p{margin:0 0 24px;color:var(--muted);font-size:12.5px}
.topology-bar{display:flex;width:100%;height:30px;overflow:hidden;border:1px solid #829195;background:#e8e4d9}
.bar-segment{display:block;min-width:0;height:100%}.segment-unique{background:var(--teal)}.segment-shape{background:var(--blue)}.segment-multiple{background:var(--amber)}.segment-incomplete{background:var(--brick)}
.bar-scale{display:flex;justify-content:space-between;margin-top:6px;color:var(--muted);font:700 9.5px/1.2 "IBM Plex Mono","Cascadia Code",monospace;text-transform:uppercase}
.topology-legend{padding:18px}.topology-legend ul{display:grid;gap:12px;margin:0;padding:0;list-style:none}.topology-legend li{display:grid;grid-template-columns:10px 1fr auto;gap:10px;align-items:start}.legend-swatch{width:10px;height:30px}.topology-legend strong,.topology-legend small,.topology-legend b{display:block}.topology-legend strong{font:750 11px/1.3 "IBM Plex Mono","Cascadia Code",monospace}.topology-legend small{margin-top:2px;color:var(--muted);font-size:10.5px}.topology-legend b{text-align:right;font:750 13px/1.2 "IBM Plex Mono","Cascadia Code",monospace}.topology-legend b small{font-size:9.5px}
.sample-bars{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px;margin-top:12px}
.sample-topology{padding:15px 16px}.sample-topology-head{display:flex;justify-content:space-between;gap:12px;margin-bottom:9px}.sample-topology-head a{font-size:12.5px}.sample-topology-head span{color:var(--muted);font:700 10px/1.3 "IBM Plex Mono","Cascadia Code",monospace}.sample-topology .topology-bar{height:16px}
.dimension-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px}.dimension-card{display:flex;min-width:0;min-height:306px;flex-direction:column;padding:20px;border-top:4px solid var(--teal)}.dimension-card:nth-child(2){border-top-color:var(--amber)}.dimension-card:nth-child(3){border-top-color:var(--blue)}
.dimension-code{display:flex;justify-content:space-between;gap:10px;color:var(--teal);font:800 10px/1.2 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.07em;text-transform:uppercase}.dimension-card:nth-child(2) .dimension-code{color:var(--amber)}.dimension-card:nth-child(3) .dimension-code{color:var(--blue)}
.dimension-card h3{margin:18px 0 8px;font-size:17px}.dimension-answer{margin:0;color:var(--ink-soft);font-size:12.5px}.dimension-answer strong{color:var(--ink);font:800 22px/1.1 "IBM Plex Mono","Cascadia Code",monospace}.dimension-list{display:grid;gap:5px;margin:15px 0}.dimension-row{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:10px;padding:7px 8px;border:1px solid #ddd9cf;background:#f6f3eb}.dimension-row span{color:var(--muted);font-size:10.5px}.dimension-row b{font:750 10.5px/1.3 "IBM Plex Mono","Cascadia Code",monospace;text-align:right}.dimension-boundary{margin:auto 0 0;color:var(--muted);font-size:11.5px}.dimension-link{display:inline-block;margin-top:13px;padding:8px 10px;border:1px solid var(--rail);background:var(--rail);color:#fff;text-decoration:none;font-size:11px;font-weight:800}.dimension-link:hover{background:var(--teal)}
.position-leaders{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:4px;margin:0 0 15px}.position-leaders a{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:6px;padding:6px 7px;border:1px solid #ddd9cf;background:#f6f3eb;text-decoration:none}.position-leaders span{min-width:0;overflow-wrap:anywhere;color:var(--ink-soft);font-size:9.5px}.position-leaders b{font:750 9.5px/1.35 "IBM Plex Mono","Cascadia Code",monospace;white-space:nowrap}
.genome-launchers{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:9px}
.genome-link{display:grid;grid-template-columns:auto 1fr;gap:11px;align-items:start;min-height:122px;padding:16px;border:1px solid #98a5a4;background:var(--panel);color:var(--ink);text-decoration:none;transition:transform .18s ease,box-shadow .18s ease,border-color .18s ease}
.genome-link:hover{transform:translateY(-2px);border-color:var(--teal);box-shadow:0 10px 24px rgba(19,45,52,.09)}
.launcher-index{display:grid;place-items:center;width:29px;height:29px;border-radius:50%;background:var(--rail);color:#fff;font:700 10px/1 "IBM Plex Mono","Cascadia Code",monospace}.genome-link strong,.genome-link small{display:block}.genome-link strong{overflow-wrap:anywhere;font:800 13px/1.25 "IBM Plex Mono","Cascadia Code",monospace}.genome-link small{margin-top:6px;color:var(--muted);font-size:10.5px}.launcher-action{grid-column:2;margin-top:auto;color:var(--teal);font-size:11px;font-weight:800}
.table-shell{max-width:100%;box-shadow:var(--shadow)}.scroll-cue{display:flex;justify-content:space-between;gap:12px;padding:9px 12px;border-bottom:1px solid var(--line);background:#e9e7df;color:var(--muted);font-size:11px}.table-scroll{width:100%;overflow-x:auto;overscroll-behavior-inline:contain}.table-scroll:focus-visible{outline-offset:-3px}
table{width:100%;min-width:1020px;border-collapse:collapse;font-size:12px}th,td{padding:11px 10px;text-align:left;vertical-align:middle;border-bottom:1px solid #e2dfd6}thead th{position:sticky;top:0;z-index:2;background:#dfe5e3;color:#2f494f;font:800 9.5px/1.3 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.03em}tbody tr:last-child>*{border-bottom:0}tbody tr:hover>*{background:#f4f5ef}th:first-child,td:first-child{position:sticky;left:0;z-index:1}thead th:first-child{z-index:3}tbody th:first-child{background:var(--panel);box-shadow:7px 0 12px -13px #000}.sample-cell{min-width:166px}.sample-cell strong,.sample-cell small,td small{display:block}.sample-cell small,td small{margin-top:3px;color:var(--muted);font-size:10px}.numeric{text-align:right;font-family:"IBM Plex Mono","Cascadia Code",monospace;font-variant-numeric:tabular-nums}.table-action{display:inline-block;padding:7px 9px;border:1px solid var(--rail);background:var(--rail);color:#fff;text-decoration:none;white-space:nowrap;font-weight:750}.table-action:hover{background:var(--teal)}
.claim-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px}.claim-card{min-height:185px;padding:20px}.claim-card:nth-child(1){border-top:4px solid var(--teal)}.claim-card:nth-child(2){border-top:4px solid var(--amber)}.claim-card:nth-child(3){border-top:4px solid var(--brick)}.claim-index{color:var(--muted);font:800 10px/1 "IBM Plex Mono","Cascadia Code",monospace}.claim-card h3{margin:25px 0 9px;font-size:16px}.claim-card p{margin:0;color:var(--ink-soft);font-size:12.5px}
.layer-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:0;border:1px solid var(--line)}.layer-card{min-height:168px;padding:19px;border:0;border-right:1px solid var(--line)}.layer-card:last-child{border-right:0}.layer-code{color:var(--teal);font:850 11px/1 "IBM Plex Mono","Cascadia Code",monospace}.layer-card h3{margin:17px 0 7px;font-size:14px}.layer-card p{margin:0;color:var(--muted);font-size:11.5px}.layer-card.pending{background:#f8eee9}.layer-card.pending .layer-code{color:var(--brick)}
.qc-note{display:grid;grid-template-columns:auto 1fr;gap:14px;align-items:start;margin-top:10px;padding:15px 18px;border-left:5px solid var(--blue);background:#e8edf0}.qc-note strong{font:800 11px/1.3 "IBM Plex Mono","Cascadia Code",monospace}.qc-note p{margin:0;color:var(--ink-soft);font-size:12px}
.evidence-drawer{margin-top:34px}.evidence-drawer summary{cursor:pointer;padding:14px 17px;font-weight:800}.evidence-drawer[open] summary{border-bottom:1px solid var(--line)}.evidence-inner{padding:16px 18px}.evidence-inner p{margin:0 0 12px;color:var(--muted);font-size:11.5px}.evidence-list{margin:0;padding-left:20px}.evidence-list li+li{margin-top:9px}.evidence-list a{font-size:11.5px}.evidence-list code{display:block;margin-top:2px;overflow-wrap:anywhere;color:var(--muted);font-size:9.5px}
footer{display:flex;justify-content:space-between;gap:20px;margin-top:18px;padding-top:13px;border-top:1px solid var(--line);color:var(--muted);font-size:10.5px}.mono{font-family:"IBM Plex Mono","Cascadia Code",monospace;overflow-wrap:anywhere}
@media(max-width:1050px){.metric-grid{grid-template-columns:repeat(2,1fr)}.genome-launchers{grid-template-columns:repeat(3,1fr)}.topology-board{grid-template-columns:1fr}.dimension-grid{grid-template-columns:1fr}.dimension-card{min-height:0}.layer-grid{grid-template-columns:repeat(2,1fr)}.layer-card:nth-child(2){border-right:0}.layer-card:nth-child(-n+2){border-bottom:1px solid var(--line)}}
@media(max-width:760px){.shell{padding:12px 11px 36px}.hero{grid-template-columns:1fr}.hero-aside{border-left:0;border-top:1px solid rgba(255,255,255,.18)}.sensitivity-banner{grid-template-columns:1fr}.sensitivity-metrics{grid-template-columns:1fr 1fr}.section-head{display:block}.section-note{margin-top:7px}.sample-bars,.claim-grid{grid-template-columns:1fr}.genome-launchers{grid-template-columns:repeat(2,1fr)}footer{display:block}footer span{display:block;margin-top:5px}}
@media(max-width:500px){.metric-grid,.genome-launchers,.layer-grid{grid-template-columns:1fr}.hero-main{padding:28px 20px 58px}.hero-aside{padding:24px 20px 48px}.metric-card{min-height:0}.position-leaders{grid-template-columns:1fr}.layer-card{min-height:0;border-right:0;border-bottom:1px solid var(--line)}.layer-card:last-child{border-bottom:0}.scroll-cue span:last-child{display:none}.topology-main{padding:18px}.topology-legend{padding:15px}}
@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}.genome-link{transition:none}}
@media print{body{background:#fff}.shell{width:100%;padding:0}.hero,.table-shell{box-shadow:none}.genome-link{break-inside:avoid}.table-scroll{overflow:visible}table{min-width:0;font-size:8px}.skip-link,.scroll-cue,.table-action{display:none}}
</style>
</head>
<body>
<a class="skip-link" href="#main-content">跳至主要內容</a>
<main class="shell" id="main-content">
  <header class="hero">
    <div class="hero-main">
      <p class="eyebrow"><span class="canonical">Canonical / 7 of 7 PASS</span><span>chr1–22 全基因範圍</span></p>
      <h1>分層重建<br>全基因指揮中心</h1>
      <p class="lede">以 LongPhase-S recalibrated FILTER=PASS 為唯一主骨幹，從 cohort 總覽進入 7 個 dataset 的全常染色體 region 巡覽；每個數字保留明確 grain、分母與候選完整性。</p>
    </div>
    <aside class="hero-aside" aria-label="Canonical authority">
      <p class="authority-label">Machine-bound authority</p>
      <dl class="authority-list">
        <div><dt class="authority-label">Run</dt><dd class="authority-value">$run_id</dd></div>
        <div><dt class="authority-label">Verified</dt><dd class="authority-value">$verified_at</dd></div>
        <div><dt class="authority-label">Summary SHA-256</dt><dd class="authority-value">$summary_short</dd></div>
      </dl>
    </aside>
  </header>

  <aside class="sensitivity-banner" aria-labelledby="sensitivity-title">
    <div><p class="section-kicker">Backbone sensitivity · 不可混合分母</p><h2 id="sensitivity-title">$backbone_verdict</h2><p>LongPhase-S recalibrated FILTER=PASS 仍是 canonical 主結果；ClairS FILTER=PASS 僅作 sensitivity。低 overlap 指標表示所有解讀都必須標示 backbone-sensitive。</p></div>
    <div class="sensitivity-metrics" aria-label="Backbone sensitivity minimum concordance metrics">
      <div class="sensitivity-metric"><span>Retained-position Jaccard · min</span><strong>$retained_jaccard</strong></div>
      <div class="sensitivity-metric"><span>Primary-unit-key Jaccard · min</span><strong>$primary_jaccard</strong></div>
      <div class="sensitivity-metric"><span>Shared topology digest · min</span><strong>$topology_concordance</strong></div>
      <div class="sensitivity-metric"><span>Max absolute proportion delta</span><strong>$max_delta_pp pp</strong></div>
    </div>
  </aside>

  <section class="section" aria-labelledby="pulse-title">
    <div class="section-head"><div><p class="section-kicker">Cohort pulse</p><h2 id="pulse-title">先鎖定 W 與候選完整性</h2></div><p class="section-note">W_tree 是 retained multilocus regions；C/Topo 只以具有 mutation-bearing HP1/HP2 primary units 的 W_primary 為分母。</p></div>
    <div class="metric-grid">
      <article class="metric-card"><span class="metric-code">W_tree · region</span><strong class="metric-number">$w_tree</strong><p>chr1–22 retained regional groups。</p><small class="metric-sub">no-primary $no_primary · $no_primary_rate%</small></article>
      <article class="metric-card"><span class="metric-code">W_primary · region</span><strong class="metric-number">$w_primary</strong><p>至少一個 mutation-bearing HP1/HP2 primary unit。</p><small class="metric-sub">C / Topo denominator</small></article>
      <article class="metric-card"><span class="metric-code">Complete · region</span><strong class="metric-number">$complete</strong><p>所有 primary units 均 non-capped、full-pass、candidate-complete。</p><small class="metric-sub">$complete_rate% of W_primary</small></article>
      <article class="metric-card"><span class="metric-code">Incomplete · region</span><strong class="metric-number">$incomplete</strong><p>任一 primary unit 未完成時，C 與 Topo 必須保持 unavailable。</p><small class="metric-sub">$incomplete_rate% of W_primary</small></article>
    </div>
  </section>

  <section class="section" aria-labelledby="dimensions-title">
    <div class="section-head"><div><p class="section-kicker">Three reading dimensions</p><h2 id="dimensions-title">先用三個問題定位每個數字</h2></div><p class="section-note">舊頁的三軸教學結構保留；所有名詞與分母改綁 canonical v5，不沿用 clone-first 類別或 archive 數字。</p></div>
    <div class="dimension-grid">
      <article class="dimension-card"><div class="dimension-code"><span>01 · 拓樸型態</span><span>Topology shape</span></div><h3>分子累積形狀有幾種？</h3><p class="dimension-answer"><strong>$shape_unique</strong> 個 complete regions 的 Topo=1。</p><div class="dimension-list"><div class="dimension-row"><span>形狀唯一 / complete</span><b>$shape_unique / $complete · $shape_unique_rate%</b></div><div class="dimension-row"><span>多種形狀相容</span><b>$multi_shape · $multi_shape_rate%</b></div></div><p class="dimension-boundary">Topo 比較無節點標籤的 regional mutation-state shape；不等於 biological ancestry 或真實時間順序。</p></article>
      <article class="dimension-card"><div class="dimension-code"><span>02 · 可辨識度</span><span>Determinacy</span></div><h3>目前能唯一辨識到哪一層？</h3><p class="dimension-answer"><strong>$exact_unique</strong> 個 regions 可辨識到 exact candidate。</p><div class="dimension-list"><div class="dimension-row"><span>Exact 唯一</span><b>$exact_unique</b></div><div class="dimension-row"><span>只到 shape 唯一</span><b>$shape_only</b></div><div class="dimension-row"><span>Multi-shape</span><b>$multi_shape</b></div><div class="dimension-row"><span>尚未評估</span><b>$incomplete</b></div></div><p class="dimension-boundary">Incomplete 是 candidate set 未完整，所以 C／Topo 不可評估；不是 0，也不是已證明多解。</p></article>
      <article class="dimension-card"><div class="dimension-code"><span>03 · 基因體位置</span><span>Genome position</span></div><h3>這些 regions 位於哪裡？</h3><p class="dimension-answer"><strong>7 × chr1–22</strong> 的 canonical 全基因巡覽。</p><div class="dimension-list"><div class="dimension-row"><span>W_tree regions</span><b>$w_tree</b></div><div class="dimension-row"><span>Dataset pages</span><b>7</b></div><div class="dimension-row"><span>Autosomes / page</span><b>22</b></div></div><div class="position-leaders" aria-label="各 dataset 的 W primary count leader chromosome">$position_leaders</div><p class="dimension-boundary">上列為各 dataset 的 W_primary count leader；點選可直達該染色體。位置只做 chr:start–end 與 region 分布巡覽；未校正 chromosome 長度或輸入密度時，不作 hotspot／enrichment 判定。</p><a class="dimension-link" href="#launch-title">選擇 dataset 看位置</a></article>
    </div>
  </section>

  <section class="section" aria-labelledby="topology-title">
    <div class="section-head"><div><p class="section-kicker">Candidate topology</p><h2 id="topology-title">Exact 組合與拓撲形狀分開讀</h2></div><p class="section-note">C_region 是 exact candidate-tree 組合乘積；Topo_region 是 exact unlabeled-shape 組合乘積。Incomplete 不以 0 代替。</p></div>
    <div class="topology-board">
      <article class="topology-main"><h3>全 cohort · W_primary $w_primary</h3><p>四類互斥且守恆；不可能狀態 C=1 / Topo&gt;1 為 0。</p><div class="topology-bar" role="img" aria-label="全 cohort 候選拓撲四類">$aggregate_segments</div><div class="bar-scale"><span>0</span><span>W_primary composition</span><span>100%</span></div></article>
      <aside class="topology-legend" aria-label="Topology legend"><ul>$topology_legend</ul></aside>
    </div>
    <div class="sample-bars" aria-label="7 dataset topology distributions">$sample_topology_rows</div>
  </section>

  <section class="section" aria-labelledby="launch-title">
    <div class="section-head"><div><p class="section-kicker">Genome launchpad</p><h2 id="launch-title">進入 chr1–22 全基因巡覽</h2></div><p class="section-note">每個入口都已由 current summary SHA 與 sample region-view SHA 雙重綁定；沒有 freshness 證據的頁面不會出現在本 index。</p></div>
    <nav class="genome-launchers" aria-label="開啟各 dataset 全基因頁">$genome_launchers</nav>
  </section>

  <section class="section" id="cohort-table" aria-labelledby="table-title">
    <div class="section-head"><div><p class="section-kicker">Cohort ledger</p><h2 id="table-title">7 dataset 簡潔盤點</h2></div><p class="section-note">LPS PASS 是全 contig tree-input records；retained sSNV 與 W metrics 是 chr1–22。不同 grain 不直接互除。</p></div>
    <div class="table-shell"><div class="scroll-cue"><span>表格可水平捲動；樣本欄固定。</span><span>Canonical main only</span></div><div class="table-scroll" role="region" aria-label="Canonical cohort table" tabindex="0"><table>
      <caption class="sr-only">Canonical v5 七個 dataset 的輸入、retained sSNV、W_tree、W_primary、候選完整性與 primary units</caption>
      <thead><tr><th scope="col">Dataset</th><th scope="col" class="numeric">LPS PASS<br>records</th><th scope="col" class="numeric">Retained<br>sSNV</th><th scope="col" class="numeric">W_tree<br>regions</th><th scope="col" class="numeric">W_primary<br>regions</th><th scope="col" class="numeric">Complete<br>regions</th><th scope="col" class="numeric">Incomplete<br>regions</th><th scope="col" class="numeric">Primary<br>units</th><th scope="col">巡覽</th></tr></thead>
      <tbody>$cohort_rows</tbody>
    </table></div></div>
  </section>

  <section class="section" aria-labelledby="boundary-title">
    <div class="section-head"><div><p class="section-kicker">Claim boundary</p><h2 id="boundary-title">觀察、推論、限制分開</h2></div><p class="section-note">這是 regional mutation-state candidate sets；不是細胞層級真相或已確認的全腫瘤演化關係。</p></div>
    <div class="claim-grid">
      <article class="claim-card"><span class="claim-index">01 / 回答</span><h3>哪些區域候選與資料相容？</h3><p>在每個 mutation-bearing HP1/HP2 primary unit 內，保留全部 analysis-complete minimal candidates，並區分 exact uniqueness 與 topology uniqueness。</p></article>
      <article class="claim-card"><span class="claim-index">02 / 證據</span><h3>read-level sSNV 共現</h3><p>LongPhase-S recalibrated PASS 提供 tree backbone；L0 family partition 與 L1 read-state constraints 承重，CN 只作 post-tree context。</p></article>
      <article class="claim-card"><span class="claim-index">03 / 限制</span><h3>不可越過資料識別度</h3><p>單一 bulk 的 regional candidates 不能單獨證明細胞比例、真實 ancestry 或正交生物學確認；多解與 incomplete 都是應保留的結果。</p></article>
    </div>
  </section>

  <section class="section" aria-labelledby="layers-title">
    <div class="section-head"><div><p class="section-kicker">Evidence stack</p><h2 id="layers-title">L0 → L3 的承重順序</h2></div><p class="section-note">後續 auxiliary layer 不得回頭重寫 read-level observation，也不得偷偷選樹。</p></div>
    <div class="layer-grid">
      <article class="layer-card"><span class="layer-code">L0</span><h3>HP family partition</h3><p>HP1/HP2 mutation-bearing units 是 primary；H3、unphased 與 reference-only 都排除於 primary denominator。</p></article>
      <article class="layer-card"><span class="layer-code">L1</span><h3>Candidate enumeration</h3><p>以 read-state constraints 枚舉 minimal candidates；capped 代表候選集合未完成。</p></article>
      <article class="layer-card"><span class="layer-code">L2</span><h3>Copy-number context</h3><p>只對 recurrence-like outcome 提供 post-tree annotation；不重新定義 observed state。</p></article>
      <article class="layer-card pending"><span class="layer-code">L3 · NOT EVALUATED</span><h3>Bounded methyl auxiliary</h3><p>本 canonical run 未評估；未來僅限 negative screen / residual flag，禁止排序或確認候選。</p></article>
    </div>
    <aside class="qc-note"><strong>PS / QC ONLY</strong><p>Phase-set 只保留作 phase-block 品質脈絡；不是 topology edge，也不是 lineage label。全 cohort 有 $mixed_ps 個 multi-PS regions。</p></aside>
  </section>

  <details class="evidence-drawer">
    <summary>機器證據與原始 JSON（預設收合）</summary>
    <div class="evidence-inner"><p>下列連結只供 provenance readback；正文所有數字均由 hash-verified canonical machine summary 產生。</p><ul class="evidence-list">$evidence_links</ul></div>
  </details>

  <footer><span>Canonical main · LongPhase-S recalibrated FILTER=PASS · chr1–22 · 7/7 PASS</span><span class="mono">Page generated $generated_at</span></footer>
</main>
</body>
</html>
""")
    return template.substitute(
        summary_sha256=escaped(authority["summary_sha256"]),
        comparison_sha256=escaped(authority["summary"]["backbone_comparison"]["sha256"]),
        summary_short=escaped(authority["summary_sha256"][:16] + "…"),
        run_id=escaped(run_id),
        backbone_verdict=escaped(str(comparison["verdict"]).replace("_", " ").upper()),
        retained_jaccard=f'{float(comparison["min_retained_position_jaccard"]):.3f}',
        primary_jaccard=f'{float(comparison["min_primary_unit_key_jaccard"]):.3f}',
        topology_concordance=f'{float(comparison["min_shared_topology_digest_concordance"]):.3f}',
        max_delta_pp=f'{float(comparison["max_abs_delta_pp"]):.2f}',
        verified_at=escaped(verified_at.isoformat(timespec="minutes")),
        generated_at=escaped(generated_at),
        w_tree=f'{int(aggregate["W_tree"]):,}',
        w_primary=f'{int(aggregate["W_primary"]):,}',
        no_primary=f'{int(aggregate["no_primary_lineage"]):,}',
        no_primary_rate=f'{no_primary_rate:.1f}',
        complete=f'{int(aggregate["complete_regions"]):,}',
        incomplete=f'{int(aggregate["incomplete_regions"]):,}',
        complete_rate=f'{complete_rate:.1f}',
        incomplete_rate=f'{incomplete_rate:.1f}',
        exact_unique=f'{exact_unique:,}',
        shape_only=f'{shape_unique_exact_multiple:,}',
        multi_shape=f'{multi_shape:,}',
        shape_unique=f'{shape_unique:,}',
        shape_unique_rate=f'{percentage(shape_unique, int(aggregate["complete_regions"])):.1f}',
        multi_shape_rate=f'{percentage(multi_shape, int(aggregate["complete_regions"])):.1f}',
        mixed_ps=f'{int(aggregate["mixed_PS_regions"]):,}',
        aggregate_segments=topology_segments(topology, int(aggregate["W_primary"])),
        topology_legend=topology_legend(aggregate),
        sample_topology_rows=sample_topology_rows(rows),
        position_leaders=sample_position_leaders(rows),
        genome_launchers=genome_launchers(rows),
        cohort_rows=cohort_rows(rows),
        evidence_links=evidence_links(authority, rows),
    )


def write_index_atomic(document: str) -> Path:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    output = OUTDIR / "index.html"
    staging = OUTDIR / f".index.html.stage.{os.getpid()}"
    with staging.open("w", encoding="utf-8") as handle:
        handle.write(document)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(staging, output)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--index-only",
        action="store_true",
        help="Verify seven existing hash-bound sample pages and rebuild only index.html.",
    )
    args = parser.parse_args()
    try:
        authority = load_authority()
        rows = collect_rows(authority, build_samples=not args.index_only)
        index_path = write_index_atomic(build_index(authority, rows))
    except AuthorityError as exc:
        print(f"FAIL CLOSED: {exc}", file=sys.stderr)
        return 2
    print(f"OK canonical index + 7/7 hash-bound sample pages -> {index_path}")
    print(f"SUMMARY_SHA256={authority['summary_sha256']}")
    print(f"CANONICAL_RUN={authority['run_root']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
