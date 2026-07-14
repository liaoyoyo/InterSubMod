#!/usr/bin/env python3
"""
Build one canonical-v5 layered reconstruction dataset workstation.

This renderer is intentionally fail-closed:
  - the current machine summary is the only authority for dataset paths/counts;
  - the region-view SHA-256 must match the canonical record;
  - region-level C/Topo is recomputed and reconciled before HTML is written.

Inputs (environment):
  SM_RV                       canonical layered_region_view JSON
  SM_OUT                      output standalone HTML
  SM_CURRENT_SUMMARY          current machine-summary JSON
  SM_CANONICAL_SAMPLE         canonical sample record encoded as JSON (optional)
  SM_REGION_VIEW_SHA256       expected region-view SHA-256 (optional)

Output:
  one offline HTML with a lightweight all-genome index and lazy JSON chunks for
  chr1-22 detail. Raw JSON links are present only inside a collapsed evidence
  drawer.
"""

from __future__ import annotations

import hashlib
import html
import json
import math
import os
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


DATASET_ORDER = [
    "HCC1395",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
]

BIOLOGICAL_SAMPLE = {
    "HCC1395": "HCC1395",
    "HCC1395_DORADO": "HCC1395",
    "COLO829": "COLO829",
    "H1437": "H1437",
    "H2009": "H2009",
    "HCC1937": "HCC1937",
    "HCC1954": "HCC1954",
}

TOPOLOGY_KEYS = (
    "exact_and_topology_unique",
    "topology_unique_exact_multiple",
    "topology_multiple_exact_multiple",
    "incomplete",
)


def fail(message: str) -> None:
    raise SystemExit(f"CANONICAL WORKSTATION BUILD FAILED: {message}")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_text(value) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")


def product(values) -> int:
    result = 1
    for value in values:
        result *= int(value)
    return result


def primary_lineages(region):
    return [line for line in region.get("lineages", []) if line.get("is_primary_lineage") is True]


def topology_profile(region):
    primary = primary_lineages(region)
    auxiliary = [line for line in region.get("lineages", []) if not line.get("is_primary_lineage")]
    # This facet means the solver required recurrence for a primary unit.  A
    # capped display prefix may itself contain recurrent edges, but that is not
    # equivalent to a recurrence-required classification and must not create a
    # false-positive region badge.
    has_recurrence = any(
        str(line.get("L1_class", "")).startswith("recurrence")
        for line in primary
    )
    if not primary:
        return {
            "topology_class": "no_primary_lineage",
            "candidate_complete": None,
            "C": None,
            "Topo": None,
            "n_primary": 0,
            "n_auxiliary": len(auxiliary),
            "evidence_mode": "not_applicable",
            "has_recurrence": has_recurrence,
            "hidden_positive": any(int(line.get("n_hidden") or 0) > 0 for line in auxiliary),
        }

    complete = all(
        line.get("analysis_candidate_set_complete") is True
        and line.get("verification_status") == "full_pass"
        for line in primary
    )
    full_units = sum(1 for line in primary if int(line.get("n_full_pops") or 0) > 0)
    partial_units = sum(1 for line in primary if int(line.get("n_partial") or 0) > 0)
    if full_units and partial_units:
        evidence_mode = "full_and_partial"
    elif partial_units:
        evidence_mode = "partial_only"
    elif full_units:
        evidence_mode = "full_only"
    else:
        evidence_mode = "no_read_state_groups"

    if not complete:
        topology_class = "incomplete"
        joint_c = None
        joint_topology = None
    else:
        joint_c = product(line.get("n_trees") for line in primary)
        joint_topology = product(line.get("n_distinct_shapes_exact") for line in primary)
        if not (joint_c >= joint_topology >= 1):
            fail(
                f"{region.get('region')}: invariant C>=Topo>=1 failed "
                f"(C={joint_c}, Topo={joint_topology})"
            )
        if joint_c == 1 and joint_topology == 1:
            topology_class = "exact_and_topology_unique"
        elif joint_c > 1 and joint_topology == 1:
            topology_class = "topology_unique_exact_multiple"
        elif joint_c > 1 and joint_topology > 1:
            topology_class = "topology_multiple_exact_multiple"
        else:
            fail(
                f"{region.get('region')}: impossible exact-unique/topology-multiple state "
                f"(C={joint_c}, Topo={joint_topology})"
            )

    return {
        "topology_class": topology_class,
        "candidate_complete": complete,
        "C": joint_c,
        "Topo": joint_topology,
        "n_primary": len(primary),
        "n_auxiliary": len(auxiliary),
        "evidence_mode": evidence_mode,
        "has_recurrence": has_recurrence,
        "hidden_positive": any(int(line.get("n_hidden") or 0) > 0 for line in primary),
    }


def compact_index(region, chunk_index, profile):
    primary = primary_lineages(region)
    return {
        "region": region["region"],
        "chrom": region["chrom"],
        "start": int(region["start"]),
        "end": int(region["end"]),
        "n_sSNV": int(region.get("n_sSNV") or 0),
        "endpoint_distance_bp": int(region["end"]) - int(region["start"]),
        "inclusive_span_bp": int(region["end"]) - int(region["start"]) + 1,
        "cn": region.get("cn", "unavailable"),
        "hp_multiplicity": int(region.get("hp_multiplicity") or 0),
        "n_primary": profile["n_primary"],
        "n_auxiliary": profile["n_auxiliary"],
        "topology_class": profile["topology_class"],
        "candidate_complete": profile["candidate_complete"],
        "C": profile["C"],
        "Topo": profile["Topo"],
        "evidence_mode": profile["evidence_mode"],
        "has_recurrence": profile["has_recurrence"],
        "hidden_positive": profile["hidden_positive"],
        "verification_full_pass": bool(primary)
        and all(line.get("verification_status") == "full_pass" for line in primary),
        "chunk_index": chunk_index,
    }


def validate_recomputed(sample_record, indices):
    counts = Counter(row["topology_class"] for row in indices)
    expected = sample_record["topology_classes"]
    for key in TOPOLOGY_KEYS:
        if counts[key] != int(expected[key]):
            fail(
                f"{sample_record['sample']}: recomputed {key}={counts[key]} "
                f"does not match canonical {expected[key]}"
            )
    if counts["no_primary_lineage"] != int(sample_record["no_primary_lineage"]):
        fail(
            f"{sample_record['sample']}: no_primary_lineage={counts['no_primary_lineage']} "
            f"does not match canonical {sample_record['no_primary_lineage']}"
        )
    if len(indices) != int(sample_record["W_tree"]):
        fail(f"{sample_record['sample']}: region count does not equal W_tree")
    if len(indices) - counts["no_primary_lineage"] != int(sample_record["W_primary"]):
        fail(f"{sample_record['sample']}: reconstructed W_primary mismatch")


def main() -> None:
    rv_value = os.environ.get("SM_RV")
    out_value = os.environ.get("SM_OUT")
    summary_value = os.environ.get("SM_CURRENT_SUMMARY")
    if not rv_value or not out_value or not summary_value:
        fail("SM_RV, SM_OUT and SM_CURRENT_SUMMARY are required")

    rv_path = Path(rv_value).resolve()
    out_path = Path(out_value).resolve()
    summary_path = Path(summary_value).resolve()
    if not rv_path.is_file() or not summary_path.is_file():
        fail("canonical region-view or machine-summary file is missing")

    summary_sha = file_sha256(summary_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("all_pass") is not True:
        fail("machine summary all_pass is not true")
    canonical = summary.get("canonical") or {}
    if canonical.get("label") != "longphase_s_recalibrated_FILTER_PASS":
        fail(f"unexpected canonical label: {canonical.get('label')}")
    aggregate = canonical.get("aggregate") or {}
    if int(aggregate.get("dataset_count", -1)) != 7:
        fail("canonical dataset_count must equal 7")
    sensitivity = summary.get("sensitivity") or {}
    comparison_meta = summary.get("backbone_comparison") or {}
    comparison = comparison_meta.get("aggregate") or {}
    if sensitivity.get("label") != "clairs_FILTER_PASS_sensitivity":
        fail(f"unexpected sensitivity label: {sensitivity.get('label')}")
    if comparison.get("verdict") != "backbone_sensitive":
        fail(f"unexpected backbone comparison verdict: {comparison.get('verdict')}")
    for metric in (
        "min_retained_position_jaccard",
        "min_primary_unit_key_jaccard",
        "min_shared_topology_digest_concordance",
    ):
        value = float(comparison.get(metric, -1))
        if not 0 <= value <= 1:
            fail(f"invalid backbone comparison metric {metric}={value}")
    if float(comparison.get("max_abs_delta_pp", -1)) < 0:
        fail("invalid backbone comparison max_abs_delta_pp")
    comparison_path = Path(comparison_meta.get("path", ""))
    if not comparison_path.is_absolute():
        comparison_path = (summary_path.parents[2] / comparison_path).resolve()
    if not comparison_path.is_file():
        fail(f"missing backbone comparison file: {comparison_path}")
    comparison_sha = file_sha256(comparison_path)
    if comparison_sha != comparison_meta.get("sha256"):
        fail("backbone comparison SHA-256 mismatch")
    comparison_payload = json.loads(comparison_path.read_text(encoding="utf-8"))

    requested_record = None
    if os.environ.get("SM_CANONICAL_SAMPLE"):
        requested_record = json.loads(os.environ["SM_CANONICAL_SAMPLE"])
    region_view = json.loads(rv_path.read_text(encoding="utf-8"))
    sample = region_view.get("sample")
    sample_record = next(
        (record for record in canonical.get("samples", []) if record.get("sample") == sample),
        None,
    )
    if not sample_record:
        fail(f"{sample}: no canonical sample record")
    sample_comparison = next(
        (record for record in comparison_payload.get("comparisons", []) if record.get("sample") == sample),
        None,
    )
    if not sample_comparison:
        fail(f"{sample}: no backbone comparison record")
    if sample_comparison.get("verdict") not in {"robust", "moderately_sensitive", "backbone_sensitive"}:
        fail(f"{sample}: invalid backbone comparison verdict")
    for metric in (
        "retained_position_jaccard",
        "primary_unit_key_jaccard",
        "shared_unit_topology_digest_concordance",
    ):
        value = float(sample_comparison.get(metric, -1))
        if not 0 <= value <= 1:
            fail(f"{sample}: invalid backbone comparison metric {metric}={value}")
    delta_values = list(sample_comparison.get("delta", {}).values())
    if not delta_values:
        fail(f"{sample}: backbone comparison delta metrics missing")
    sample_max_delta = max(abs(float(value)) for value in delta_values)
    if requested_record is not None and requested_record != sample_record:
        fail(f"{sample}: SM_CANONICAL_SAMPLE differs from machine summary")
    if sample_record.get("pass") is not True or sample_record.get("all_invariants_pass") is not True:
        fail(f"{sample}: canonical verifier did not pass")

    recorded_path = Path(sample_record["paths"]["layered_region_view"]).resolve()
    if recorded_path != rv_path:
        fail(f"{sample}: SM_RV is not the canonical region-view path")
    actual_rv_sha = file_sha256(rv_path)
    expected_rv_sha = os.environ.get("SM_REGION_VIEW_SHA256") or sample_record["sha256"][
        "layered_region_view"
    ]
    if actual_rv_sha != expected_rv_sha or actual_rv_sha != sample_record["sha256"][
        "layered_region_view"
    ]:
        fail(f"{sample}: region-view SHA-256 mismatch")

    census = region_view.get("census") or {}
    if census.get("U1_backbone_source") != "longphase_s_recalibrated_filter_pass":
        fail(f"{sample}: region-view backbone is not LongPhase-S recalibrated FILTER=PASS")
    if census.get("analysis_scope") != "chr1-22 primary; chrX/chrY out-of-scope census only":
        fail(f"{sample}: unexpected analysis scope")
    if region_view.get("analysis_contract", {}).get("primary_lineage") != "mutation-bearing HP1/HP2 only":
        fail(f"{sample}: primary-lineage contract mismatch")
    l3 = region_view.get("L3_methyl") or {}
    if l3.get("status") != "not_evaluated" or l3.get("role") != "bounded_auxiliary":
        fail(f"{sample}: L3 contract mismatch")

    regions = region_view.get("regions") or []
    by_chrom = defaultdict(list)
    indices = []
    for region in regions:
        chrom = region.get("chrom")
        if chrom not in {f"chr{i}" for i in range(1, 23)}:
            fail(f"{sample}: out-of-scope region entered display payload: {region.get('region')}")
        profile = topology_profile(region)
        detail = dict(region)
        detail["_ui"] = profile
        chunk_index = len(by_chrom[chrom])
        by_chrom[chrom].append(detail)
        indices.append(compact_index(region, chunk_index, profile))
    validate_recomputed(sample_record, indices)

    chunk_scripts = []
    chunk_manifest = {}
    for chromosome_number in range(1, 23):
        chrom = f"chr{chromosome_number}"
        payload = json_text(by_chrom.get(chrom, []))
        digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
        chunk_manifest[chrom] = {
            "count": len(by_chrom.get(chrom, [])),
            "sha256": digest,
        }
        chunk_scripts.append(
            '<script type="application/json" id="chunk-{}" '
            'data-sha256="{}">{}</script>'.format(chrom, digest, payload)
        )

    page_data = {
        "sample": sample,
        "biological_sample": BIOLOGICAL_SAMPLE[sample],
        "dataset_order": DATASET_ORDER,
        "canonical_sample": sample_record,
        "canonical_aggregate": aggregate,
        "sensitivity": {
            "label": sensitivity.get("label"),
            "run_root": sensitivity.get("run_root"),
            "success_sha256": sensitivity.get("success_sha256"),
            "verification_sha256": sensitivity.get("verification_sha256"),
        },
        "backbone_comparison": {
            "aggregate": comparison,
            "sample": sample_comparison,
            "path": str(comparison_path),
            "sha256": comparison_sha,
        },
        "summary": {
            "schema_name": summary.get("schema_name"),
            "schema_version": summary.get("schema_version"),
            "generated_at": summary.get("generated_at"),
            "claim_scope": summary.get("claim_scope"),
            "canonical_label": canonical.get("label"),
            "run_root": canonical.get("run_root"),
            "success_sha256": canonical.get("success_sha256"),
            "verification_sha256": canonical.get("verification_sha256"),
            "all_pass": summary.get("all_pass"),
        },
        "analysis_contract": region_view.get("analysis_contract"),
        "copy_number_contract": region_view.get("copy_number_contract"),
        "l3": l3,
        "read_tag_census": census.get("read_tag_census"),
        "census": census,
        "region_index": indices,
        "chunk_manifest": chunk_manifest,
        "source": {
            "region_view": str(rv_path),
            "region_view_sha256": actual_rv_sha,
            "machine_summary": str(summary_path),
            "machine_summary_sha256": summary_sha,
            "backbone_comparison": str(comparison_path),
            "backbone_comparison_sha256": comparison_sha,
            "layered_reconstruction": sample_record["paths"]["layered_reconstruction"],
            "layered_reconstruction_sha256": sample_record["sha256"][
                "layered_reconstruction"
            ],
        },
        "built_at": datetime.now().astimezone().isoformat(timespec="minutes"),
    }

    sample_options = "".join(
        '<option value="{}.html"{}>{}</option>'.format(
            html.escape(name),
            " selected" if name == sample else "",
            html.escape(name),
        )
        for name in DATASET_ORDER
    )
    run_id = Path(canonical["run_root"]).name

    template = r"""<!doctype html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light">
<link rel="icon" href="data:,">
<meta name="intersubmod-current-summary-sha256" content="__SUMMARY_SHA__">
<meta name="intersubmod-region-view-sha256" content="__REGION_SHA__">
<meta name="intersubmod-backbone-comparison-sha256" content="__COMPARISON_SHA__">
<meta name="intersubmod-canonical-sample" content="__SAMPLE__">
<title>全基因分層拓撲工作站 · __SAMPLE__</title>
<style>
:root{
 color-scheme:light;--canvas:#edf2f6;--paper:#fff;--ink:#152231;--muted:#536476;--line:#ced8e2;
 --blue:#075ea8;--cyan:#087d8f;--magenta:#a9236d;--green:#167452;--amber:#9a5200;--violet:#7040a0;
 --red:#a83227;--slate:#5c6875;--soft:#f5f8fb;--shadow:0 14px 36px rgba(20,39,58,.08);
 --mono:ui-monospace,"SFMono-Regular","Cascadia Code","Roboto Mono",monospace;
 --sans:"Aptos","Noto Sans TC","Microsoft JhengHei",system-ui,sans-serif
}
*{box-sizing:border-box}
html{background:var(--canvas);scroll-behavior:smooth}
body{margin:0;color:var(--ink);background:linear-gradient(180deg,#e6eef4 0,#f5f8fa 420px,var(--canvas) 100%);font:14px/1.58 var(--sans)}
a{color:var(--blue);text-underline-offset:3px}
button,select,input{font:inherit}
button,select,input,summary,a{touch-action:manipulation}
:focus-visible{outline:3px solid #ef9f20;outline-offset:3px}
.skip-link{position:fixed;z-index:100;left:14px;top:10px;transform:translateY(-180%);padding:9px 13px;border:2px solid #111;background:#fff;color:#111}
.skip-link:focus{transform:none}
.wrap{width:min(100%,1320px);margin:0 auto;padding:22px clamp(13px,2.7vw,36px) 48px}
.hero{overflow:hidden;border:1px solid #c4d0dc;border-top:5px solid var(--blue);background:rgba(255,255,255,.97);box-shadow:var(--shadow)}
.hero-main{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:24px;padding:clamp(22px,3vw,38px)}
.crumbs{display:flex;gap:9px;align-items:center;margin-bottom:13px;color:var(--muted);font:700 11px/1.2 var(--mono);letter-spacing:.05em;text-transform:uppercase}
.crumbs a{text-decoration:none}.crumbs a:hover{text-decoration:underline}
.eyebrow{display:flex;gap:9px;align-items:center;margin:0 0 12px;color:#365064;font:750 10.5px/1.2 var(--mono);letter-spacing:.11em;text-transform:uppercase}
.signal{width:9px;height:9px;border-radius:50%;background:var(--green);box-shadow:0 0 0 4px rgba(22,116,82,.13)}
h1{margin:0;font-size:clamp(27px,4vw,48px);font-weight:690;line-height:1.06;letter-spacing:-.035em}
.lede{max-width:850px;margin:14px 0 0;color:#3e5162;font-size:clamp(14px,1.5vw,17px)}
.hero-actions{display:grid;align-content:start;gap:10px;min-width:230px}
.control-label{display:grid;gap:5px;color:var(--muted);font-size:11px;font-weight:750}
select,input,.button{min-height:44px;border:1px solid #bcc9d6;border-radius:4px;background:#fff;color:var(--ink);padding:9px 11px}
.button{cursor:pointer;font-weight:750}.button:hover{border-color:var(--blue);color:var(--blue)}
.button.primary{border-color:var(--blue);background:var(--blue);color:#fff}.button.primary:hover{background:#044d8b;color:#fff}
.scope-ribbon{display:grid;grid-template-columns:1.2fr .9fr 1fr 1fr;border-top:1px solid var(--line);background:#f7fafc}
.scope-item{min-width:0;padding:12px 16px;border-right:1px solid var(--line)}
.scope-item:last-child{border-right:0}
.scope-item span{display:block;color:var(--muted);font:750 9.5px/1.2 var(--mono);letter-spacing:.08em;text-transform:uppercase}
.scope-item strong{display:block;margin-top:5px;overflow-wrap:anywhere;font-size:12px}
.scope-item.scope strong{color:var(--blue);font-size:14px}
.claim-row{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px;margin-top:14px}
.claim{padding:14px 16px;border:1px solid var(--line);background:#fff}
.claim:nth-child(1){border-top:3px solid var(--blue)}.claim:nth-child(2){border-top:3px solid var(--cyan)}.claim:nth-child(3){border-top:3px solid var(--amber)}
.claim b{display:block;margin-bottom:4px;color:var(--blue);font-size:11px}.claim p{margin:0;color:#45586a;font-size:12.5px}
.sensitivity-banner{display:grid;grid-template-columns:minmax(220px,.75fr) minmax(0,1.8fr);gap:18px;margin-top:12px;padding:14px 16px;border:1px solid #dfb27c;border-left:6px solid var(--amber);background:#fff8ee}
.sensitivity-banner h2{font-size:17px}.sensitivity-banner p{margin:5px 0 0;color:#624522;font-size:12px}.sensitivity-metrics{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:6px}.sensitivity-metric{padding:8px 9px;border:1px solid #ead1b2;background:#fff}.sensitivity-metric span{display:block;color:#7a5b35;font-size:9px}.sensitivity-metric strong{display:block;margin-top:4px;font:750 14px/1.1 var(--mono)}
.section{margin-top:26px}
.section-heading{display:flex;align-items:end;justify-content:space-between;gap:18px;margin-bottom:10px}
.kicker{margin:0 0 3px;color:var(--blue);font:750 10px/1.2 var(--mono);letter-spacing:.1em;text-transform:uppercase}
h2{margin:0;font-size:clamp(20px,2.4vw,28px);line-height:1.2;letter-spacing:-.02em}
h3{margin:0;font-size:17px;line-height:1.3}
.section-note{max-width:660px;margin:0;color:var(--muted);font-size:12.5px}
.metrics{display:grid;grid-template-columns:repeat(7,minmax(0,1fr));border:1px solid #c8d3de;background:#fff;box-shadow:var(--shadow)}
.metric{min-width:0;padding:14px;border-right:1px solid var(--line)}.metric:last-child{border-right:0}
.metric .value{display:block;font:720 clamp(17px,2vw,24px)/1.05 var(--mono);font-variant-numeric:tabular-nums}
.metric .label{display:block;margin-top:6px;color:var(--muted);font-size:11px}
.topology-strip{display:flex;overflow:hidden;height:12px;border:1px solid #b9c7d4;background:#e7edf2}
.seg-exact{background:var(--green)}.seg-shape{background:var(--blue)}.seg-multiple{background:var(--magenta)}.seg-incomplete{background:repeating-linear-gradient(135deg,#7a8793 0 5px,#a9b3bc 5px 10px)}
.topology-legend{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:8px;margin-top:8px}
.legend-card{padding:10px 12px;border:1px solid var(--line);background:#fff}
.legend-card b{display:block;font:750 13px/1.2 var(--mono)}.legend-card small{display:block;margin-top:3px;color:var(--muted)}
.genome-panel{border:1px solid #c7d2dd;background:#fff;box-shadow:var(--shadow)}
.genome-toolbar{display:flex;align-items:center;justify-content:space-between;gap:14px;padding:12px 14px;border-bottom:1px solid var(--line);background:#f8fafc}
.genome-toolbar p{margin:0;color:var(--muted);font-size:12px}
.chrom-grid{display:grid;grid-template-columns:repeat(11,minmax(0,1fr));gap:1px;background:var(--line)}
.chrom-button{min-height:86px;padding:10px;border:0;background:#fff;color:var(--ink);cursor:pointer;text-align:left}
.chrom-button:hover{background:#f2f7fb}.chrom-button.active{position:relative;background:#eaf4fc;box-shadow:inset 0 0 0 3px var(--blue)}
.chrom-button strong{display:flex;align-items:baseline;justify-content:space-between;gap:6px;font:750 12px/1.2 var(--mono)}
.chrom-button strong span{color:var(--muted);font-size:10px}
.mini-bar{display:flex;overflow:hidden;height:6px;margin:9px 0 6px;background:#e4eaf0}
.chrom-button small{display:block;color:var(--muted);font-size:9.5px;line-height:1.35}
.dimension-guide{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px}
.dimension-card{display:flex;min-width:0;min-height:315px;flex-direction:column;padding:18px;border:1px solid var(--line);border-top:4px solid var(--blue);background:#fff;box-shadow:var(--shadow)}
.dimension-card:nth-child(2){border-top-color:var(--magenta)}.dimension-card:nth-child(3){border-top-color:var(--cyan)}
.dimension-label{display:flex;align-items:center;justify-content:space-between;gap:10px;color:var(--blue);font:750 10px/1.2 var(--mono);letter-spacing:.08em;text-transform:uppercase}
.dimension-card:nth-child(2) .dimension-label{color:var(--magenta)}.dimension-card:nth-child(3) .dimension-label{color:var(--cyan)}
.dimension-card h3{margin:14px 0 7px;font-size:18px}.dimension-answer{margin:0;color:#30485c;font-size:12.5px}
.dimension-answer strong{color:var(--ink);font:750 21px/1.1 var(--mono)}
.dimension-stats{display:grid;gap:5px;margin:14px 0}.dimension-stat{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:10px;padding:7px 8px;border:1px solid #dae2e9;background:#f8fafc}
.dimension-stat span{min-width:0;color:var(--muted);font-size:10.5px}.dimension-stat b{font:750 11px/1.25 var(--mono);text-align:right}
.dimension-rule{margin:0 0 14px;color:var(--muted);font-size:11.5px}.dimension-card .button{width:100%;margin-top:auto;text-align:center}
.dimension-glossary{margin-top:10px}.glossary-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:8px}
.glossary-item{padding:10px;border:1px solid var(--line);background:#f8fafc}.glossary-item b{display:block;color:var(--blue);font-size:11px}.glossary-item p{margin:4px 0 0!important}
.glossary-boundary{margin-top:10px!important;padding:9px 10px;border-left:4px solid var(--amber);background:#fff8ee;color:#624522!important}
.browser-shell{border-top:4px solid var(--blue)}
.filters{display:grid;grid-template-columns:1.15fr 1fr 1fr 1fr 1.35fr auto;gap:9px;align-items:end;padding:14px;border:1px solid var(--line);background:#f8fafc}
.filter{display:grid;gap:4px;min-width:0}.filter span{color:var(--muted);font-size:10px;font-weight:750;letter-spacing:.04em}
.filter select,.filter input{width:100%;min-width:0}
.filter-boundary{margin:0;padding:9px 13px;border:1px solid var(--line);border-top:0;background:#fff8ee;color:#624522;font-size:11px}
.count{display:flex;align-items:center;min-height:44px;padding:0 10px;border:1px solid var(--line);background:#fff;color:var(--muted);font:700 11px/1.2 var(--mono);white-space:nowrap}
.workspace{display:grid;grid-template-columns:minmax(320px,390px) minmax(0,1fr);gap:12px;align-items:start;margin-top:12px}
.results{min-width:0;border:1px solid var(--line);background:#fff}
.results-head{display:flex;justify-content:space-between;gap:8px;padding:10px 12px;border-bottom:1px solid var(--line);background:#f8fafc;font-size:11.5px}
.result-list{max-height:72vh;overflow-y:auto}
.result-row{display:block;width:100%;min-height:70px;padding:10px 11px;border:0;border-bottom:1px solid #e2e8ee;background:#fff;color:var(--ink);cursor:pointer;text-align:left}
.result-row:hover{background:#f5f9fc}.result-row[aria-current="true"]{background:#e9f3fb;box-shadow:inset 4px 0 0 var(--blue)}
.result-title{display:flex;align-items:center;justify-content:space-between;gap:8px;font:730 11.5px/1.35 var(--mono)}
.result-dimensions{display:flex;flex-wrap:wrap;gap:5px;margin-top:6px;color:var(--muted);font-size:11px}.result-dimensions span{padding-right:6px;border-right:1px solid var(--line)}.result-dimensions span:last-child{padding-right:0;border-right:0}
.badges{display:flex;flex-wrap:wrap;gap:4px;margin-top:6px}
.badge{display:inline-flex;align-items:center;min-height:22px;padding:2px 7px;border:1px solid #c8d3dd;border-radius:999px;background:#f8fafc;color:#48596a;font-size:10px;line-height:1.2}
.badge.exact{border-color:#9acbb9;background:#eef8f4;color:#0d6246}.badge.shape{border-color:#95bce0;background:#edf5fc;color:#0a548f}
.badge.multiple{border-color:#d5a3c2;background:#fbf0f7;color:#842057}.badge.incomplete{border-style:dashed;border-color:#8995a0;background:#f3f5f7;color:#44515d}
.badge.none{border-color:#c6cdd4;background:#f3f5f6;color:#56616b}.badge.warn{border-color:#e1b77c;background:#fff7eb;color:#764000}
.load-more{display:block;width:calc(100% - 20px);margin:10px;min-height:44px}
[hidden]{display:none!important}
.detail{min-width:0;border:1px solid var(--line);background:#fff;box-shadow:var(--shadow)}
.detail-empty{min-height:280px;padding:22px;color:var(--muted)}
.detail-toolbar{position:sticky;z-index:10;top:0;display:flex;align-items:center;justify-content:space-between;gap:8px;padding:9px 12px;border-bottom:1px solid var(--line);background:rgba(255,255,255,.97);backdrop-filter:blur(8px)}
.detail-toolbar .routes{display:flex;gap:7px;flex-wrap:wrap}.detail-toolbar button{min-height:40px}
.detail-body{padding:16px}
.verdict{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:16px;padding:15px 16px;border:1px solid #bdd0e0;background:#f2f8fc}
.verdict.incomplete{border-style:dashed;background:repeating-linear-gradient(135deg,#fafbfc 0 12px,#f1f4f6 12px 24px)}
.verdict h3{font-size:19px}.verdict p{margin:6px 0 0;color:#46596b}
.ct-readout{display:grid;grid-template-columns:1fr 1fr;gap:6px;align-content:start}
.ct-box{min-width:88px;padding:9px;border:1px solid var(--line);background:#fff;text-align:center}.ct-box span{display:block;color:var(--muted);font:700 9px/1 var(--mono)}.ct-box strong{display:block;margin-top:5px;font:750 20px/1 var(--mono)}
.facet-row{display:flex;flex-wrap:wrap;gap:6px;margin-top:9px}
.assertion-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:8px;margin-top:10px}
.assertion{padding:11px 12px;border:1px solid var(--line);background:#fff}.assertion b{display:block;color:var(--blue);font-size:11px}.assertion p{margin:4px 0 0;color:var(--muted);font-size:12px}
.region-dimensions{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:8px;margin-top:10px}
.region-dimension{min-width:0;padding:11px 12px;border:1px solid var(--line);background:#fff}.region-dimension:nth-child(1){border-top:3px solid var(--blue)}.region-dimension:nth-child(2){border-top:3px solid var(--magenta)}.region-dimension:nth-child(3){border-top:3px solid var(--cyan)}
.region-dimension span{display:block;color:var(--muted);font:750 9.5px/1.2 var(--mono);letter-spacing:.06em;text-transform:uppercase}.region-dimension strong{display:block;margin-top:6px;overflow-wrap:anywhere;font-size:14px}.region-dimension p{margin:5px 0 0;color:var(--muted);font-size:11px}
.subsection{margin-top:18px}.subsection-head{display:flex;align-items:end;justify-content:space-between;gap:12px;margin-bottom:8px}.subsection h4{margin:0;font-size:15px}.subsection-note{color:var(--muted);font-size:11.5px}
.scroll-cue{margin:0 0 5px;color:var(--muted);font-size:10.5px}
.scroll-region{width:100%;max-width:100%;overflow-x:auto;overscroll-behavior-inline:contain;-webkit-overflow-scrolling:touch}
.scroll-region:focus-visible{outline-offset:-3px}
table{width:100%;border-collapse:collapse;background:#fff;font-size:11.5px}
caption{padding:6px 0;color:var(--muted);text-align:left;font-weight:700}
th,td{padding:7px 8px;border:1px solid var(--line);vertical-align:top}th{background:#f2f6f9;color:#3d5061;text-align:left}td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}
.site-table{min-width:680px}.read-table{min-width:520px}
.site-state{font:750 10px/1.2 var(--mono);text-transform:uppercase}.site-state.observed{color:var(--blue)}.site-state.absent{color:var(--red)}
.lane{margin-top:12px;border:1px solid var(--line);border-left:5px solid var(--blue);background:#fff}
.lane.family-2{border-left-color:var(--magenta)}.lane.auxiliary{border-left-style:dashed;border-left-color:#7c8995;background:#fafbfc}
.lane-head{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:10px;padding:12px 13px;border-bottom:1px solid var(--line);background:#f8fafc}
.lane-title{display:flex;flex-wrap:wrap;align-items:center;gap:6px}.lane-title b{font-size:14px}
.lane-body{padding:12px}
.candidate-summary{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:7px;margin-bottom:10px}
.candidate-stat{padding:8px;border:1px solid var(--line);background:#fff}.candidate-stat span{display:block;color:var(--muted);font-size:9.5px}.candidate-stat strong{display:block;margin-top:4px;font:730 15px/1.1 var(--mono)}
.scope-warning{padding:9px 10px;border:1px dashed #8b96a0;background:repeating-linear-gradient(135deg,#fafbfc 0 10px,#f1f4f6 10px 20px);color:#455464;font-size:11.5px}
.shape-scope,.recurrence-note{display:grid;gap:3px;margin-top:8px;padding:9px 10px;border:1px solid var(--line);background:#fff;color:#455464;font-size:11px}.shape-scope b,.recurrence-note b{color:var(--blue)}
.recurrence-note{border-color:#d6a3c2;background:#fff6fb}.recurrence-note b{color:var(--magenta)}
.network-card{margin-top:10px;border:1px solid var(--line);background:#f7f9fb}
.network-head{display:flex;align-items:flex-start;justify-content:space-between;gap:10px;padding:10px 11px;border-bottom:1px solid var(--line)}
.network-head p{margin:3px 0 0;color:var(--muted);font-size:11px}
.network-scroll{width:100%;max-width:100%;overflow-x:auto;padding:10px;overscroll-behavior-inline:contain}
.network-scroll svg{display:block;max-width:none;height:auto;margin:0 auto}
.network-legend{display:flex;flex-wrap:wrap;gap:8px;padding:8px 10px;border-top:1px solid var(--line);background:#fff;color:var(--muted);font-size:10.5px}
.line-swatch{display:inline-block;width:24px;margin-right:4px;border-top:2px solid #536476;vertical-align:middle}.line-swatch.variable{border-top-style:dashed;border-color:var(--amber)}.line-swatch.unevaluated{border-top-style:dotted;border-color:#8b97a2}
.shape-tabs,.tree-tabs{display:flex;flex-wrap:wrap;gap:5px;margin-top:8px}
.shape-tab,.tree-tab{min-height:38px;padding:6px 9px;border:1px solid var(--line);background:#fff;color:var(--muted);cursor:pointer}
.shape-tab[aria-pressed="true"],.tree-tab[aria-pressed="true"]{border-color:var(--blue);background:var(--blue);color:#fff}
.candidate-view{margin-top:8px}
.tree-caption{padding:7px 10px;color:var(--muted);font-size:11px}
.sidecars{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:8px}
.sidecar{padding:12px;border:1px solid var(--line);background:#f8fafc}.sidecar b{display:block}.sidecar p{margin:5px 0 0;color:var(--muted);font-size:11.5px}
details.drawer{margin-top:10px;border:1px solid var(--line);background:#fff}
details.drawer>summary{cursor:pointer;padding:11px 12px;font-weight:730}
details.drawer[open]>summary{border-bottom:1px solid var(--line)}
.drawer-body{padding:11px 12px}.drawer-body p{margin:5px 0;color:var(--muted);font-size:11.5px}
.trace{padding:8px 10px;background:#f3f6f8;color:#4f5f6d;font:10.5px/1.5 var(--mono);overflow-wrap:anywhere}
.raw-links{display:grid;gap:8px}.raw-link{padding:9px 10px;border:1px solid var(--line);overflow-wrap:anywhere}.raw-link small{display:block;color:var(--muted)}
.footer{margin-top:22px;padding-top:12px;border-top:1px solid var(--line);color:var(--muted);font-size:11.5px}
.mono{font-family:var(--mono);font-variant-numeric:tabular-nums}
.muted{color:var(--muted)}
.sr-only{position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0}
.error{padding:15px;border:2px solid var(--red);background:#fff0ee;color:#6e1d17}
@media(max-width:1080px){
 .metrics{grid-template-columns:repeat(4,minmax(0,1fr))}.metric:nth-child(4){border-right:0}.metric:nth-child(-n+4){border-bottom:1px solid var(--line)}
 .chrom-grid{grid-template-columns:repeat(6,minmax(0,1fr))}.filters{grid-template-columns:repeat(3,minmax(0,1fr))}.workspace{grid-template-columns:340px minmax(0,1fr)}
}
@media(max-width:860px){
 .wrap{padding:14px 12px 36px}.hero-main{grid-template-columns:1fr}.hero-actions{min-width:0;grid-template-columns:1fr 1fr}
 .scope-ribbon{grid-template-columns:1fr 1fr}.scope-item:nth-child(2){border-right:0}.scope-item:nth-child(-n+2){border-bottom:1px solid var(--line)}
 .claim-row,.dimension-guide,.assertion-grid,.region-dimensions,.sidecars,.sensitivity-banner{grid-template-columns:1fr}.sensitivity-metrics{grid-template-columns:1fr 1fr}
 .glossary-grid{grid-template-columns:1fr 1fr}
 .metrics{grid-template-columns:repeat(2,minmax(0,1fr))}.metric{border-bottom:1px solid var(--line)}.metric:nth-child(even){border-right:0}.metric:last-child{border-bottom:0}
 .topology-legend{grid-template-columns:1fr 1fr}.chrom-grid{grid-template-columns:repeat(4,minmax(0,1fr))}
 .workspace{grid-template-columns:1fr}.result-list{max-height:44vh}.detail{scroll-margin-top:12px}.filters{grid-template-columns:1fr 1fr}
}
@media(max-width:520px){
 body{font-size:13.5px}.wrap{padding:10px 9px 30px}.hero-main{padding:20px 16px}.hero-actions{grid-template-columns:1fr}
 .scope-ribbon{grid-template-columns:1fr}.scope-item{border-right:0;border-bottom:1px solid var(--line)}.scope-item:last-child{border-bottom:0}
 .section-heading{display:block}.section-note{margin-top:5px}
 .metrics{grid-template-columns:1fr 1fr}.metric .value{font-size:17px}.topology-legend{grid-template-columns:1fr}
 .genome-toolbar{display:block}.genome-toolbar .button{width:100%;margin-top:8px}.chrom-grid{grid-template-columns:1fr 1fr}
 .filters{grid-template-columns:1fr 1fr}.filter.search{grid-column:1/-1}.count{grid-column:1/-1}
 .result-list{max-height:38vh}.detail-body{padding:10px}.detail-toolbar{align-items:flex-start}.detail-toolbar .routes{gap:5px}
 .verdict{grid-template-columns:1fr}.ct-readout{grid-template-columns:1fr 1fr}.candidate-summary{grid-template-columns:1fr 1fr}
 .network-scroll{padding:8px}.site-table{min-width:660px}.scroll-cue{font-weight:650}
 .dimension-card{min-height:0}.glossary-grid{grid-template-columns:1fr}
}
@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}}
</style>
</head>
<body>
<a class="skip-link" href="#genome-overview">跳至 chr1–22 全基因視圖</a>
<main class="wrap">
 <header class="hero">
  <div class="hero-main">
   <div>
    <nav class="crumbs" aria-label="麵包屑"><a href="index.html">Cohort console</a><span aria-hidden="true">/</span><span>Dataset __SAMPLE__</span></nav>
    <p class="eyebrow"><span class="signal" aria-hidden="true"></span>Canonical v5 · verifier pass · offline workstation</p>
    <h1>全基因分層拓撲工作站<br><span class="mono">__SAMPLE__</span></h1>
    <p class="lede">從 chr1–22 全景縮放到 region，再檢視 HP1／HP2 主重建單位的 exact candidate 與 topology shape。所有結果都是區域 mutation-state 候選樹，不是已確認的細胞 clone 或祖先關係。</p>
   </div>
   <div class="hero-actions">
    <label class="control-label" for="dataset-switch">切換 Dataset<select id="dataset-switch" onchange="location.href=this.value">__SAMPLE_OPTIONS__</select></label>
    <button class="button" type="button" id="copy-link">複製目前檢視連結</button>
   </div>
  </div>
  <div class="scope-ribbon" aria-label="資料範圍與權威來源">
   <div class="scope-item scope"><span>Whole-genome scope</span><strong>chr1–22 · W_primary 主分析</strong></div>
   <div class="scope-item"><span>Cohort / current dataset</span><strong>7 datasets / 6 biological samples · __SAMPLE__ / __BIO_SAMPLE__</strong></div>
   <div class="scope-item"><span>Backbone</span><strong>LongPhase-S recalibrated · FILTER=PASS</strong></div>
   <div class="scope-item"><span>Canonical run</span><strong class="mono">__RUN_ID__</strong></div>
  </div>
 </header>

 <div class="claim-row" aria-label="工作站解讀邊界">
  <article class="claim"><b>回答什麼</b><p>目前 read-state 約束下，哪些區域候選集合完整，以及 exact candidates 與 topology shapes 各有多少。</p></article>
  <article class="claim"><b>證據如何分層</b><p>Observed read states 與 inferred hidden states 分開；HP1／HP2 是主線，H3／H4／none／reference control 收在輔助層。</p></article>
  <article class="claim"><b>不能宣稱什麼</b><p>read ALT fraction 不等於 CCF；CN、PS、L3 不排序候選；partial-only 代表 overlap-constrained，不代表沒有訊號。</p></article>
 </div>

 <aside class="sensitivity-banner" aria-labelledby="sensitivity-title">
  <div><p class="kicker">Dataset backbone sensitivity · cohort __BACKBONE_VERDICT__</p><h2 id="sensitivity-title">__SAMPLE_BACKBONE_VERDICT__</h2><p>Canonical v5 仍是唯一主分母；ClairS FILTER=PASS 只作 sensitivity，不可混入本頁 C／Topo。</p></div>
  <div class="sensitivity-metrics" aria-label="Dataset backbone sensitivity comparison metrics">
   <div class="sensitivity-metric"><span>Dataset retained-position Jaccard</span><strong>__RETAINED_JACCARD__</strong></div>
   <div class="sensitivity-metric"><span>Dataset primary-unit-key Jaccard</span><strong>__PRIMARY_JACCARD__</strong></div>
   <div class="sensitivity-metric"><span>Dataset shared topology digest</span><strong>__TOPOLOGY_CONCORDANCE__</strong></div>
   <div class="sensitivity-metric"><span>Max absolute proportion delta</span><strong>__MAX_DELTA_PP__ pp</strong></div>
  </div>
 </aside>

 <section class="section" aria-labelledby="sample-funnel-title">
  <div class="section-heading">
   <div><p class="kicker">Canonical disposition</p><h2 id="sample-funnel-title">Dataset 全範圍漏斗</h2></div>
   <p class="section-note">各數字直接取自 canonical v5 machine summary；W_tree 含無主線區，W_primary 僅含 mutation-bearing HP1／HP2。</p>
  </div>
  <div class="metrics" id="sample-metrics"></div>
 <div class="topology-strip" id="topology-strip" aria-hidden="true"></div>
 <div class="topology-legend" id="topology-legend"></div>
 </section>

 <section class="section" id="dimension-guide" aria-labelledby="dimension-title">
  <div class="section-heading">
   <div><p class="kicker">Three reading dimensions</p><h2 id="dimension-title">先回答三個問題，再看候選網路</h2></div>
   <p class="section-note">每張卡都把名詞、目前 dataset 觀察、下鑽入口與推論邊界放在一起；舊 clone-first 類別與 archive 數字不沿用。</p>
  </div>
  <div class="dimension-guide" id="dimension-cards"></div>
  <details class="drawer dimension-glossary">
   <summary>名詞與概念解釋：canonical 欄位如何對應三個問題</summary>
   <div class="drawer-body">
    <div class="glossary-grid">
     <article class="glossary-item"><b>C · exact candidate 組合數</b><p>各 mutation-bearing HP1／HP2 primary unit 的 <span class="mono">n_trees</span> 乘積；C=1 才是 exact candidate 唯一。</p></article>
     <article class="glossary-item"><b>Topo · 無標籤形狀組合數</b><p>各 primary unit 的 <span class="mono">n_distinct_shapes_exact</span> 乘積；Topo=1 只表示 shape 唯一。</p></article>
     <article class="glossary-item"><b>可辨識度 · determinacy</b><p>依 candidate-complete 與 C／Topo 分成 exact 唯一、shape 唯一、multi-shape、未評估、N/A 五層。</p></article>
     <article class="glossary-item"><b>基因體位置 · coordinate</b><p>顯示 canonical sSNV 的 <span class="mono">chr:start–end</span>、端點距離 <span class="mono">end-start</span> 與位點數。</p></article>
    </div>
    <p class="glossary-boundary"><b>固定邊界：</b>Topo=1 不是 biological truth、confirmed clone 或唯一時間順序；Incomplete 是「尚不可評估」而不是「不唯一」；位置分布只作描述與巡覽，不代表 hotspot 或 enrichment。</p>
   </div>
  </details>
 </section>

 <section class="section" id="genome-overview" aria-labelledby="genome-title">
  <div class="section-heading">
   <div><p class="kicker">Genome → chromosome</p><h2 id="genome-title">chr1–22 全基因視圖</h2></div>
   <p class="section-note">22 條染色體永遠保留在 scope；點擊只改變下方 active view，不會改寫全基因分母。</p>
  </div>
  <div class="genome-panel">
   <div class="genome-toolbar"><p id="genome-status">全 chr1–22 · 正在建立 topology composition</p><button type="button" class="button" id="all-genome">顯示全部 chr1–22</button></div>
   <div class="chrom-grid" id="chrom-grid" aria-label="22 條常染色體導覽"></div>
  </div>
 </section>

 <section class="section browser-shell" id="region-browser" aria-labelledby="browser-title">
  <div class="section-heading">
   <div><p class="kicker">Chromosome → region → HP lane</p><h2 id="browser-title">Region 候選集合檢視</h2></div>
   <p class="section-note">先用獨立 facets 找區域，再以 assertion-first 順序閱讀 verdict、evidence、candidate network 與 sidecars。</p>
  </div>
  <div class="filters" aria-label="Region 篩選">
   <label class="filter" for="fchr"><span>染色體</span><select id="fchr"><option value="">全部 chr1–22</option></select></label>
   <label class="filter" for="ftopo"><span>可辨識度（C / Topo）</span><select id="ftopo"><option value="">全部層次</option><option value="exact_and_topology_unique">Exact 唯一 · C=1 / Topo=1</option><option value="topology_unique_exact_multiple">Shape 唯一 · C&gt;1 / Topo=1</option><option value="topology_multiple_exact_multiple">Multi-shape · C&gt;1 / Topo&gt;1</option><option value="incomplete">未評估 · Incomplete</option><option value="no_primary_lineage">不適用 · 無 primary</option></select></label>
   <label class="filter" for="fevidence"><span>Read-state 約束</span><select id="fevidence"><option value="">全部</option><option value="full_and_partial">full + partial</option><option value="partial_only">partial-only</option><option value="full_only">full-only</option></select></label>
   <label class="filter" for="fsignal"><span>獨立 facet</span><select id="fsignal"><option value="">全部</option><option value="recurrence">有 recurrence</option><option value="hidden">有 hidden state</option><option value="multi_hp">HP1 + HP2</option></select></label>
   <label class="filter search" for="fq"><span>基因體位置／區間 overlap</span><input id="fq" type="search" placeholder="chr8:34220481 或 chr8:34200000-34300000"></label>
   <div class="count" id="result-count" role="status" aria-live="polite">0 regions</div>
  </div>
  <p class="filter-boundary"><b>PS QC boundary：</b>此 dataset 有 __MIXED_PS__ 個 mixed-PS regions，但 canonical region payload 沒有逐 region PS membership；因此本頁刻意不提供 PS filter，且不產生任何 PS-derived edge。</p>
  <div class="workspace">
   <aside class="results" id="results-panel" aria-label="Region 結果">
    <div class="results-head"><strong>符合條件的 regions</strong><span>每次顯示 80 筆</span></div>
    <div class="result-list" id="result-list"></div>
    <button class="button load-more" id="load-more" type="button">顯示更多 regions</button>
   </aside>
   <article class="detail" id="detail" tabindex="-1" aria-labelledby="detail-title">
    <div class="detail-empty"><h3 id="detail-title">尚未選擇 region</h3><p>從左側結果或上方 chromosome overview 選取；detail 會先回答拓樸型態、可辨識度與基因體位置，再展開 read evidence。</p></div>
   </article>
   <p class="sr-only" id="detail-status" role="status" aria-live="polite"></p>
  </div>
 </section>

 <details class="drawer section evidence-drawer">
  <summary>方法、驗證雜湊與原始資料連結</summary>
  <div class="drawer-body raw-links" id="source-links"></div>
 </details>
 <footer class="footer">Offline standalone artifact · 7 datasets / 6 biological samples · chr1–22 · canonical v5 · PS only for QC · L3 not evaluated and does not rank candidates.</footer>
</main>
<script type="application/json" id="workstation-data">__PAGE_DATA__</script>
__CHUNK_SCRIPTS__
<script>
"use strict";
const DATA=JSON.parse(document.getElementById("workstation-data").textContent);
const INDEX=DATA.region_index;
const SAMPLE=DATA.sample;
const SAMPLE_REC=DATA.canonical_sample;
const CHUNK_CACHE={};
const CANDIDATE_CACHE={};
const state={filtered:[],limit:80,selectedRegion:null,selectedRow:null};
const FILTER_PARAMS={fchr:"chr",ftopo:"topo",fevidence:"evidence",fsignal:"signal",fq:"q"};
let networkCounter=0;

const esc=value=>String(value==null?"":value).replace(/[&<>"']/g,char=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[char]));
const fmt=value=>Number(value||0).toLocaleString();
const pct=(value,total)=>total?(100*Number(value||0)/Number(total)).toFixed(1)+"%":"—";
const fmtSpan=value=>{
 const span=Number(value||0);
 if(span>=1000000)return (span/1000000).toFixed(span>=10000000?1:2)+" Mb";
 if(span>=1000)return (span/1000).toFixed(span>=100000?0:1)+" kb";
 return fmt(span)+" bp";
};
const chromRank=chrom=>Number(String(chrom).replace("chr",""))||99;
const classInfo={
 exact_and_topology_unique:{short:"C=1 · Topo=1",label:"exact 與 topology 皆唯一",badge:"exact"},
 topology_unique_exact_multiple:{short:"C>1 · Topo=1",label:"topology 唯一；exact 不唯一",badge:"shape"},
 topology_multiple_exact_multiple:{short:"C>1 · Topo>1",label:"topology 與 exact 皆不唯一",badge:"multiple"},
 incomplete:{short:"Incomplete",label:"候選集合未完整",badge:"incomplete"},
 no_primary_lineage:{short:"No primary",label:"無 mutation-bearing HP1／HP2",badge:"none"}
};

function identifiabilityInfo(row){
 const mapping={
  exact_and_topology_unique:{value:"Exact 唯一",note:"C=1；exact candidate 與 shape 都唯一"},
  topology_unique_exact_multiple:{value:"Shape 唯一",note:"C>1、Topo=1；只辨識到無標籤形狀"},
  topology_multiple_exact_multiple:{value:"Multi-shape",note:"C>1、Topo>1；仍有多種相容形狀"},
  incomplete:{value:"尚未評估",note:"candidate set 未完整；C／Topo unavailable"},
  no_primary_lineage:{value:"不適用",note:"沒有 mutation-bearing HP1／HP2 primary unit"}
 };
 return mapping[row.topology_class]||mapping.no_primary_lineage;
}

function topologyShapeText(row){
 if(row.topology_class==="no_primary_lineage")return "不適用";
 if(row.topology_class==="incomplete")return "尚未評估";
 return fmt(row.Topo)+" 種 shape";
}

function topologyClassBadge(key){
 const info=classInfo[key]||classInfo.no_primary_lineage;
 return '<span class="badge '+info.badge+'">'+info.short+'</span>';
}

function loadChrom(chrom){
 if(CHUNK_CACHE[chrom])return CHUNK_CACHE[chrom];
 const node=document.getElementById("chunk-"+chrom);
 if(!node)throw new Error("Missing embedded chromosome chunk: "+chrom);
 try{
  const rows=JSON.parse(node.textContent);
  if(rows.length!==DATA.chunk_manifest[chrom].count)throw new Error("count mismatch");
  CHUNK_CACHE[chrom]=rows;
  return rows;
 }catch(error){
  throw new Error("Cannot decode "+chrom+" detail chunk: "+error.message);
 }
}

function getDetail(indexRow){
 const rows=loadChrom(indexRow.chrom);
 const region=rows[indexRow.chunk_index];
 if(!region||region.region!==indexRow.region)throw new Error("Chunk/index identity mismatch for "+indexRow.region);
 return region;
}

function renderMetrics(){
 const t=SAMPLE_REC.topology_classes;
 const metrics=[
  [SAMPLE_REC.tree_input_records,"Tree inputs"],
  [SAMPLE_REC.autosomal_biallelic_sSNV,"Autosomal biallelic"],
  [SAMPLE_REC.retained_sSNV,"Retained sSNV"],
  [SAMPLE_REC.W_tree,"W_tree"],
  [SAMPLE_REC.W_primary,"W_primary"],
  [SAMPLE_REC.complete_regions,"Complete"],
  [SAMPLE_REC.incomplete_regions,"Incomplete"]
 ];
 document.getElementById("sample-metrics").innerHTML=metrics.map(item=>'<div class="metric"><span class="value">'+fmt(item[0])+'</span><span class="label">'+item[1]+'</span></div>').join("");
 const parts=[
  ["exact_and_topology_unique","seg-exact"],
  ["topology_unique_exact_multiple","seg-shape"],
  ["topology_multiple_exact_multiple","seg-multiple"],
  ["incomplete","seg-incomplete"]
 ];
 document.getElementById("topology-strip").innerHTML=parts.map(item=>'<span class="'+item[1]+'" style="width:'+pct(t[item[0]],SAMPLE_REC.W_primary)+'" title="'+classInfo[item[0]].label+' '+fmt(t[item[0]])+'"></span>').join("");
 document.getElementById("topology-legend").innerHTML=parts.map(item=>'<div class="legend-card"><b>'+classInfo[item[0]].short+' · '+fmt(t[item[0]])+'</b><small>'+classInfo[item[0]].label+' · '+pct(t[item[0]],SAMPLE_REC.W_primary)+'</small></div>').join("");
}

function summarizeChromosomes(){
 const summaries={};
 for(let number=1;number<=22;number++){
  const chrom="chr"+number;
  summaries[chrom]={chrom:chrom,W_tree:0,W_primary:0,exact:0,shape:0,multiple:0,incomplete:0,noPrimary:0};
 }
 INDEX.forEach(row=>{
  const out=summaries[row.chrom];out.W_tree++;
  if(row.topology_class==="no_primary_lineage"){out.noPrimary++;return;}
  out.W_primary++;
  if(row.topology_class==="exact_and_topology_unique")out.exact++;
  else if(row.topology_class==="topology_unique_exact_multiple")out.shape++;
  else if(row.topology_class==="topology_multiple_exact_multiple")out.multiple++;
  else if(row.topology_class==="incomplete")out.incomplete++;
 });
 return summaries;
}

const CHROM_SUMMARIES=summarizeChromosomes();

function openDimensionView(options){
 const chrom=options.chrom||"";const topology=options.topology||"";
 document.getElementById("fevidence").value="";
 document.getElementById("fsignal").value="";
 document.getElementById("fq").value="";
 document.getElementById("fchr").value=chrom;
 document.getElementById("ftopo").value=topology;
 updateChromButtonState(chrom);
 document.getElementById("genome-status").textContent=chrom
  ? chrom+" active view · 卡片統計以 W_primary="+fmt(CHROM_SUMMARIES[chrom].W_primary)+"；下方列出該染色體全部 W_tree="+fmt(CHROM_SUMMARIES[chrom].W_tree)
  : "全部 chr1–22 · W_primary "+fmt(SAMPLE_REC.W_primary);
 applyFilters();renderEmptyDetail(false);syncHash("push");
 document.getElementById("region-browser").scrollIntoView({block:"start"});
 (topology?document.getElementById("ftopo"):document.getElementById("fchr")).focus({preventScroll:true});
}

function renderDimensions(){
 const t=SAMPLE_REC.topology_classes;
 const exact=Number(t.exact_and_topology_unique||0);
 const shapeOnly=Number(t.topology_unique_exact_multiple||0);
 const multiple=Number(t.topology_multiple_exact_multiple||0);
 const incomplete=Number(t.incomplete||0);
 const complete=exact+shapeOnly+multiple;
 const shapeUnique=exact+shapeOnly;
 const chromosomes=Object.values(CHROM_SUMMARIES);
 const countLeader=[...chromosomes].sort((a,b)=>(b.W_primary-a.W_primary)||(chromRank(a.chrom)-chromRank(b.chrom)))[0];
 const incompleteLeader=[...chromosomes].filter(row=>row.W_primary>0).sort((a,b)=>
  (b.incomplete/b.W_primary)-(a.incomplete/a.W_primary)||chromRank(a.chrom)-chromRank(b.chrom)
 )[0];
 document.getElementById("dimension-cards").innerHTML=
  '<article class="dimension-card"><div class="dimension-label"><span>01 · 拓樸型態</span><span>Topology shape</span></div><h3>分子累積形狀有幾種？</h3><p class="dimension-answer"><strong>'+fmt(shapeUnique)+'</strong> 個 complete regions 的 Topo=1。</p><div class="dimension-stats"><div class="dimension-stat"><span>形狀唯一 / complete</span><b>'+fmt(shapeUnique)+' / '+fmt(complete)+' · '+pct(shapeUnique,complete)+'</b></div><div class="dimension-stat"><span>多種形狀相容</span><b>'+fmt(multiple)+' · '+pct(multiple,complete)+'</b></div></div><p class="dimension-rule">Topo 比較無節點標籤的 regional mutation-state shape；不代表 biological ancestry 或真實時間順序。</p><button class="button dimension-action" type="button" data-dimension-topology="topology_multiple_exact_multiple">查看 multi-shape regions</button></article>'+
  '<article class="dimension-card"><div class="dimension-label"><span>02 · 可辨識度</span><span>Determinacy</span></div><h3>目前能唯一辨識到哪一層？</h3><p class="dimension-answer"><strong>'+fmt(exact)+'</strong> 個 regions 可辨識到 exact candidate。</p><div class="dimension-stats"><div class="dimension-stat"><span>Exact 唯一</span><b>'+fmt(exact)+' · '+pct(exact,SAMPLE_REC.W_primary)+'</b></div><div class="dimension-stat"><span>只到 shape 唯一</span><b>'+fmt(shapeOnly)+' · '+pct(shapeOnly,SAMPLE_REC.W_primary)+'</b></div><div class="dimension-stat"><span>Multi-shape</span><b>'+fmt(multiple)+' · '+pct(multiple,SAMPLE_REC.W_primary)+'</b></div><div class="dimension-stat"><span>尚未評估</span><b>'+fmt(incomplete)+' · '+pct(incomplete,SAMPLE_REC.W_primary)+'</b></div></div><p class="dimension-rule">Incomplete 是 candidate set 未完整，因此「未評估」；不是已知存在多解，也不能以 0 代替。</p><button class="button dimension-action" type="button" data-dimension-topology="exact_and_topology_unique">查看 exact-unique regions</button></article>'+
  '<article class="dimension-card"><div class="dimension-label"><span>03 · 基因體位置</span><span>Genome position</span></div><h3>這些 regions 位於哪裡？</h3><p class="dimension-answer"><strong>'+esc(countLeader.chrom)+'</strong> 的 W_primary 數量最多（'+fmt(countLeader.W_primary)+'）。</p><div class="dimension-stats"><div class="dimension-stat"><span>W_primary count leader</span><b>'+esc(countLeader.chrom)+' · '+fmt(countLeader.W_primary)+'</b></div><div class="dimension-stat"><span>Incomplete rate leader</span><b>'+esc(incompleteLeader.chrom)+' · '+fmt(incompleteLeader.incomplete)+' / '+fmt(incompleteLeader.W_primary)+' · '+pct(incompleteLeader.incomplete,incompleteLeader.W_primary)+'</b></div></div><p class="dimension-rule">卡片比較 W_primary；CTA 會列出該染色體全部 '+fmt(countLeader.W_tree)+' 個 W_tree regions。純描述且未校正 chromosome 長度、輸入 sSNV density 或 region construction，不作 hotspot／enrichment 判定。</p><button class="button dimension-action" type="button" data-dimension-chrom="'+esc(countLeader.chrom)+'">前往 '+esc(countLeader.chrom)+' · 全部 '+fmt(countLeader.W_tree)+' 個 W_tree regions</button></article>';
 document.querySelectorAll("[data-dimension-topology]").forEach(button=>button.addEventListener("click",()=>openDimensionView({topology:button.dataset.dimensionTopology})));
 document.querySelectorAll("[data-dimension-chrom]").forEach(button=>button.addEventListener("click",()=>openDimensionView({chrom:button.dataset.dimensionChrom})));
}

function updateChromButtonState(chrom){
 document.querySelectorAll(".chrom-button").forEach(button=>{
  const active=button.dataset.chrom===chrom;
  button.classList.toggle("active",active);
  button.setAttribute("aria-pressed",active?"true":"false");
 });
}

function setActiveChrom(chrom,moveToBrowser){
 document.getElementById("fchr").value=chrom||"";
 updateChromButtonState(chrom);
 document.getElementById("genome-status").textContent=chrom
  ? chrom+" active view · 全基因 denominator 仍為 chr1–22 / W_primary "+fmt(SAMPLE_REC.W_primary)
  : "全部 chr1–22 · W_primary "+fmt(SAMPLE_REC.W_primary)+" · "+fmt(SAMPLE_REC.complete_regions)+" complete / "+fmt(SAMPLE_REC.incomplete_regions)+" incomplete";
 applyFilters();
 if(moveToBrowser){syncHash("push");document.getElementById("region-browser").scrollIntoView({block:"start"})}
}

function renderGenome(){
 const host=document.getElementById("chrom-grid");
 host.innerHTML=Object.values(CHROM_SUMMARIES).map(row=>{
  const total=row.W_primary||1;
  const mini='<div class="mini-bar" aria-hidden="true"><span class="seg-exact" style="width:'+pct(row.exact,total)+'"></span><span class="seg-shape" style="width:'+pct(row.shape,total)+'"></span><span class="seg-multiple" style="width:'+pct(row.multiple,total)+'"></span><span class="seg-incomplete" style="width:'+pct(row.incomplete,total)+'"></span></div>';
  const label=row.chrom+"，W primary "+fmt(row.W_primary)+"，complete "+fmt(row.exact+row.shape+row.multiple)+"，incomplete "+fmt(row.incomplete);
  return '<button type="button" class="chrom-button" data-chrom="'+row.chrom+'" aria-pressed="false" aria-label="'+label+'"><strong>'+row.chrom+'<span>'+fmt(row.W_primary)+'</span></strong>'+mini+'<small>'+fmt(row.exact+row.shape+row.multiple)+' complete<br>'+fmt(row.incomplete)+' incomplete</small></button>';
 }).join("");
 host.querySelectorAll(".chrom-button").forEach(button=>button.addEventListener("click",()=>setActiveChrom(button.dataset.chrom,true)));
 document.getElementById("all-genome").addEventListener("click",()=>setActiveChrom("",true));
 const select=document.getElementById("fchr");
 for(let number=1;number<=22;number++){
  const option=document.createElement("option");option.value="chr"+number;option.textContent="chr"+number+" ("+fmt(CHROM_SUMMARIES["chr"+number].W_primary)+" W primary)";select.appendChild(option);
 }
}

function independentSignal(row,signal){
 if(signal==="recurrence")return row.has_recurrence;
 if(signal==="hidden")return row.hidden_positive;
 if(signal==="multi_hp")return row.hp_multiplicity===2;
 return true;
}

function matchesCoordinateQuery(row,rawQuery){
 const query=rawQuery.trim().toLowerCase().replace(/,/g,"");
 if(!query)return true;
 const coordinate=query.match(/^(chr(?:[1-9]|1[0-9]|2[0-2]))(?::(\d+)(?:-(\d+))?)?$/);
 if(!coordinate)return row.region.toLowerCase().includes(query);
 if(row.chrom.toLowerCase()!==coordinate[1])return false;
 if(!coordinate[2])return true;
 const start=Number(coordinate[2]);const end=Number(coordinate[3]||coordinate[2]);
 const queryStart=Math.min(start,end),queryEnd=Math.max(start,end);
 return row.end>=queryStart&&row.start<=queryEnd;
}

function applyFilters(){
 const chrom=document.getElementById("fchr").value;
 const topology=document.getElementById("ftopo").value;
 const evidence=document.getElementById("fevidence").value;
 const signal=document.getElementById("fsignal").value;
 const query=document.getElementById("fq").value;
 state.filtered=INDEX.filter(row=>
  (!chrom||row.chrom===chrom)&&
  (!topology||row.topology_class===topology)&&
  (!evidence||row.evidence_mode===evidence)&&
  independentSignal(row,signal)&&
  matchesCoordinateQuery(row,query)
 ).sort((a,b)=>(chromRank(a.chrom)-chromRank(b.chrom))||(a.start-b.start)||(a.end-b.end));
 state.limit=80;
 renderResultList();
}

function evidenceBadge(mode){
 const labels={full_and_partial:"full + partial",partial_only:"partial-only",full_only:"full-only",no_read_state_groups:"no read-state groups",not_applicable:"not applicable"};
 const cls=mode==="partial_only"?"warn":"";
 return '<span class="badge '+cls+'">'+(labels[mode]||mode)+'</span>';
}

function renderResultList(){
 let shown=state.filtered.slice(0,state.limit);
 const selectedRow=state.selectedRegion?state.filtered.find(row=>row.region===state.selectedRegion):null;
 const pinnedRegion=selectedRow&&!shown.some(row=>row.region===selectedRow.region)?selectedRow.region:null;
 if(pinnedRegion)shown=[selectedRow,...shown];
 document.getElementById("result-count").textContent=fmt(state.filtered.length)+" regions";
 const host=document.getElementById("result-list");
 host.innerHTML=shown.map(row=>
  '<button type="button" class="result-row" data-region="'+esc(row.region)+'" aria-current="'+(state.selectedRegion===row.region?"true":"false")+'">'+
   '<span class="result-title"><span>'+esc(row.region)+'</span><span>'+row.n_sSNV+' sSNV</span></span>'+
   '<span class="result-dimensions"><span>拓樸 '+esc(topologyShapeText(row))+'</span><span>辨識 '+esc(identifiabilityInfo(row).value)+'</span><span>端點距離 '+fmtSpan(row.endpoint_distance_bp)+'</span></span>'+
   '<span class="badges">'+topologyClassBadge(row.topology_class)+
    '<span class="badge">HP primary '+row.n_primary+'</span>'+evidenceBadge(row.evidence_mode)+
    (row.has_recurrence?'<span class="badge warn">recurrence</span>':"")+
    (row.region===pinnedRegion?'<span class="badge shape">pinned selection</span>':"")+
    '<span class="badge">CN '+esc(row.cn)+'</span></span>'+
  '</button>'
 ).join("");
 host.querySelectorAll(".result-row").forEach(button=>button.addEventListener("click",()=>selectRegion(button.dataset.region,button,true)));
 const load=document.getElementById("load-more");
 load.hidden=state.limit>=state.filtered.length;
 load.textContent="顯示更多 regions（"+fmt(Math.min(80,Math.max(0,state.filtered.length-state.limit)))+"）";
}

function syncHash(mode){
 if(mode==="none")return;
 const params=new URLSearchParams();
 Object.entries(FILTER_PARAMS).forEach(([id,key])=>{
  const value=document.getElementById(id).value.trim();if(value)params.set(key,value);
 });
 if(state.selectedRegion)params.set("region",state.selectedRegion);
 const hash=params.toString();const base=location.href.split("#")[0];
 const next=base+(hash?"#"+hash:"");if(next===location.href)return;
 history[mode==="push"?"pushState":"replaceState"](null,"",next);
}

function renderEmptyDetail(restoreFocus){
 state.selectedRegion=null;
 state.selectedRow=null;
 renderResultList();
 const detail=document.getElementById("detail");
 detail.innerHTML='<div class="detail-empty"><h3 id="detail-title">尚未選擇 region</h3><p>從左側結果或上方 chromosome overview 選取；detail 會先回答拓樸型態、可辨識度與基因體位置，再展開 read evidence。</p></div>';
 document.getElementById("detail-status").textContent="尚未選擇 region";
 if(restoreFocus)document.getElementById("fchr").focus({preventScroll:true});
}

function selectRegion(regionId,rowElement,restoreFocus,historyMode){
 const indexRow=INDEX.find(row=>row.region===regionId);
 if(!indexRow)return;
 try{
  if(!document.getElementById("fchr").value)setActiveChrom(indexRow.chrom,false);
  const region=getDetail(indexRow);
  state.selectedRegion=regionId;
  state.selectedRow=indexRow;
  renderResultList();
  renderDetail(region,indexRow);
  syncHash(historyMode||"push");
  if(restoreFocus){
   const detail=document.getElementById("detail");detail.focus({preventScroll:true});
   if(matchMedia("(max-width:860px)").matches)detail.scrollIntoView({block:"start"});
  }
 }catch(error){
  document.getElementById("detail").innerHTML='<div class="error"><h3>Region detail 無法載入</h3><p>'+esc(error.message)+'</p></div>';
 }
}

function primaryOf(region){return (region.lineages||[]).filter(line=>line.is_primary_lineage===true)}
function auxiliaryOf(region){return (region.lineages||[]).filter(line=>line.is_primary_lineage!==true)}

function verdictText(indexRow){
 if(indexRow.topology_class==="no_primary_lineage")return "此區沒有 mutation-bearing HP1／HP2 主重建單位；H3／H4／none／reference 僅能作輔助或控制。";
 if(indexRow.topology_class==="incomplete")return "此區至少一個 primary HP candidate set 未完整，因此 joint C／Topo unavailable；未完整 unit 不評估 edge stability，完整 sibling unit 仍可在其 unit scope 評估。";
 if(indexRow.topology_class==="exact_and_topology_unique")return "此區 candidate set 完整，exact candidate 與 topology shape 都唯一。";
 if(indexRow.topology_class==="topology_unique_exact_multiple")return "此區 candidate set 完整；exact candidates 不唯一，但 topology shape 唯一。";
 return "此區 candidate set 完整；exact candidates 與 topology shapes 都有多種相容解。";
}

function limitationText(region,indexRow){
 const limits=[];
 if(indexRow.evidence_mode==="partial_only")limits.push("所有 primary units 只由 overlapping partial reads 約束");
 if(indexRow.topology_class==="incomplete")limits.push("枚舉未完整");
 if(indexRow.has_recurrence)limits.push("存在 recurrence facet");
 if(region.cn&&region.cn!=="neutral")limits.push("CN "+region.cn+" sidecar");
 if(!limits.length)limits.push("仍是單一 bulk 的區域候選，不是細胞 lineage 確認");
 return limits.join("；");
}

function siteEvidence(region){
 const positions=region.positions||[];
 if(!positions.length)return '<p class="muted">此 region-view 沒有可展示的位點座標。</p>';
 const primaryLines=primaryOf(region);
 const evidenceLines=primaryLines.length?primaryLines:(region.lineages||[]);
 const evidenceScope=primaryLines.length
  ? "mutation-bearing HP1／HP2 primary units 加總"
  : "no-primary region：auxiliary／control units 加總；僅描述觀察，不進 C／Topo";
 const columns=positions.map(position=>{
  let ref=0,alt=0;
  evidenceLines.forEach(line=>{const value=(line.obs_col_coverage||{})[String(position)];if(value){ref+=Number(value[0]||0);alt+=Number(value[1]||0)}});
  const total=ref+alt;
  return {position:position,ref:ref,alt:alt,fraction:total?alt/total:0};
 });
 const headers=columns.map((item,index)=>'<th scope="col">S'+(index+1)+'<br><span class="mono">'+item.position+'</span></th>').join("");
 const states=columns.map(item=>'<td><span class="site-state '+(item.alt?"observed":"absent")+'">'+(item.alt?"ALT observed":"0 ALT")+'</span><br><span class="mono">'+fmt(item.alt)+' A / '+fmt(item.ref)+' R</span><br><span class="muted">read ALT '+(100*item.fraction).toFixed(1)+'%</span></td>').join("");
 return '<p class="scroll-cue" id="site-scroll-cue">可水平捲動；read ALT fraction 僅描述 read 證據，不是 CCF。</p>'+
  '<div class="scroll-region" role="region" tabindex="0" aria-label="位點 read 證據表" aria-describedby="site-scroll-cue">'+
  '<table class="site-table"><caption>Observed 位點證據（'+evidenceScope+'）</caption><thead><tr><th scope="col">Evidence</th>'+headers+'</tr></thead><tbody><tr><th scope="row">ALT / REF reads</th>'+states+'</tr></tbody></table></div>';
}

function nodeKind(node){
 if(node==="ROOT")return "root";
 return String(node).startsWith("H_")?"hidden":"observed";
}

function stateLabel(node){
 if(node==="ROOT")return "all-R";
 return String(node).startsWith("H_")?String(node).slice(2)+" H":String(node);
}

function shapeSignature(edges){
 const children={};const childNodes=new Set();
 (edges||[]).forEach(edge=>{(children[edge[0]]=children[edge[0]]||[]).push(edge[1]);childNodes.add(edge[1])});
 const roots=Object.keys(children).filter(node=>!childNodes.has(node));
 const walk=node=>"("+((children[node]||[]).map(walk).sort().join(""))+")";
 return (roots.length?roots.map(walk).sort().join("|"):"()");
}

function layoutNetwork(edges){
 const nodes=new Set();const incoming={};const outgoing={};
 (edges||[]).forEach(edge=>{nodes.add(edge[0]);nodes.add(edge[1]);(outgoing[edge[0]]=outgoing[edge[0]]||[]).push(edge[1]);incoming[edge[1]]=(incoming[edge[1]]||0)+1});
 if(!nodes.size)nodes.add("ROOT");
 const depth={};[...nodes].filter(node=>!incoming[node]).forEach(node=>depth[node]=0);
 for(let pass=0;pass<nodes.size+2;pass++){
  let changed=false;(edges||[]).forEach(edge=>{if(depth[edge[0]]!=null){const next=depth[edge[0]]+1;if(depth[edge[1]]==null||next>depth[edge[1]]){depth[edge[1]]=next;changed=true}}});if(!changed)break;
 }
 [...nodes].forEach(node=>{if(depth[node]==null)depth[node]=node==="ROOT"?0:(String(node).replace(/^H_/,"").match(/A/g)||[]).length});
 const levels={};[...nodes].forEach(node=>(levels[depth[node]]=levels[depth[node]]||[]).push(node));
 Object.values(levels).forEach(level=>level.sort());
 const maxAcross=Math.max(1,...Object.values(levels).map(level=>level.length));
 const width=Math.max(720,maxAcross*118+80);
 const maxDepth=Math.max(0,...Object.keys(levels).map(Number));
 const height=Math.max(170,(maxDepth+1)*90+45);
 const position={};
 Object.entries(levels).forEach(([level,nodesAtLevel])=>{
  const y=35+Number(level)*90;const gap=width/(nodesAtLevel.length+1);
  nodesAtLevel.forEach((node,index)=>position[node]={x:gap*(index+1),y:y});
 });
 return {nodes:[...nodes],position:position,width:width,height:height};
}

function edgeAcquisition(parent,child){
 const c=String(child).replace(/^H_/,"");
 const p=parent==="ROOT"?"R".repeat(c.length):String(parent).replace(/^H_/,"");
 const changes=[];for(let index=0;index<Math.min(p.length,c.length);index++)if(p[index]!==c[index]&&c[index]==="A")changes.push("S"+(index+1));
 return changes.length?"+"+changes.join(","):"";
}

function networkSvg(edges,edgeStates,title,recurrenceIndices){
 const layout=layoutNetwork(edges);
 const markerBase="network-arrow-"+(++networkCounter);
 const descriptionId=markerBase+"-description";
 const exactCandidate=Array.isArray(recurrenceIndices);
 const recurrenceSites=new Set((recurrenceIndices||[]).map(index=>"S"+(Number(index)+1)));
 const accessibleEdges=(edges||[]).map(edge=>{
  const key=edge[0]+">"+edge[1];const status=edgeStates[key]||"unevaluated";const acquisition=edgeAcquisition(edge[0],edge[1]);
  const repeated=[...recurrenceSites].filter(site=>acquisition.split(/[^A-Za-z0-9]+/).includes(site));
  return stateLabel(edge[0])+" → "+stateLabel(edge[1])+"；"+status+(acquisition?"；"+acquisition:"")+(repeated.length?"；repeated-acquisition "+repeated.join(","):"");
 });
 let accessibleSummary=accessibleEdges.length
  ? "Directed mutation-state edges: "+accessibleEdges.join("。")+"。"
  : "No stored candidate edge; ROOT-only state.";
 if(recurrenceSites.size)accessibleSummary+=" Candidate recurrence sites: "+[...recurrenceSites].join(", ")+".";
 let svg='<svg role="img" aria-label="'+esc(title)+'" aria-describedby="'+descriptionId+'" width="'+layout.width+'" height="'+layout.height+'" viewBox="0 0 '+layout.width+' '+layout.height+'">';
 svg+='<title>'+esc(title)+'</title>';
 svg+='<defs><marker id="'+markerBase+'-forced" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#536476"/></marker><marker id="'+markerBase+'-variable" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#9a5200"/></marker><marker id="'+markerBase+'-unevaluated" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#8b97a2"/></marker></defs>';
 (edges||[]).forEach(edge=>{
  const p=layout.position[edge[0]],c=layout.position[edge[1]];if(!p||!c)return;
  const key=edge[0]+">"+edge[1];const status=edgeStates[key]||"unevaluated";
  const stroke=status==="forced"?"#536476":status==="variable"?"#9a5200":"#8b97a2";
  const dash=status==="variable"?' stroke-dasharray="8 5"':status==="unevaluated"?' stroke-dasharray="2 6"':"";
  const dx=c.x-p.x,dy=c.y-p.y,length=Math.max(1,Math.hypot(dx,dy));
  const x1=p.x+dx/length*12,y1=p.y+dy/length*12,x2=c.x-dx/length*15,y2=c.y-dy/length*15;
  svg+='<line x1="'+x1+'" y1="'+y1+'" x2="'+x2+'" y2="'+y2+'" stroke="'+stroke+'" stroke-width="2" marker-end="url(#'+markerBase+'-'+status+')"'+dash+'><title>'+esc(status)+" directed edge "+esc(key)+'</title></line>';
  const label=edgeAcquisition(edge[0],edge[1]);const repeated=[...recurrenceSites].filter(site=>label.split(/[^A-Za-z0-9]+/).includes(site));
  if(label&&exactCandidate){
   const diagonal=Math.abs(dx)>12;
   const exactLabelX=c.x+(diagonal?(dx>0?14:-14):14);
   const exactAnchor=diagonal&&dx<0?"end":"start";
   const title=repeated.length?label+" repeated acquisition in this candidate":label+" acquisition";
   svg+='<text x="'+exactLabelX+'" y="'+(c.y-13)+'" text-anchor="'+exactAnchor+'" fill="'+(repeated.length?"#a9236d":"#536476")+'" font-size="10" font-weight="'+(repeated.length?"700":"400")+'"><title>'+esc(title)+'</title>'+esc(label)+'</text>';
  }
 });
 layout.nodes.forEach(node=>{
  const point=layout.position[node];const kind=nodeKind(node);
  if(kind==="root")svg+='<rect x="'+(point.x-8)+'" y="'+(point.y-8)+'" width="16" height="16" rx="2" fill="#5c6875"><title>all-reference mutation state</title></rect>';
  else if(kind==="hidden")svg+='<circle cx="'+point.x+'" cy="'+point.y+'" r="9" fill="#fff" stroke="#7040a0" stroke-width="2" stroke-dasharray="4 3"><title>inferred hidden mutation state; not directly observed</title></circle>';
  else svg+='<circle cx="'+point.x+'" cy="'+point.y+'" r="9" fill="#075ea8"><title>observed read-state genotype</title></circle>';
  svg+='<text x="'+point.x+'" y="'+(point.y+25)+'" text-anchor="middle" fill="'+(kind==="hidden"?"#7040a0":"#3d5061")+'" font-size="11" font-family="ui-monospace,monospace">'+esc(stateLabel(node))+'</text>';
 });
 return svg+'</svg><p class="sr-only" id="'+descriptionId+'">'+esc(accessibleSummary)+'</p>';
}

function candidateExplorer(line,laneId){
 const trees=line.trees||[];
 if(!trees.length)return '<div class="scope-warning">此 unit 沒有可畫的 stored candidate edges；這是資料狀態，不是空白錯誤。</div>';
 const analysisComplete=line.analysis_candidate_set_complete===true;
 const fullDisplay=analysisComplete&&line.display_trees_complete===true&&Number(line.n_trees_stored)===Number(line.n_trees);
 const edgeCounts={};trees.forEach(tree=>(tree.edges||[]).forEach(edge=>{const key=edge[0]+">"+edge[1];edgeCounts[key]=(edgeCounts[key]||0)+1}));
 const unionEdges=Object.keys(edgeCounts).map(key=>key.split(">"));
 const edgeStates={};Object.keys(edgeCounts).forEach(key=>edgeStates[key]=fullDisplay?(edgeCounts[key]===trees.length?"forced":"variable"):"unevaluated");
 const groups={};trees.forEach((tree,index)=>{const signature=shapeSignature(tree.edges||[]);(groups[signature]=groups[signature]||[]).push(index)});
 const signatures=Object.keys(groups);
 const storedShapeCount=signatures.length;const declaredStoredShapes=Number(line.n_distinct_shapes_stored);
 if(Number.isFinite(declaredStoredShapes)&&storedShapeCount!==declaredStoredShapes)throw new Error("Stored topology shape count mismatch for "+line.fam_label);
 const exactShapeCount=Number(line.n_distinct_shapes_exact);
 if(analysisComplete&&Number.isFinite(exactShapeCount)&&storedShapeCount>exactShapeCount)throw new Error("Stored topology shapes exceed exact total for "+line.fam_label);
 if(fullDisplay&&storedShapeCount!==exactShapeCount)throw new Error("Exact topology shape count mismatch for "+line.fam_label);
 const allExactShapesStored=analysisComplete&&Number.isFinite(exactShapeCount)&&storedShapeCount===exactShapeCount;
 const exactShapeText=line.n_distinct_shapes_exact==null?"unavailable":fmt(line.n_distinct_shapes_exact);
 const shapePrefix=allExactShapesStored?"Shape":analysisComplete?"Stored shape":"Prefix shape";
 const shapeScope=fullDisplay
  ? "Shape controls cover the complete exact set ("+storedShapeCount+" / "+exactShapeText+")."
  : analysisComplete
   ? allExactShapesStored
    ? "Shape controls cover all "+storedShapeCount+" exact topology groups; exact candidate trees within a shape may still be outside this stored display subset."
    : "Shape controls cover "+storedShapeCount+" stored groups from an exact total of "+exactShapeText+"; unshown exact shapes remain outside this display."
   : "Shape controls cover "+storedShapeCount+" groups in the enumerated prefix; the exact universe total is unknown.";
 const recurrenceCounts={};trees.forEach(tree=>(tree.recurrence||[]).forEach(index=>{const site="S"+(Number(index)+1);recurrenceCounts[site]=(recurrenceCounts[site]||0)+1}));
 const recurrenceSummary=Object.entries(recurrenceCounts).map(([site,count])=>site+" in "+count+" / "+trees.length+" stored candidates").join("；");
 const candidatePrefix=fullDisplay?"Exact":analysisComplete?"Stored exact":"Prefix";
 CANDIDATE_CACHE[laneId]={trees:trees,groups:signatures.map(signature=>groups[signature]),fullDisplay:fullDisplay,edgeStates:edgeStates,shapePrefix:shapePrefix,candidatePrefix:candidatePrefix,famLabel:line.fam_label};
 const networkTitle="Candidate-space network for "+line.fam_label+"; "+(fullDisplay?"forced and variable edges evaluated across the complete displayed exact set":"edge stability not evaluated because displayed candidates are incomplete");
 let out='<div class="network-card"><div class="network-head"><div><b>候選空間 network</b><p>'+esc(fullDisplay?"完整 exact set 的聯集；實線為全候選共有，虛線為候選間變動。":"只顯示 stored candidates 的聯集；所有邊以點線標為 not evaluated，不宣稱 forced。")+'</p></div><span class="badge '+(fullDisplay?"exact":"incomplete")+'">'+(fullDisplay?"edge scope complete":"edge scope incomplete")+'</span></div>';
 out+='<div class="network-scroll" role="region" tabindex="0" aria-label="'+esc(line.fam_label)+' 候選空間 network，可水平捲動">'+networkSvg(unionEdges,edgeStates,networkTitle)+'</div>';
 out+='<div class="network-legend"><span>箭頭＝mutation-state transition 方向</span><span><span class="line-swatch"></span>forced across complete exact set</span><span><span class="line-swatch variable"></span>varies between candidates</span><span><span class="line-swatch unevaluated"></span>not evaluated</span><span>實心節點 observed</span><span>空心虛線 H inferred</span><span>聯集圖省略 +S_i label；逐邊 acquisition 見下方 single candidate</span></div></div>';
 out+='<div class="shape-scope"><b>Shape display scope</b><span>'+esc(shapeScope)+'</span></div>';
 if(recurrenceSummary)out+='<div class="recurrence-note"><b>Stored candidate recurrence annotation</b><span>'+esc(recurrenceSummary)+'。這是 candidate-level repeated acquisition；是否構成 recurrence facet 仍依 primary L1 class。</span></div>';
 out+='<div class="shape-tabs" role="group" aria-label="'+esc(line.fam_label+' '+shapePrefix+' topology groups')+'">'+signatures.map((signature,index)=>'<button class="shape-tab" type="button" aria-label="'+esc(line.fam_label+' '+shapePrefix+' '+(index+1)+', '+groups[signature].length+' stored exact candidates')+'" aria-pressed="'+(index===0?"true":"false")+'" data-lane="'+laneId+'" data-shape="'+index+'">'+shapePrefix+' '+(index+1)+' · '+groups[signature].length+' stored exact</button>').join("")+'</div>';
 out+='<div class="candidate-view" id="'+laneId+'-candidate"></div>';
 return out;
}

function renderCandidateShape(laneId,shapeIndex,treeIndex,restoreTreeFocus){
 const payload=CANDIDATE_CACHE[laneId];if(!payload)return;
 const candidateIndices=payload.groups[shapeIndex]||[];const activeIndex=candidateIndices[Math.min(treeIndex||0,candidateIndices.length-1)]||0;
 const tree=payload.trees[activeIndex]||{edges:[]};
 const recurrenceSites=(tree.recurrence||[]).map(index=>"S"+(Number(index)+1));
 const recurrenceNote=recurrenceSites.length?'<div class="recurrence-note"><b>Repeated acquisition in this stored candidate</b><span>'+esc(recurrenceSites.join(", "))+'；圖中相關 acquisition label 以洋紅粗體標出，accessible edge description 亦逐邊明示 repeated acquisition。</span></div>':"";
 const localStates={};(tree.edges||[]).forEach(edge=>localStates[edge[0]+">"+edge[1]]=payload.fullDisplay?payload.edgeStates[edge[0]+">"+edge[1]]:"unevaluated");
 const host=document.getElementById(laneId+"-candidate");
 host.innerHTML='<div class="tree-tabs" role="group" aria-label="'+esc(payload.famLabel+' '+payload.candidatePrefix+' candidates in '+payload.shapePrefix+' '+(shapeIndex+1))+'">'+candidateIndices.map((index,position)=>'<button class="tree-tab" type="button" aria-label="'+esc(payload.famLabel+' '+payload.shapePrefix+' '+(shapeIndex+1)+' '+payload.candidatePrefix+' candidate '+(index+1))+'" aria-pressed="'+(index===activeIndex?"true":"false")+'" data-lane="'+laneId+'" data-shape="'+shapeIndex+'" data-tree="'+position+'">'+payload.candidatePrefix+' #'+(index+1)+'</button>').join("")+'</div>'+recurrenceNote+
  '<div class="network-card"><div class="network-head"><div><b>'+payload.candidatePrefix+' candidate #'+(activeIndex+1)+'</b><p>此圖不做候選排序；stored order 只供逐一檢視。</p></div></div>'+
  '<div class="network-scroll" role="region" tabindex="0" aria-label="'+esc(payload.famLabel+' '+payload.shapePrefix+' '+(shapeIndex+1)+' '+payload.candidatePrefix+' candidate '+(activeIndex+1)+' mutation-state network，可水平捲動')+'">'+networkSvg(tree.edges||[],localStates,payload.famLabel+" "+payload.shapePrefix+" "+(shapeIndex+1)+" "+payload.candidatePrefix+" mutation-state candidate "+(activeIndex+1),tree.recurrence||[])+'</div>'+
  '<div class="tree-caption">'+payload.shapePrefix+' '+(shapeIndex+1)+' · '+payload.candidatePrefix.toLowerCase()+' candidate '+(activeIndex+1)+' / stored '+payload.trees.length+(payload.fullDisplay?" · complete display":" · display subset")+'</div></div>';
 host.querySelectorAll(".tree-tab").forEach(button=>button.addEventListener("click",()=>renderCandidateShape(button.dataset.lane,Number(button.dataset.shape),Number(button.dataset.tree),true)));
 document.querySelectorAll('.shape-tab[data-lane="'+laneId+'"]').forEach((button,index)=>button.setAttribute("aria-pressed",index===shapeIndex?"true":"false"));
 if(restoreTreeFocus){const active=host.querySelector('.tree-tab[aria-pressed="true"]');if(active)active.focus()}
}

function readGroups(line){
 const full=line.obs_populations||{};const partial=line.obs_subreads||{};
 const rows=[];
 Object.entries(full).sort((a,b)=>b[1]-a[1]).forEach(item=>rows.push([item[0],"full-span observed",item[1]]));
 Object.entries(partial).sort((a,b)=>b[1]-a[1]).forEach(item=>rows.push([item[0],"partial overlap",item[1]]));
 if(!rows.length)return '<p>沒有 stored read-state groups。</p>';
 return '<p class="scroll-cue">可水平捲動；X 表示該 read 未覆蓋該位點。</p><div class="scroll-region" role="region" tabindex="0" aria-label="'+esc(line.fam_label)+' raw read-state groups"><table class="read-table"><caption>'+esc(line.fam_label)+' raw read-state groups</caption><thead><tr><th scope="col">State</th><th scope="col">Evidence mode</th><th scope="col" class="num">Reads</th></tr></thead><tbody>'+rows.map(row=>'<tr><td class="mono">'+esc(row[0])+'</td><td>'+row[1]+'</td><td class="num">'+fmt(row[2])+'</td></tr>').join("")+'</tbody></table></div>';
}

let laneCounter=0;
function laneHtml(line,isAuxiliary){
 const laneId="lane-"+(++laneCounter);
 const family=String(line.family||"none");const familyClass=family==="2"?"family-2":"";
 const role=isAuxiliary?"auxiliary":family==="1"?"primary HP1":"primary HP2";
 const complete=line.analysis_candidate_set_complete===true;
 const displayComplete=line.display_trees_complete===true;
 const candidateLabel=complete?"Exact candidates":"Enumerated prefix candidates";
 const shapeLabel=complete?"Exact topology shapes":"Stored prefix shapes";
 const shapeValue=complete?line.n_distinct_shapes_exact:line.n_distinct_shapes_stored;
 const storageLabel=complete?"Stored / exact total":"Stored / enumerated prefix";
 let out='<section class="lane '+familyClass+(isAuxiliary?" auxiliary":"")+'" id="'+laneId+'"><div class="lane-head"><div class="lane-title"><b>'+esc(line.fam_label)+'</b><span class="badge">'+esc(role)+'</span><span class="badge '+(complete?"exact":"incomplete")+'">'+(complete?"candidate complete":"candidate incomplete")+'</span></div><span class="badge">'+esc(line.verification_status)+'</span></div><div class="lane-body">';
 out+='<div class="candidate-summary"><div class="candidate-stat"><span>'+candidateLabel+'</span><strong>'+fmt(line.n_trees)+'</strong></div><div class="candidate-stat"><span>'+shapeLabel+'</span><strong>'+(shapeValue==null?"—":fmt(shapeValue))+'</strong></div><div class="candidate-stat"><span>'+storageLabel+'</span><strong>'+fmt(line.n_trees_stored)+' / '+fmt(line.n_trees)+'</strong></div><div class="candidate-stat"><span>Read states</span><strong>'+fmt(line.n_full_pops)+' full · '+fmt(line.n_partial)+' partial</strong></div></div>';
 if(!complete||!displayComplete)out+='<div class="scope-warning"><b>候選／展示範圍未完整。</b> '+(complete?"Exact universe 已知，但本頁只存展示子集。":"Exact universe total unknown；x / x 只表示 stored rows 與 enumerated prefix 相等，不代表候選已收齊。")+' 這個 unit 不宣稱 forced edge；network 中全部 edge stability 標為 not evaluated。</div>';
 if(!isAuxiliary)out+=candidateExplorer(line,laneId);
 else out+='<p class="muted">Auxiliary unit 不納入 W_primary、C 或 Topo；不畫入 primary candidate network。</p>';
 out+='<details class="drawer"><summary>Raw read-state groups 與 L0–L3 trace</summary><div class="drawer-body">'+readGroups(line)+'<div class="trace">'+(line.trace||[]).map(item=>'<div>'+esc(item)+'</div>').join("")+'</div><p class="mono">analysis digest: '+esc(line.analysis_tree_digest_sha256||"not available")+'</p></div></details>';
 return out+"</div></section>";
}

function sidecars(region,indexRow){
 const cn=DATA.copy_number_contract||{};const l3=DATA.l3||{};const readTag=DATA.read_tag_census||{};
 const primaryCnVerdicts=[...new Set(primaryOf(region).map(line=>line.L2_cn_verdict).filter(Boolean))];
 const verdictText=primaryCnVerdicts.length?primaryCnVerdicts.join("；"):"not applicable (no primary units)";
 const cnText=cn.availability==="measured"
  ?"Region CN="+esc(region.cn||"unavailable")+"；source="+esc(cn.source||"unspecified")+"；availability=measured；primary L2 verdict="+esc(verdictText)+"。只作 recurrence/confound sidecar，不排序候選。"
  :"CN source="+esc(cn.source||"unspecified")+"；availability="+esc(cn.availability||"unavailable")+"；primary L2 verdict="+esc(verdictText)+"。缺失不得默認 neutral。";
 return '<div class="sidecars"><article class="sidecar"><b>CN sidecar</b><p>'+cnText+'</p></article>'+
  '<article class="sidecar"><b>PS QC sidecar</b><p>Dataset mixed-PS regions '+fmt(SAMPLE_REC.mixed_PS_regions)+'；read-tag exact '+fmt(SAMPLE_REC.read_tag_exact_matches)+' / '+fmt(SAMPLE_REC.read_tag_exposures)+'。PS 不作 topology edge；此 region-view 未提供 region-level PS 展開。</p></article>'+
  '<article class="sidecar"><b>L3 methylation sidecar</b><p>Status '+esc(l3.status)+'；role '+esc(l3.role)+'。只允許 negative screen / residual flag；禁止 tree ranking、lineage 或 clone confirmation。</p></article></div>';
}

function regionDimensionCards(indexRow){
 const identifiability=identifiabilityInfo(indexRow);
 const topologyNote=indexRow.topology_class==="no_primary_lineage"
  ? "沒有 mutation-bearing HP1／HP2，因此不進 regional shape 分母。"
  : indexRow.topology_class==="incomplete"
   ? "Candidate set 未完整，Topo 保持 unavailable。"
   : "Topo="+fmt(indexRow.Topo)+"；比較無節點標籤的 regional mutation-state shapes。";
 return '<div class="region-dimensions" aria-label="此 region 的三個閱讀維度">'+
  '<article class="region-dimension"><span>01 · 拓樸型態</span><strong>'+esc(topologyShapeText(indexRow))+'</strong><p>'+esc(topologyNote)+'</p></article>'+
  '<article class="region-dimension"><span>02 · 可辨識度</span><strong>'+esc(identifiability.value)+'</strong><p>'+esc(identifiability.note)+'</p></article>'+
  '<article class="region-dimension"><span>03 · 基因體位置</span><strong class="mono">'+esc(indexRow.region)+'</strong><p>'+fmt(indexRow.n_sSNV)+' retained sSNV；端點距離 '+fmtSpan(indexRow.endpoint_distance_bp)+'（end - start）；inclusive span '+fmtSpan(indexRow.inclusive_span_bp)+'。</p></article>'+
 '</div>';
}

function renderDetail(region,indexRow){
 const detail=document.getElementById("detail");
 Object.keys(CANDIDATE_CACHE).forEach(key=>delete CANDIDATE_CACHE[key]);
 laneCounter=0;
 const info=classInfo[indexRow.topology_class]||classInfo.no_primary_lineage;
 const primary=primaryOf(region),auxiliary=auxiliaryOf(region);
 const currentIndex=state.filtered.findIndex(row=>row.region===indexRow.region);
 const prev=currentIndex>0?state.filtered[currentIndex-1]:null;
 const next=currentIndex>=0&&currentIndex<state.filtered.length-1?state.filtered[currentIndex+1]:null;
 const cValue=indexRow.C==null?"—":fmt(indexRow.C);const topoValue=indexRow.Topo==null?"—":fmt(indexRow.Topo);
 let out='<div class="detail-toolbar"><div class="routes"><button class="button" type="button" id="back-results">返回結果</button><button class="button" type="button" id="prev-region" '+(prev?"":"disabled")+'>上一區</button><button class="button" type="button" id="next-region" '+(next?"":"disabled")+'>下一區</button></div><span class="mono">'+esc(indexRow.chrom)+'</span></div><div class="detail-body">';
 out+='<div class="verdict '+(indexRow.topology_class==="incomplete"?"incomplete":"")+'"><div><p class="kicker">Region verdict</p><h3 id="detail-title" class="mono">'+esc(region.region)+'</h3><p>'+verdictText(indexRow)+'</p><div class="facet-row">'+topologyClassBadge(indexRow.topology_class)+evidenceBadge(indexRow.evidence_mode)+'<span class="badge">Primary HP '+indexRow.n_primary+'</span><span class="badge">Auxiliary '+indexRow.n_auxiliary+'</span>'+(indexRow.has_recurrence?'<span class="badge warn">recurrence facet</span>':'')+'<span class="badge">CN '+esc(region.cn||"unavailable")+'</span></div></div><div class="ct-readout"><div class="ct-box"><span>C exact</span><strong>'+cValue+'</strong></div><div class="ct-box"><span>Topo shapes</span><strong>'+topoValue+'</strong></div></div></div>';
 out+=regionDimensionCards(indexRow);
 out+='<div class="assertion-grid"><article class="assertion"><b>Observed</b><p>'+fmt(region.n_sSNV)+' retained sSNV；'+fmt(region.n_full_cov_reads)+' reads 橫跨整個 region 位點集合。</p></article><article class="assertion"><b>Inferred</b><p>'+fmt(primary.reduce((sum,line)=>sum+Number(line.n_hidden||0),0))+' hidden states across primary units；只代表 solver 所需中間 mutation states。</p></article><article class="assertion"><b>Limit</b><p>'+esc(limitationText(region,indexRow))+'。</p></article></div>';
 out+='<section class="subsection" aria-labelledby="evidence-title"><div class="subsection-head"><h4 id="evidence-title">Observed site evidence</h4><span class="subsection-note">read ALT fraction is descriptive, not CCF</span></div>'+siteEvidence(region)+'</section>';
 out+='<section class="subsection" aria-labelledby="candidate-title"><div class="subsection-head"><h4 id="candidate-title">Primary HP candidate networks</h4><span class="subsection-note">HP1 / HP2 independently reconstructed</span></div>';
 out+=primary.length?primary.map(line=>laneHtml(line,false)).join(""):'<div class="scope-warning">此區沒有 primary mutation-bearing HP1／HP2；不計入 W_primary，也不計 C／Topo。</div>';
 out+='</section>';
 if(auxiliary.length)out+='<details class="drawer"><summary>Auxiliary／control units（'+auxiliary.length+'）</summary><div class="drawer-body">'+auxiliary.map(line=>laneHtml(line,true)).join("")+'</div></details>';
 out+='<section class="subsection" aria-labelledby="sidecar-title"><div class="subsection-head"><h4 id="sidecar-title">Independent sidecars</h4><span class="subsection-note">never topology edges; never candidate ranking</span></div>'+sidecars(region,indexRow)+'</section>';
 out+='</div>';
 detail.innerHTML=out;
 detail.querySelector("#back-results").addEventListener("click",()=>{const target=document.querySelector('.result-row[aria-current="true"]')||document.getElementById("fq");target.focus({preventScroll:true});document.getElementById("results-panel").scrollIntoView({block:"start"})});
 if(prev)detail.querySelector("#prev-region").addEventListener("click",()=>selectRegion(prev.region,null,true));
 if(next)detail.querySelector("#next-region").addEventListener("click",()=>selectRegion(next.region,null,true));
 detail.querySelectorAll(".shape-tab").forEach(button=>button.addEventListener("click",()=>renderCandidateShape(button.dataset.lane,Number(button.dataset.shape),0,false)));
 primary.forEach((line,index)=>{
  const lane=detail.querySelectorAll('.lane:not(.auxiliary)')[index];if(!lane)return;
  const laneId=lane.id;if(CANDIDATE_CACHE[laneId])renderCandidateShape(laneId,0,0,false);
 });
 document.getElementById("detail-status").textContent="已載入 "+region.region+"；"+info.label;
}

function renderSourceLinks(){
 const src=DATA.source;
 const fileHref=path=>"file://"+encodeURI(path);
 document.getElementById("source-links").innerHTML=
  '<div class="raw-link"><a href="'+fileHref(src.region_view)+'">開啟原始 layered region-view JSON</a><small>SHA-256 <span class="mono">'+esc(src.region_view_sha256)+'</span></small></div>'+
  '<div class="raw-link"><a href="'+fileHref(src.layered_reconstruction)+'">開啟原始 layered reconstruction JSON</a><small>SHA-256 <span class="mono">'+esc(src.layered_reconstruction_sha256)+'</span></small></div>'+
  '<div class="raw-link"><a href="'+fileHref(src.machine_summary)+'">開啟 canonical machine summary JSON</a><small>SHA-256 <span class="mono">'+esc(src.machine_summary_sha256)+'</span></small></div>'+
  '<div class="raw-link"><a href="'+fileHref(src.backbone_comparison)+'">開啟 backbone sensitivity comparison JSON</a><small>SHA-256 <span class="mono">'+esc(src.backbone_comparison_sha256)+'</span></small></div>'+
  '<div class="raw-link"><strong>Verification</strong><small>all_pass='+DATA.summary.all_pass+' · success <span class="mono">'+esc(DATA.summary.success_sha256)+'</span> · verifier <span class="mono">'+esc(DATA.summary.verification_sha256)+'</span></small></div>';
}

function setupControls(){
 document.querySelectorAll("#fevidence option[value]").forEach(option=>{
  if(!option.value)return;
  const count=INDEX.filter(row=>row.evidence_mode===option.value).length;
  option.textContent=option.textContent+"（"+fmt(count)+"）";
  option.disabled=count===0;
 });
 ["fchr","ftopo","fevidence","fsignal"].forEach(id=>document.getElementById(id).addEventListener("change",()=>{
  if(id==="fchr"){
   const chrom=document.getElementById("fchr").value;updateChromButtonState(chrom);
  }
  applyFilters();syncHash("push");
 }));
 document.getElementById("fq").addEventListener("input",()=>{applyFilters();syncHash("replace")});
 document.getElementById("load-more").addEventListener("click",()=>{state.limit+=80;renderResultList()});
 document.getElementById("copy-link").addEventListener("click",async event=>{
  const button=event.currentTarget;let ok=false;
  try{await navigator.clipboard.writeText(location.href);ok=true}catch(error){
   const area=document.createElement("textarea");area.value=location.href;area.style.position="fixed";area.style.opacity="0";document.body.appendChild(area);area.select();ok=document.execCommand("copy");area.remove();
  }
  button.textContent=ok?"已複製檢視連結":"請複製網址列";setTimeout(()=>button.textContent="複製目前檢視連結",1600);
 });
}

function restoreFromHash(restoreFocus){
 const params=new URLSearchParams(location.hash.slice(1));
 [["ftopo","topo"],["fevidence","evidence"],["fsignal","signal"]].forEach(([id,key])=>{
  const value=params.get(key);const select=document.getElementById(id);
  select.value=value&&[...select.options].some(option=>option.value===value)?value:"";
 });
 document.getElementById("fq").value=params.get("q")||"";
 const wanted=params.get("region");const row=wanted?INDEX.find(item=>item.region===wanted):null;
 const requestedChrom=params.get("chr")||(row?row.chrom:"");
 state.selectedRegion=null;state.selectedRow=null;
 setActiveChrom(/^chr(?:[1-9]|1[0-9]|2[0-2])$/.test(requestedChrom)?requestedChrom:"",false);
 if(row)selectRegion(wanted,null,restoreFocus,"none");else renderEmptyDetail(restoreFocus);
}

function init(){
 renderMetrics();renderGenome();renderDimensions();renderSourceLinks();setupControls();
 restoreFromHash(false);
 syncHash("replace");
 window.addEventListener("popstate",()=>restoreFromHash(true));
}
init();
</script>
</body>
</html>
"""

    rendered = (
        template.replace("__SUMMARY_SHA__", summary_sha)
        .replace("__REGION_SHA__", actual_rv_sha)
        .replace("__COMPARISON_SHA__", comparison_sha)
        .replace("__SAMPLE__", html.escape(sample))
        .replace("__BIO_SAMPLE__", html.escape(BIOLOGICAL_SAMPLE[sample]))
        .replace("__RUN_ID__", html.escape(run_id))
        .replace("__MIXED_PS__", f'{int(sample_record["mixed_PS_regions"]):,}')
        .replace("__BACKBONE_VERDICT__", html.escape(str(comparison["verdict"]).replace("_", " ").upper()))
        .replace("__SAMPLE_BACKBONE_VERDICT__", html.escape(str(sample_comparison["verdict"]).replace("_", " ").upper()))
        .replace("__RETAINED_JACCARD__", f'{float(sample_comparison["retained_position_jaccard"]):.3f}')
        .replace("__PRIMARY_JACCARD__", f'{float(sample_comparison["primary_unit_key_jaccard"]):.3f}')
        .replace("__TOPOLOGY_CONCORDANCE__", f'{float(sample_comparison["shared_unit_topology_digest_concordance"]):.3f}')
        .replace("__MAX_DELTA_PP__", f'{sample_max_delta:.2f}')
        .replace("__SAMPLE_OPTIONS__", sample_options)
        .replace("__PAGE_DATA__", json_text(page_data))
        .replace("__CHUNK_SCRIPTS__", "\n".join(chunk_scripts))
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(rendered, encoding="utf-8")
    print(
        "OK wrote {} ({} bytes; {} regions; summary_sha={}; region_sha={})".format(
            out_path,
            len(rendered.encode("utf-8")),
            len(regions),
            summary_sha,
            actual_rv_sha,
        )
    )


if __name__ == "__main__":
    main()
