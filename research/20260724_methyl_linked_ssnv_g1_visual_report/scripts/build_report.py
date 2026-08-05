#!/usr/bin/env python3
"""Build the all-sample methyl-group × linked-sSNV G1 evidence atlas.

The report intentionally separates three evidence layers:
1. stable focal-ALT methyl groups,
2. methyl-group × linked-partner allele co-segregation,
3. undirected read linkage within exact PS×HP W containers.

No layer is promoted to a cellular clone, causal methylation effect, or
identifiable mutation order.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
WORK = REPO / "research/20260724_methyl_linked_ssnv_g1_visual_report"
DATA_DIR = WORK / "data"
RESULTS_DIR = WORK / "results"

COOCCURRENCE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
    "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity"
)
PAIR_RESULTS = COOCCURRENCE_ROOT / "methyl_ssnv_pair_results.tsv.gz"
SUMMARY = COOCCURRENCE_ROOT / "summary.json"
ASSIGNMENTS = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
    "all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/"
    "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
)
CLAIM_CONTRACT = (
    REPO
    / "research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md"
)
STRICT_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260723_production_exact_ps_strict_read_linkage"
)
STRICT_JOIN_AUDIT = (
    REPO
    / "research/20260723_production_exact_ps_strict_read_linkage/results/"
    "20260724_exactPS_strict與methylM1_join盤點_01.md"
)

OUTPUT_DATA = DATA_DIR / "report_data.json"
OUTPUT_HTML = WORK / "20260724_methyl_linked_sSNV_G1全樣本關聯圖譜_01.standalone.html"
OUTPUT_RECEIPT = RESULTS_DIR / "20260724_G1全樣本HTML建置收據_01.json"

SAMPLES = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]

# Independent whole-scope recount from site/pair authorities. These values are
# asserted again against the global source summary and pair-table totals.
SAMPLE_FUNNEL = {
    "HCC1395": {
        "all_sites": 79687,
        "m1": 12838,
        "m2_evaluable": 37,
        "m2_pass": 4,
        "m2_with_partner": 2,
        "exact_focal_sites": 2,
        "exact_pairs": 2,
        "by": 1,
        "g1": 1,
    },
    "HCC1395_DORADO": {
        "all_sites": 79739,
        "m1": 14789,
        "m2_evaluable": 68,
        "m2_pass": 13,
        "m2_with_partner": 1,
        "exact_focal_sites": 0,
        "exact_pairs": 0,
        "by": 0,
        "g1": 0,
    },
    "COLO829": {
        "all_sites": 37788,
        "m1": 3579,
        "m2_evaluable": 0,
        "m2_pass": 0,
        "m2_with_partner": 0,
        "exact_focal_sites": 0,
        "exact_pairs": 0,
        "by": 0,
        "g1": 0,
    },
    "H1437": {
        "all_sites": 77080,
        "m1": 10187,
        "m2_evaluable": 41,
        "m2_pass": 16,
        "m2_with_partner": 5,
        "exact_focal_sites": 1,
        "exact_pairs": 2,
        "by": 0,
        "g1": 0,
    },
    "H2009": {
        "all_sites": 154465,
        "m1": 54644,
        "m2_evaluable": 1603,
        "m2_pass": 816,
        "m2_with_partner": 507,
        "exact_focal_sites": 122,
        "exact_pairs": 138,
        "by": 8,
        "g1": 5,
    },
    "HCC1937": {
        "all_sites": 18690,
        "m1": 1938,
        "m2_evaluable": 86,
        "m2_pass": 50,
        "m2_with_partner": 19,
        "exact_focal_sites": 2,
        "exact_pairs": 4,
        "by": 0,
        "g1": 0,
    },
    "HCC1954": {
        "all_sites": 22400,
        "m1": 4867,
        "m2_evaluable": 32,
        "m2_pass": 20,
        "m2_with_partner": 14,
        "exact_focal_sites": 1,
        "exact_pairs": 1,
        "by": 1,
        "g1": 1,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-html", type=Path, default=OUTPUT_HTML)
    parser.add_argument("--output-data", type=Path, default=OUTPUT_DATA)
    parser.add_argument("--output-receipt", type=Path, default=OUTPUT_RECEIPT)
    return parser.parse_args()


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() == "true"


def as_float(value: Any) -> float | None:
    text = str(value).strip()
    return None if text == "" else float(text)


def as_int(value: Any) -> int | None:
    text = str(value).strip()
    return None if text == "" else int(float(text))


def parse_json(value: Any, default: Any) -> Any:
    text = str(value).strip()
    return default if text == "" else json.loads(text)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_pair_authority() -> tuple[list[dict[str, str]], dict[str, dict[str, int]]]:
    formal_rows: list[dict[str, str]] = []
    census = {
        sample: {
            "pair_rows": 0,
            "testable": 0,
            "exact_pairs": 0,
            "by": 0,
            "g1": 0,
        }
        for sample in SAMPLES
    }
    with gzip.open(PAIR_RESULTS, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            census[sample]["pair_rows"] += 1
            census[sample]["testable"] += int(as_bool(row["endpoint_a_testable"]))
            census[sample]["exact_pairs"] += int(str(row["endpoint_a_q_global_by"]).strip() != "")
            census[sample]["by"] += int(as_bool(row["endpoint_a_exact_by_discovery"]))
            is_g1 = as_bool(row["endpoint_a_formal_pair_by_confirmed"])
            census[sample]["g1"] += int(is_g1)
            if is_g1:
                formal_rows.append(row)
    formal_rows.sort(
        key=lambda r: (
            SAMPLES.index(r["sample"]),
            int(r["focal_chrom"].replace("chr", "")),
            int(r["focal_pos"]),
        )
    )
    return formal_rows, census


def load_assignments(
    target_keys: set[tuple[str, str, int]],
) -> tuple[dict[tuple[str, str, int], dict[str, Any]], Counter[str]]:
    targets: dict[tuple[str, str, int], dict[str, Any]] = {}
    counts: Counter[str] = Counter()
    with gzip.open(ASSIGNMENTS, "rt") as handle:
        for line in handle:
            record = json.loads(line)
            key = (record["dataset"], record["chrom"], int(record["pos"]))
            counts[record["dataset"]] += 1
            if key in target_keys:
                targets[key] = record
    missing = sorted(target_keys - targets.keys())
    if missing:
        raise AssertionError(f"Missing target assignments: {missing}")
    return targets, counts


def matrix_for_record(record: dict[str, Any]) -> dict[str, Any]:
    matrix_path = Path(record["primary_artifacts"]["methylation_matrix"]["path"])
    expected_sha = record["primary_artifacts"]["methylation_matrix"]["sha256"]
    observed_sha = sha256_file(matrix_path)
    if observed_sha != expected_sha:
        raise AssertionError(f"Matrix SHA mismatch: {matrix_path}")

    with matrix_path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        positions = [int(value) for value in header[1:]]
        raw: dict[str, list[float | None]] = {}
        for row in reader:
            raw[str(row[0])] = [
                None if value == "NA" or value == "" else float(value) for value in row[1:]
            ]

    core_reads = record["core_reads"]
    missing_core = [read["read_id"] for read in core_reads if str(read["read_id"]) not in raw]
    if missing_core:
        raise AssertionError(f"Core read IDs absent from matrix: {missing_core[:5]}")

    labels = sorted({read["label"] for read in core_reads})
    group_reads: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for read in core_reads:
        group_reads[read["label"]].append(read)

    column_scores: list[tuple[float, int, dict[str, float], dict[str, int]]] = []
    for column_index, _ in enumerate(positions):
        means: dict[str, float] = {}
        called: dict[str, int] = {}
        valid = True
        for label in labels:
            values = [
                raw[str(read["read_id"])][column_index]
                for read in group_reads[label]
                if raw[str(read["read_id"])][column_index] is not None
            ]
            minimum = max(2, min(5, math.ceil(len(group_reads[label]) * 0.15)))
            if len(values) < minimum:
                valid = False
                break
            means[label] = sum(values) / len(values)
            called[label] = len(values)
        if valid:
            score = max(means.values()) - min(means.values())
            column_scores.append((score, column_index, means, called))

    if not column_scores:
        raise AssertionError(f"No displayable CpGs for {record['dataset']} {record['chrom']}:{record['pos']}")

    selected = sorted(column_scores, key=lambda item: (-item[0], item[1]))[:40]
    selected_indices = sorted(item[1] for item in selected)
    score_by_index = {item[1]: item[0] for item in selected}

    rows: list[dict[str, Any]] = []
    for read in sorted(
        core_reads,
        key=lambda item: (
            labels.index(item["label"]),
            int(item["read_id"]) if str(item["read_id"]).isdigit() else str(item["read_id"]),
        ),
    ):
        values = raw[str(read["read_id"])]
        rows.append(
            {
                "read_id": str(read["read_id"]),
                "read_name_short": str(read["read_name"])[:8],
                "group": f"MG-{read['label']}",
                "raw_group": read["label"],
                "hp": str(read["latest_hp"]),
                "ps": read["latest_ps"],
                "strand": read["strand"],
                "values": [
                    None if values[index] is None else round(float(values[index]), 4)
                    for index in selected_indices
                ],
            }
        )

    group_summaries: list[dict[str, Any]] = []
    for label in labels:
        reads = group_reads[label]
        all_values = [
            value
            for read in reads
            for value in raw[str(read["read_id"])]
            if value is not None
        ]
        selected_values = [
            raw[str(read["read_id"])][index]
            for read in reads
            for index in selected_indices
            if raw[str(read["read_id"])][index] is not None
        ]
        group_summaries.append(
            {
                "group": f"MG-{label}",
                "raw_group": label,
                "n_reads": len(reads),
                "all_cpg_mean_5mc": round(sum(all_values) / len(all_values), 4),
                "selected_cpg_mean_5mc": round(
                    sum(selected_values) / len(selected_values), 4
                ),
                "called_cells_all": len(all_values),
            }
        )

    return {
        "source": str(matrix_path),
        "source_sha256": observed_sha,
        "raw_read_rows": len(raw),
        "raw_cpg_columns": len(positions),
        "core_read_rows": len(rows),
        "groups": group_summaries,
        "display_contract": (
            "all stable focal-ALT core reads × top 40 CpGs ranked by between-MG "
            "mean range after per-group coverage gate; display-only selection"
        ),
        "display_cpg_positions": [positions[index] for index in selected_indices],
        "display_cpg_scores": [round(score_by_index[index], 4) for index in selected_indices],
        "rows": rows,
    }


def strict_chrom_dir(sample: str, chrom: str) -> Path:
    if sample == "HCC1395":
        return STRICT_ROOT / "hcc1395_strict_regions_v2/chromosomes" / chrom
    return (
        STRICT_ROOT
        / "all7_production_v1/samples"
        / sample
        / "strict_regions_v1/chromosomes"
        / chrom
    )


def strict_w_evidence(sample: str, chrom: str, focal_pos: int, partner_pos: int) -> list[dict[str, Any]]:
    directory = strict_chrom_dir(sample, chrom)
    prefix = directory / f"{sample}.{chrom}"
    edges_path = Path(str(prefix) + ".endpoint_edges.tsv.gz")
    members_path = Path(str(prefix) + ".site_component_membership.tsv.gz")
    components_path = Path(str(prefix) + ".components.tsv.gz")

    edge_rows: list[dict[str, str]] = []
    with gzip.open(edges_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            positions = {int(row["pos_i1"]), int(row["pos_j1"])}
            if positions == {focal_pos, partner_pos} and as_bool(row["passes_primary_threshold"]):
                edge_rows.append(row)

    members: dict[tuple[str, str, str, int], str] = {}
    with gzip.open(members_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if int(row["threshold"]) != 3 or int(row["pos1"]) not in {focal_pos, partner_pos}:
                continue
            key = (
                row["linkage_basis"],
                row["phase_set"],
                row["pos1"],
                int(row["threshold"]),
            )
            members[key] = row["component_id"]

    component_by_id: dict[str, dict[str, str]] = {}
    with gzip.open(components_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if int(row["threshold"]) == 3:
                component_by_id[row["component_id"]] = row

    output: list[dict[str, Any]] = []
    for edge in edge_rows:
        base = edge["linkage_basis"]
        phase_set = edge["phase_set"]
        focal_component = members.get((base, phase_set, str(focal_pos), 3))
        partner_component = members.get((base, phase_set, str(partner_pos), 3))
        if focal_component is None or focal_component != partner_component:
            continue
        component = component_by_id[focal_component]
        output.append(
            {
                "phase_set": int(phase_set),
                "hp": int(edge["hp_family"]),
                "linkage_basis": base,
                "k": int(component["k"]),
                "span_bp": int(component["span_bp"]),
                "direct_support": int(edge["support_total"]),
                "support_coordinate_order": {
                    "RR": int(edge["support_RR"]),
                    "RA": int(edge["support_RA"]),
                    "AR": int(edge["support_AR"]),
                    "AA": int(edge["support_AA"]),
                },
                "component_id": focal_component,
            }
        )
    if not output:
        raise AssertionError(
            f"No threshold-3 direct W evidence: {sample} {chrom}:{focal_pos}/{partner_pos}"
        )
    output.sort(key=lambda item: (item["phase_set"], item["hp"]))
    return output


def effect_label(groups: list[dict[str, Any]]) -> str:
    lowest = min(groups, key=lambda item: item["alt_fraction"])
    highest = max(groups, key=lambda item: item["alt_fraction"])
    return (
        f"{highest['group']} 的 linked ALT 比例為 {highest['alt_fraction']:.1%}，"
        f"{lowest['group']} 為 {lowest['alt_fraction']:.1%}，相差 "
        f"{highest['alt_fraction'] - lowest['alt_fraction']:.1%}。"
    )


def endpoint_b_reason(status: str) -> str:
    if status.endswith("INSUFFICIENT_FOUR_STATE"):
        return "四種聯合 allele state 的資訊不足，無法排除多種關係。"
    if status.endswith("NO_FOCAL_REF"):
        return "缺少 focal-REF 分子，祖先違反率無法穩定估計。"
    if status.endswith("FIXED_ERROR_CEILING"):
        return "觀察到的違反狀態超過固定錯誤模型可安全解釋的範圍。"
    return "目前證據不足以識別唯一突變順序。"


def build_case(
    row: dict[str, str],
    assignment: dict[str, Any],
) -> dict[str, Any]:
    groups_raw = parse_json(row["endpoint_a_groups"], [])
    table = parse_json(row["endpoint_a_table"], [])
    if len(groups_raw) != len(table):
        raise AssertionError("Endpoint-A group/table length mismatch")
    groups: list[dict[str, Any]] = []
    for raw_group, counts in zip(groups_raw, table):
        ref_count, alt_count = [int(value) for value in counts]
        total = ref_count + alt_count
        groups.append(
            {
                "group": f"MG-{raw_group}",
                "raw_group": raw_group,
                "partner_ref_n": ref_count,
                "partner_alt_n": alt_count,
                "callable_n": total,
                "alt_fraction": alt_count / total if total else 0.0,
            }
        )

    focal_pos = int(row["focal_pos"])
    partner_pos = int(row["partner_pos"])
    state_counts = parse_json(row["endpoint_b_state_counts"], {})
    core_calls = parse_json(row["core_partner_call_counts"], {})
    methyl = matrix_for_record(assignment)
    strict_w = strict_w_evidence(
        row["sample"], row["focal_chrom"], focal_pos, partner_pos
    )
    endpoint_status = row["endpoint_b_status"]

    return {
        "id": f"{row['sample']}-{row['focal_chrom']}-{focal_pos}",
        "dataset": row["sample"],
        "biological_id": row["biological_id"],
        "evidence_level": "G1",
        "focal": {
            "chrom": row["focal_chrom"],
            "pos": focal_pos,
            "ref": row["focal_ref"],
            "alt": row["focal_alt"],
            "truth": row["focal_truth_label"],
        },
        "partner": {
            "chrom": row["partner_chrom"],
            "pos": partner_pos,
            "ref": row["partner_ref"],
            "alt": row["partner_alt"],
            "truth": row["partner_truth_label"],
        },
        "distance_bp": int(row["distance_bp"]),
        "core_focal_alt_reads": int(row["n_core_focal_alt_reads"]),
        "all_focal_ref_alt_reads": int(row["n_all_focal_ref_alt_reads"]),
        "partner_calls_in_core": {
            "ref": int(core_calls.get("R", 0)),
            "alt": int(core_calls.get("A", 0)),
            "noncallable": int(core_calls.get("X", 0)),
        },
        "groups": groups,
        "statistics": {
            "exact_p": as_float(row["endpoint_a_p_fixed_margins_exact"]),
            "q_global_by": float(row["endpoint_a_q_global_by"]),
            "cramers_v": float(row["endpoint_a_cramers_v"]),
            "delta_alt_fraction": float(row["endpoint_a_delta_alt_fraction"]),
            "conditional_p": float(row["endpoint_a_p_conditional_perm"]),
            "permutations": int(row["endpoint_a_permutations"]),
            "fdr_family_n": 147,
            "callability_q_global_by": as_float(row["callability_q_global_by"]),
            "callability_cramers_v": as_float(row["callability_cramers_v"]),
        },
        "callability": {
            "status": row["callability_gate_status"],
            "pass": as_bool(row["callability_gate_pass"]),
            "table": parse_json(row["callability_table"], []),
            "noncallable_core_reads": int(row["callability_noncallable_core_reads"]),
        },
        "association_statement": effect_label(groups),
        "state_graph": {
            "order": "focal-first",
            "RR": int(state_counts.get("RR", 0)),
            "RA": int(state_counts.get("RA", 0)),
            "AR": int(state_counts.get("AR", 0)),
            "AA": int(state_counts.get("AA", 0)),
            "X": int(state_counts.get("X", 0)),
            "O": int(state_counts.get("O", 0)),
        },
        "endpoint_b": {
            "status": endpoint_status,
            "reason_zh": endpoint_b_reason(endpoint_status),
            "compatible_relation_models": parse_json(
                row["endpoint_b_compatible_relation_models"], []
            ),
            "n_compatible_relation_models": int(
                row["endpoint_b_n_compatible_relation_models"]
            ),
            "claim_guardrail": row["endpoint_b_claim_guardrail"],
        },
        "legacy_posthoc_topology": {
            "scope": row["topology_scope"],
            "region": row["topology_region"],
            "order_status": row["topology_order_status"],
            "role": "posthoc context only; not current directed topology evidence",
        },
        "cross_platform": {
            "status": row["cross_platform_replication_status"],
            "exact_pair_present": as_bool(row["cross_platform_exact_pair_present"]),
            "biological_n": as_int(row["cross_platform_biological_n"]),
            "claim": "未達跨平台正式重現；不可視為獨立生物重複。",
        },
        "strict_w": strict_w,
        "methylation": methyl,
        "claim": {
            "supported": (
                "在 focal-ALT molecules 內，5mC pattern 所定義的 MG 與 linked partner "
                "sSNV allele state 顯著共分離，且兩位點位於同一 exact-PS×HP W 並有直接 read support。"
            ),
            "not_supported": (
                "不能由此確認甲基造成突變、突變造成甲基、cellular clone/subclone、"
                "父子節點或唯一演化樹。"
            ),
        },
    }


def summed_funnel() -> dict[str, int]:
    keys = next(iter(SAMPLE_FUNNEL.values())).keys()
    return {key: sum(values[key] for values in SAMPLE_FUNNEL.values()) for key in keys}


def build_report_data() -> dict[str, Any]:
    summary = json.loads(SUMMARY.read_text())
    formal_rows, pair_census = load_pair_authority()
    if len(formal_rows) != 7:
        raise AssertionError(f"Expected 7 formal G1 pairs, observed {len(formal_rows)}")

    target_keys = {
        (row["sample"], row["focal_chrom"], int(row["focal_pos"])) for row in formal_rows
    }
    assignments, assignment_counts = load_assignments(target_keys)
    cases = [
        build_case(
            row,
            assignments[(row["sample"], row["focal_chrom"], int(row["focal_pos"]))],
        )
        for row in formal_rows
    ]

    funnel = summed_funnel()
    expected_global = {
        "all_sites": 469849,
        "m1": int(summary["n_stable_sites"]),
        "m2_evaluable": int(summary["global_exact_fdr_family"]["n_m2_evaluable_sites"]),
        "m2_pass": int(summary["global_exact_fdr_family"]["n_m2_eligible_sites"]),
        "exact_pairs": int(summary["global_exact_fdr_family"]["n_exact_family_pairs"]),
        "by": int(summary["n_endpoint_a_exact_by_discoveries"]),
        "g1": int(summary["n_pair_by_confirmed"]),
    }
    for key, expected in expected_global.items():
        if funnel[key] != expected:
            raise AssertionError(f"Funnel mismatch for {key}: {funnel[key]} != {expected}")
    if sum(assignment_counts.values()) != funnel["m1"]:
        raise AssertionError("Assignment-line M1 census does not match funnel")
    for sample in SAMPLES:
        if assignment_counts[sample] != SAMPLE_FUNNEL[sample]["m1"]:
            raise AssertionError(f"M1 sample mismatch: {sample}")
        if pair_census[sample]["exact_pairs"] != SAMPLE_FUNNEL[sample]["exact_pairs"]:
            raise AssertionError(f"Exact-pair sample mismatch: {sample}")
        if pair_census[sample]["by"] != SAMPLE_FUNNEL[sample]["by"]:
            raise AssertionError(f"BY sample mismatch: {sample}")
        if pair_census[sample]["g1"] != SAMPLE_FUNNEL[sample]["g1"]:
            raise AssertionError(f"G1 sample mismatch: {sample}")

    fixed_states = {
        key: sum(case["state_graph"][key] for case in cases)
        for key in ("RR", "RA", "AR", "AA", "X")
    }
    strict_containers = sum(len(case["strict_w"]) for case in cases)
    strict_support = sum(
        evidence["direct_support"] for case in cases for evidence in case["strict_w"]
    )
    if strict_containers != 10 or strict_support != 790:
        raise AssertionError("Strict-W reconciliation failed")

    dataset_rows = []
    for sample in SAMPLES:
        row = {"dataset": sample, **SAMPLE_FUNNEL[sample], **pair_census[sample]}
        row["m1_rate"] = row["m1"] / row["all_sites"]
        row["g1_exact_rate"] = (
            row["g1"] / row["exact_pairs"] if row["exact_pairs"] else None
        )
        dataset_rows.append(row)

    return {
        "schema_name": "intersubmod.methyl_linked_ssnv_g1_evidence_atlas",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_type": "B_comprehensive_validation",
        "scope": {
            "technical_datasets": 7,
            "biological_samples": 6,
            "whole_genome": True,
            "partial": False,
            "note": "HCC1395 and HCC1395_DORADO count as one biological sample.",
        },
        "release_status": {
            "producer_execution_integrity": summary["status"],
            "semantics": summary["pass_semantics"],
            "scientific_release": "NOT_SIGNED_TASK_B_FINAL_RELEASE",
        },
        "headline": {
            "formal_g1_pairs": 7,
            "formal_g1_focal_sites": 7,
            "datasets_with_g1": 3,
            "biological_samples_with_g1": 3,
            "g1_per_exact_pair": 7 / 147,
            "g2_global_by": 0,
            "identified_mutation_orders": 0,
            "strict_w_pair_coverage": 7,
            "strict_w_containers": strict_containers,
            "strict_direct_support": strict_support,
            "fixed_state_totals": fixed_states,
        },
        "funnel": {
            **funnel,
            "m1_rate_all": funnel["m1"] / funnel["all_sites"],
            "m2_pass_rate_all": funnel["m2_pass"] / funnel["all_sites"],
            "m2_pass_rate_evaluable": funnel["m2_pass"] / funnel["m2_evaluable"],
            "g1_rate_exact": funnel["g1"] / funnel["exact_pairs"],
            "g2_family_sites": int(
                summary["joint_signature_global_fdr_audit"]["n_family_sites"]
            ),
        },
        "dataset_rows": dataset_rows,
        "cases": cases,
        "definitions": {
            "MG": "由 focal-ALT reads 的 5mC pattern 形成的穩定 methyl group；不是 HP tag。",
            "G1": (
                "M2 residual methyl partition + linked-sSNV fixed-margins global BY + "
                "V≥0.30 + ΔALT fraction≥0.50 + conditional HP-family×PS×strand + callability PASS。"
            ),
            "W": (
                "在 exact phase-set × HP family 容器中，以 direct read support≥3 建立的"
                "無向連結 component。"
            ),
            "tree": (
                "本報告只顯示 observed read-state graph 與 W component；"
                "7/7 都無法識別唯一 mutation order，因此不是已確認演化樹。"
            ),
        },
        "sources": [
            {
                "id": "pair-v8",
                "path": str(PAIR_RESULTS),
                "sha256": sha256_file(PAIR_RESULTS),
                "role": "G1 pair statistics and focal-first state counts",
            },
            {
                "id": "summary-v8",
                "path": str(SUMMARY),
                "sha256": sha256_file(SUMMARY),
                "role": "whole-scope denominators and execution status",
            },
            {
                "id": "assignment-v10",
                "path": str(ASSIGNMENTS),
                "sha256": sha256_file(ASSIGNMENTS),
                "role": "stable focal-ALT core-read MG labels and matrix identity",
            },
            {
                "id": "claim-contract-v5",
                "path": str(CLAIM_CONTRACT),
                "sha256": sha256_file(CLAIM_CONTRACT),
                "role": "M1/M2/G1/G2/L1 claim ceiling",
            },
            {
                "id": "strict-w-join-audit",
                "path": str(STRICT_JOIN_AUDIT),
                "sha256": sha256_file(STRICT_JOIN_AUDIT),
                "role": "exact-PS×HP W linkage reconciliation",
            },
        ],
    }


HTML_TEMPLATE = r"""<!doctype html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light">
<title>甲基群 × linked-sSNV：G1 全樣本關聯圖譜</title>
<style>
:root{
  --ink:#102a32;--muted:#5d6e72;--paper:#f5f1e8;--panel:#fffdf8;--line:#c9c3b7;
  --teal:#006d77;--cyan:#65bfc3;--orange:#e76f51;--amber:#e9b949;--navy:#17324d;
  --green:#4a7c59;--red:#a23e48;--gray:#9ba6a8;--group1:#006d77;--group2:#d67b3f;
  --group3:#7a5ea8;--group4:#4a7c59;--shadow:5px 5px 0 rgba(16,42,50,.10);
}
*{box-sizing:border-box}
html{scroll-behavior:smooth;background:var(--paper)}
body{margin:0;color:var(--ink);background:var(--paper);font-family:Inter,ui-sans-serif,system-ui,-apple-system,"Noto Sans TC",sans-serif;line-height:1.55}
a{color:var(--teal)}
.shell{max-width:1440px;margin:auto;padding:0 28px 80px}
.topbar{position:sticky;top:0;z-index:20;background:rgba(245,241,232,.96);border-bottom:1px solid var(--line);backdrop-filter:blur(10px)}
.topbar-inner{max-width:1440px;margin:auto;padding:10px 28px;display:flex;align-items:center;gap:12px;overflow:auto}
.brand{font-family:ui-serif,Georgia,serif;font-weight:800;white-space:nowrap;margin-right:auto}
.navlink{font-size:12px;text-decoration:none;color:var(--muted);white-space:nowrap;padding:5px 8px;border:1px solid transparent}
.navlink:hover{border-color:var(--line);color:var(--ink)}
.hero{padding:66px 0 42px;border-bottom:2px solid var(--ink);display:grid;grid-template-columns:1.25fr .75fr;gap:50px;align-items:end}
.kicker{font:700 12px/1.3 ui-monospace,SFMono-Regular,monospace;letter-spacing:.12em;color:var(--orange);text-transform:uppercase}
h1{font-family:ui-serif,Georgia,"Noto Serif TC",serif;font-size:clamp(40px,6vw,82px);line-height:.98;letter-spacing:-.045em;margin:16px 0 24px;max-width:960px}
.hero p{font-size:18px;max-width:820px;margin:0;color:#31464b}
.status-ribbon{display:inline-flex;gap:8px;flex-wrap:wrap;margin-top:25px}
.chip{display:inline-flex;align-items:center;gap:6px;border:1px solid currentColor;padding:5px 9px;font:700 11px/1.2 ui-monospace,SFMono-Regular,monospace;letter-spacing:.02em;background:var(--panel)}
.chip.ok{color:var(--green)}.chip.warn{color:#9a5a00}.chip.stop{color:var(--red)}.chip.info{color:var(--teal)}
.hero-score{border-left:1px solid var(--line);padding-left:34px}
.score-big{font:800 clamp(80px,12vw,160px)/.8 ui-serif,Georgia,serif;color:var(--teal);letter-spacing:-.08em}
.score-label{font:750 16px/1.35 ui-monospace,SFMono-Regular,monospace;margin-top:18px}
.score-sub{font-size:13px;color:var(--muted);margin-top:8px}
.section{padding:58px 0;border-bottom:1px solid var(--line)}
.section-head{display:grid;grid-template-columns:190px 1fr;gap:28px;margin-bottom:30px}
.section-no{font:700 13px/1 ui-monospace,SFMono-Regular,monospace;color:var(--orange)}
h2{font:800 clamp(30px,4vw,50px)/1.05 ui-serif,Georgia,"Noto Serif TC",serif;margin:0;letter-spacing:-.03em}
.lede{font-size:17px;color:#334a4f;max-width:900px;margin:12px 0 0}
.answer-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:18px}
.answer-card{background:var(--panel);border:1px solid var(--ink);padding:24px;box-shadow:var(--shadow);min-height:210px}
.answer-card .n{font:800 54px/.9 ui-serif,Georgia,serif;color:var(--teal)}
.answer-card h3{font-size:16px;margin:17px 0 8px}
.answer-card p{font-size:13px;color:var(--muted);margin:0}
.answer-card.danger{border-color:var(--red)}.answer-card.danger .n{color:var(--red)}
.funnel{display:grid;grid-template-columns:repeat(6,1fr);gap:2px;border:1px solid var(--ink);background:var(--ink)}
.funnel-step{background:var(--panel);padding:18px 12px;min-height:132px;position:relative}
.funnel-step strong{display:block;font:800 29px/1 ui-serif,Georgia,serif}
.funnel-step span{display:block;font-size:12px;margin-top:9px;color:var(--muted)}
.funnel-step em{display:block;font:700 10px/1.3 ui-monospace,SFMono-Regular,monospace;margin-top:12px;color:var(--teal);font-style:normal}
.table-wrap{overflow:auto;border:1px solid var(--line);background:var(--panel)}
table{width:100%;border-collapse:collapse;font-size:13px}
th,td{padding:12px 13px;border-bottom:1px solid #ddd7cc;text-align:right;white-space:nowrap}
th:first-child,td:first-child{text-align:left;position:sticky;left:0;background:var(--panel);z-index:1}
thead th{font:700 11px/1.25 ui-monospace,SFMono-Regular,monospace;color:var(--muted);background:#ece7dc;position:sticky;top:0}
tbody tr:hover td{background:#f1eee6}tbody tr:hover td:first-child{background:#f1eee6}
.zero{color:#a0a5a4}.hit{color:var(--teal);font-weight:800}
.evidence-ladder{display:grid;grid-template-columns:repeat(5,1fr);gap:12px}
.rung{border-top:5px solid var(--gray);background:var(--panel);padding:17px;min-height:180px}
.rung.active{border-color:var(--teal)}.rung.future{border-style:dashed}
.rung b{font:800 21px/1 ui-monospace,SFMono-Regular,monospace}.rung p{font-size:12px;color:var(--muted)}
.rung .rung-count{font:800 32px/1 ui-serif,Georgia,serif;color:var(--teal)}
.filter-row{display:flex;gap:8px;flex-wrap:wrap;margin:0 0 25px}
.filter{border:1px solid var(--ink);background:var(--panel);color:var(--ink);font:700 12px/1 ui-monospace,SFMono-Regular,monospace;padding:9px 12px;cursor:pointer}
.filter.active,.filter:hover{background:var(--ink);color:var(--paper)}
.case{background:var(--panel);border:1px solid var(--ink);margin:0 0 38px;box-shadow:var(--shadow)}
.case[hidden]{display:none}
.case-head{display:grid;grid-template-columns:1fr auto;border-bottom:1px solid var(--ink)}
.case-title{padding:23px 26px}
.case-index{font:700 11px/1 ui-monospace,SFMono-Regular,monospace;color:var(--orange);letter-spacing:.12em}
.case h3{font:800 clamp(21px,3vw,34px)/1.15 ui-serif,Georgia,serif;margin:8px 0 5px;letter-spacing:-.02em}
.locus-sub{font:650 13px/1.4 ui-monospace,SFMono-Regular,monospace;color:var(--muted)}
.case-verdict{min-width:210px;border-left:1px solid var(--ink);padding:22px;background:#edf4f2}
.case-verdict b{display:block;color:var(--teal);font-size:14px}.case-verdict span{font-size:11px;color:var(--muted)}
.locus-track{height:110px;margin:24px 26px 0;position:relative}
.track-line{position:absolute;left:4%;right:4%;top:51px;height:2px;background:var(--ink)}
.track-node{position:absolute;top:36px;width:30px;height:30px;border-radius:50%;transform:translateX(-50%);background:var(--panel);border:5px solid var(--teal);z-index:2}
.track-node.partner{border-color:var(--orange)}
.track-label{position:absolute;top:0;transform:translateX(-50%);font:700 10px/1.35 ui-monospace,SFMono-Regular,monospace;text-align:center;white-space:nowrap}
.track-distance{position:absolute;left:50%;top:65px;transform:translateX(-50%);font:700 11px/1 ui-monospace,SFMono-Regular,monospace;color:var(--muted);background:var(--panel);padding:4px 7px}
.case-grid{display:grid;grid-template-columns:1.08fr .92fr;border-top:1px solid var(--line)}
.pane{padding:25px 26px;border-right:1px solid var(--line)}.pane:nth-child(even){border-right:0}
.pane.full{grid-column:1/-1;border-right:0;border-top:1px solid var(--line)}
.pane h4{font:800 15px/1.25 ui-monospace,SFMono-Regular,monospace;margin:0 0 6px}
.pane-note{font-size:12px;color:var(--muted);margin:0 0 18px}
.metric-strip{display:grid;grid-template-columns:repeat(4,1fr);border:1px solid var(--line);margin-bottom:20px}
.metric{padding:11px;border-right:1px solid var(--line)}.metric:last-child{border:0}
.metric b{display:block;font:800 18px/1.1 ui-serif,Georgia,serif}.metric span{font-size:9px;color:var(--muted);font-family:ui-monospace,SFMono-Regular,monospace}
.assoc-row{display:grid;grid-template-columns:100px minmax(130px,1fr) 120px;gap:10px;align-items:center;margin:11px 0}
.group-name{font:700 11px/1 ui-monospace,SFMono-Regular,monospace}
.stack{height:22px;display:flex;border:1px solid var(--ink);overflow:hidden;background:#e8e5de}
.stack .ref{background:var(--cyan)}.stack .alt{background:var(--orange)}
.assoc-count{font:650 11px/1.35 ui-monospace,SFMono-Regular,monospace;text-align:right}
.statement{border-left:4px solid var(--orange);padding:10px 13px;background:#f6ebe4;font-size:13px;margin-top:18px}
.heatmap-wrap{position:relative;border:1px solid var(--line);background:#f0ede6;min-height:240px}
.heatmap{display:block;width:100%;height:280px}
.heat-tooltip{display:none;position:absolute;pointer-events:none;background:var(--ink);color:#fff;padding:7px 9px;font:11px/1.4 ui-monospace,SFMono-Regular,monospace;z-index:5;max-width:280px}
.legend{display:flex;align-items:center;gap:7px;font-size:10px;color:var(--muted);margin-top:8px}
.legend-ramp{display:grid;grid-template-columns:repeat(5,18px);gap:1px}.legend-ramp i{height:9px;display:block}
.group-summary{display:flex;gap:8px;flex-wrap:wrap;margin-top:11px}
.group-pill{border-left:4px solid var(--group1);padding:4px 7px;background:#f0ede6;font:650 10px/1.4 ui-monospace,SFMono-Regular,monospace}
.state-layout{display:grid;grid-template-columns:1fr .9fr;gap:22px}
.state-matrix{display:grid;grid-template-columns:74px 1fr 1fr;grid-template-rows:34px 1fr 1fr;min-height:215px}
.axis-corner{font:650 9px/1.3 ui-monospace,SFMono-Regular,monospace;color:var(--muted)}
.axis-head{font:700 10px/1.2 ui-monospace,SFMono-Regular,monospace;text-align:center;padding:8px}
.axis-side{font:700 10px/1.2 ui-monospace,SFMono-Regular,monospace;padding:25px 6px}
.state-cell{border:1px solid var(--line);display:flex;align-items:center;justify-content:center;flex-direction:column;background:#f0ede6}
.state-cell.alt{background:#fae8df}.state-cell b{font:800 30px/1 ui-serif,Georgia,serif}.state-cell span{font:700 10px/1 ui-monospace,SFMono-Regular,monospace;color:var(--muted)}
.topology-box{border:1px solid var(--line);padding:13px;background:#f7f4ed}
.topology-box svg{width:100%;height:82px;display:block}
.topology-status{font:750 11px/1.4 ui-monospace,SFMono-Regular,monospace;color:var(--red)}
.xbox{margin-top:9px;border:1px dashed var(--gray);padding:8px;font-size:11px;color:var(--muted)}
.w-list{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:10px}
.w-card{border:1px solid var(--teal);padding:12px;background:#eef6f4}
.w-card b{font:800 13px/1.3 ui-monospace,SFMono-Regular,monospace}.w-card p{font-size:11px;margin:5px 0 0;color:var(--muted)}
.claim-grid{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-top:18px}
.claim{padding:14px;border:1px solid var(--green);background:#eef4ea;font-size:12px}.claim.no{border-color:var(--red);background:#f8eaea}
.audit{margin-top:14px}.audit summary{cursor:pointer;font:700 11px/1.3 ui-monospace,SFMono-Regular,monospace;color:var(--muted)}
.audit pre{white-space:pre-wrap;font:10px/1.5 ui-monospace,SFMono-Regular,monospace;background:#efebe2;padding:12px;overflow:auto}
.method-grid{display:grid;grid-template-columns:repeat(2,1fr);gap:18px}
.method-card{border-top:4px solid var(--teal);background:var(--panel);padding:20px}
.method-card h3{font-size:16px;margin:0 0 8px}.method-card p,.method-card li{font-size:13px;color:var(--muted)}
.source-list{display:grid;gap:8px}.source{padding:11px;border:1px solid var(--line);background:var(--panel);font:10px/1.5 ui-monospace,SFMono-Regular,monospace;overflow-wrap:anywhere}
.footer{padding:35px 0;font-size:12px;color:var(--muted);display:flex;justify-content:space-between;gap:20px}
@media(max-width:900px){
  .hero,.section-head,.case-grid,.state-layout{grid-template-columns:1fr}
  .hero-score{border-left:0;border-top:1px solid var(--line);padding:25px 0 0}
  .answer-grid,.evidence-ladder{grid-template-columns:1fr 1fr}
  .funnel{grid-template-columns:repeat(3,1fr)}
  .case-head{grid-template-columns:1fr}.case-verdict{border-left:0;border-top:1px solid var(--ink)}
  .pane{border-right:0;border-bottom:1px solid var(--line)}.pane.full{border-bottom:0}
}
@media(max-width:560px){
  .shell{padding:0 14px 50px}.topbar-inner{padding:9px 14px}.hero{padding-top:42px}
  .answer-grid,.evidence-ladder,.method-grid{grid-template-columns:1fr}.funnel{grid-template-columns:1fr 1fr}
  .metric-strip{grid-template-columns:1fr 1fr}.metric:nth-child(2){border-right:0}
  .assoc-row{grid-template-columns:78px 1fr}.assoc-count{grid-column:2}
  .claim-grid{grid-template-columns:1fr}.case-title,.pane{padding:20px 17px}.locus-track{margin-left:17px;margin-right:17px}
}
@media print{
  .topbar,.filter-row{display:none}.shell{max-width:none;padding:0}.case{break-inside:avoid;box-shadow:none}
  body{background:white}.section{padding:28px 0}.heatmap{height:220px}
}
</style>
</head>
<body>
<nav class="topbar"><div class="topbar-inner">
  <div class="brand">InterSubMod · G1 evidence atlas</div>
  <a class="navlink" href="#answer">重點答案</a>
  <a class="navlink" href="#scope">全樣本分母</a>
  <a class="navlink" href="#cases">7 個位置</a>
  <a class="navlink" href="#methods">定義與限制</a>
</div></nav>
<main class="shell">
  <header class="hero">
    <div>
      <div class="kicker">Whole-genome · 7 technical datasets · 6 biological samples</div>
      <h1>甲基群 × linked-sSNV<br>全樣本 G1 關聯圖譜</h1>
      <p>回答「哪一個 linked 位點與甲基分群最有關係」：正式可確認的是同一批 focal-ALT molecules 中，methyl group（MG）與 partner sSNV allele 的共分離；不是因果、clone 身分或唯一演化樹。</p>
      <div class="status-ribbon">
        <span class="chip ok">FULL SCOPE</span>
        <span class="chip info">7 FORMAL G1 PAIRS</span>
        <span class="chip stop">0 IDENTIFIED ORDERS</span>
        <span class="chip warn">PRODUCER INTERMEDIATE</span>
      </div>
    </div>
    <div class="hero-score">
      <div class="score-big">7</div>
      <div class="score-label">個 linked-sSNV pair 通過完整 G1 gate</div>
      <div class="score-sub">分布於 HCC1395 ×1、H2009 ×5、HCC1954 ×1；其他四個 dataset 為 0。</div>
    </div>
  </header>

  <section class="section" id="answer">
    <div class="section-head"><div class="section-no">01 / ANSWER</div><div>
      <h2>可以確認「哪個 linked 位點和甲基群有關」；不能確認誰造成誰</h2>
      <p class="lede">判讀時要把 association、read linkage、evolutionary order 三層分開。這是本報告最重要的敘述邊界。</p>
    </div></div>
    <div class="answer-grid">
      <article class="answer-card"><div class="n">7</div><h3>正式 methyl–allele 共分離</h3><p>7/147 exact-family pairs（4.76%）同時通過 global BY、效果量、conditional permutation 與 callability gate。</p></article>
      <article class="answer-card"><div class="n">10</div><h3>exact-PS×HP W containers</h3><p>7/7 pair 都位於同一 W 且有直接無向 read support；部分 pair 同時出現在 HP1 與 HP2，因此共 10 個 W。</p></article>
      <article class="answer-card danger"><div class="n">0</div><h3>可確認的突變先後</h3><p>Endpoint-B 全部不可識別。下方只畫 observed state graph，不畫已確認的父子箭頭或真實演化樹。</p></article>
    </div>
  </section>

  <section class="section" id="scope">
    <div class="section-head"><div class="section-no">02 / SCOPE</div><div>
      <h2>完整樣本漏斗：訊號存在，但通過嚴格 G1 的比例小</h2>
      <p class="lede">分母由全 469,849 個 sSNV 開始；G1 是嚴格、低產率的局部 association evidence，不代表其餘位置沒有甲基訊號。</p>
    </div></div>
    <div class="funnel" id="funnel"></div>
    <div style="height:22px"></div>
    <div class="table-wrap"><table id="dataset-table"><thead><tr>
      <th>Dataset</th><th>全 sSNV</th><th>M1 stable</th><th>M2 evaluable</th><th>M2 PASS</th><th>有 partner</th><th>Exact focal</th><th>Exact pairs</th><th>Global BY</th><th>Formal G1</th>
    </tr></thead><tbody></tbody></table></div>
  </section>

  <section class="section">
    <div class="section-head"><div class="section-no">03 / CLAIM LADDER</div><div>
      <h2>證據階梯：G1 停在局部共分離，不自動升級成 clone</h2>
    </div></div>
    <div class="evidence-ladder">
      <div class="rung active"><b>M1</b><div class="rung-count">102,842</div><p>focal-ALT reads 內有穩定 5mC 多群結構。</p></div>
      <div class="rung active"><b>M2</b><div class="rung-count">919</div><p>八個已量測技術／讀段軸不能充分解釋分群。</p></div>
      <div class="rung active"><b>G1</b><div class="rung-count">7</div><p>MG 與一個 linked partner allele 通過嚴格共分離標準。</p></div>
      <div class="rung future"><b>G2</b><div class="rung-count">0</div><p>需要多 marker joint signature；目前 global BY 為 0/58。</p></div>
      <div class="rung future"><b>L1</b><div class="rung-count">—</div><p>細胞 clone 需 single-cell、colony、spatial 或獨立多區域證據。</p></div>
    </div>
  </section>

  <section class="section" id="cases">
    <div class="section-head"><div class="section-no">04 / LOCUS ATLAS</div><div>
      <h2>所有 7 個正式位置與可視證據</h2>
      <p class="lede">每張卡依序顯示 locus、MG×partner allele、全 stable-core read×CpG 熱圖、focal-first state graph、exact-PS×HP W 與 claim ceiling。</p>
    </div></div>
    <div class="table-wrap" style="margin-bottom:22px"><table id="case-index"><thead><tr>
      <th>強度序</th><th>Dataset</th><th>Focal</th><th>Linked partner</th><th>距離</th><th>BY q</th><th>Cramér's V</th><th>ΔALT</th><th>W</th><th>Order</th>
    </tr></thead><tbody></tbody></table></div>
    <p class="pane-note" style="margin-bottom:22px">上表依 Cramér's V 由高至低；「關聯最強」是統計效果量的描述，不表示因果更強或演化先後更確定。</p>
    <div class="filter-row" id="filters"></div>
    <div id="case-list"></div>
  </section>

  <section class="section" id="methods">
    <div class="section-head"><div class="section-no">05 / METHODS</div><div>
      <h2>如何閱讀數字、熱圖與「樹」</h2>
    </div></div>
    <div class="method-grid">
      <article class="method-card"><h3>MG 是什麼</h3><p>MG-1-1、MG-1-2 等是 focal-ALT reads 依 5mC pattern 形成的穩定甲基群標籤，<strong>不是 HP1/HP2</strong>。HP 另用於 conditional test 與 exact-PS×HP W。</p></article>
      <article class="method-card"><h3>G1 確認什麼</h3><p>在 focal-ALT molecules 內，MG 與 linked partner 的 REF/ALT 比例不同，且通過全域 BY、V≥0.30、ΔALT≥0.50、conditional permutation 與 callability。</p></article>
      <article class="method-card"><h3>熱圖顯示什麼</h3><p>每例包含全部 stable focal-ALT core reads；欄位為通過 per-group coverage gate 後，依群間平均差排序的前 40 個 CpG。這是差異 pattern 的可視化，不是額外的 discovery test。</p></article>
      <article class="method-card"><h3>為何不是演化樹</h3><p>W 線段只代表兩位點有直接 read linkage。2×2 state counts 全以 focal-first 解讀；目前 7/7 Endpoint-B 都不能識別唯一 mutation order，legacy posthoc 欄位也不得升級成箭頭。</p></article>
      <article class="method-card"><h3>5mC 範圍</h3><p>矩陣呈現 5mC probability（C+m）；不混入 5hmC。灰色代表該 read 的該 CpG 未呼叫。</p></article>
      <article class="method-card"><h3>跨平台限制</h3><p>HCC1395 與 DORADO 是同一 biological sample 的不同來源，這 7 個 G1 沒有正式跨平台重現；因此本報告支持局部訊號偵測，不支持獨立生物重複。</p></article>
    </div>
    <div style="height:28px"></div>
    <h3 style="font-family:ui-monospace,SFMono-Regular,monospace">SOURCE AUTHORITY</h3>
    <div class="source-list" id="sources"></div>
  </section>
  <footer class="footer"><span>InterSubMod · Task B comprehensive validation · G3/G4</span><span id="build-time"></span></footer>
</main>
<script id="report-data" type="application/json">__REPORT_DATA__</script>
<script>
const DATA=JSON.parse(document.getElementById('report-data').textContent);
const n=v=>Number(v).toLocaleString('en-US');
const pct=(v,d=1)=>v==null?'—':`${(100*v).toFixed(d)}%`;
const qfmt=v=>v<0.001?v.toExponential(2):v.toFixed(4);
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));

function renderFunnel(){
  const f=DATA.funnel;
  const steps=[
    [f.all_sites,'全 sSNV','100%'],
    [f.m1,'M1 stable',pct(f.m1_rate_all,2)],
    [f.m2_evaluable,'M2 evaluable',pct(f.m2_evaluable/f.all_sites,3)],
    [f.m2_pass,'M2 PASS',pct(f.m2_pass_rate_evaluable,1)+' of evaluable'],
    [f.exact_pairs,'Exact-family pairs','128 focal sites'],
    [f.g1,'Formal G1',pct(f.g1_rate_exact,2)+' of exact pairs']
  ];
  document.getElementById('funnel').innerHTML=steps.map(([v,l,r])=>
    `<div class="funnel-step"><strong>${n(v)}</strong><span>${l}</span><em>${r}</em></div>`).join('');
}
function renderDatasetTable(){
  document.querySelector('#dataset-table tbody').innerHTML=DATA.dataset_rows.map(r=>{
    const cells=[r.all_sites,r.m1,r.m2_evaluable,r.m2_pass,r.m2_with_partner,r.exact_focal_sites,r.exact_pairs,r.by,r.g1];
    return `<tr><td><b>${esc(r.dataset)}</b></td>${cells.map((v,i)=>`<td class="${v===0?'zero':i===8?'hit':''}">${n(v)}</td>`).join('')}</tr>`;
  }).join('');
}
function renderCaseIndex(){
  const ranked=[...DATA.cases].sort((a,b)=>b.statistics.cramers_v-a.statistics.cramers_v);
  document.querySelector('#case-index tbody').innerHTML=ranked.map((c,i)=>`<tr>
    <td>${i+1}</td><td><b>${esc(c.dataset)}</b></td>
    <td>${c.focal.chrom}:${n(c.focal.pos)} ${c.focal.ref}&gt;${c.focal.alt}</td>
    <td>${n(c.partner.pos)} ${c.partner.ref}&gt;${c.partner.alt}${c.partner.truth==='FP'?' · <span style="color:var(--red);font-weight:800">FP</span>':''}</td>
    <td>${c.distance_bp>0?'+':''}${n(c.distance_bp)} bp</td>
    <td>${qfmt(c.statistics.q_global_by)}</td><td class="hit">${c.statistics.cramers_v.toFixed(3)}</td>
    <td>${pct(c.statistics.delta_alt_fraction,1)}</td><td>${c.strict_w.length} container${c.strict_w.length>1?'s':''}</td>
    <td style="color:var(--red)">未定</td></tr>`).join('');
}
function track(c){
  const min=Math.min(c.focal.pos,c.partner.pos),max=Math.max(c.focal.pos,c.partner.pos),span=Math.max(1,max-min);
  const x=p=>8+84*(p-min)/span;
  return `<div class="locus-track">
    <div class="track-line"></div>
    <div class="track-label" style="left:${x(c.focal.pos)}%">Focal<br>${n(c.focal.pos)} ${c.focal.ref}&gt;${c.focal.alt}</div>
    <div class="track-node" style="left:${x(c.focal.pos)}%"></div>
    <div class="track-label" style="left:${x(c.partner.pos)}%;top:78px">Partner<br>${n(c.partner.pos)} ${c.partner.ref}&gt;${c.partner.alt}</div>
    <div class="track-node partner" style="left:${x(c.partner.pos)}%"></div>
    <div class="track-distance">${c.distance_bp>0?'+':''}${n(c.distance_bp)} bp · 無向 linkage</div>
  </div>`;
}
function association(c){
  return c.groups.map((g,i)=>{
    const ref=100*(1-g.alt_fraction),alt=100*g.alt_fraction;
    return `<div class="assoc-row">
      <div class="group-name" style="color:var(--group${i+1})">${g.group}</div>
      <div class="stack" aria-label="${g.group}: partner REF ${g.partner_ref_n}, ALT ${g.partner_alt_n}">
        <div class="ref" style="width:${ref}%"></div><div class="alt" style="width:${alt}%"></div>
      </div>
      <div class="assoc-count">R ${g.partner_ref_n} · A ${g.partner_alt_n}<br>ALT ${pct(g.alt_fraction,1)}</div>
    </div>`;
  }).join('');
}
function topologySvg(){
  return `<svg viewBox="0 0 280 82" role="img" aria-label="三個候選關係都未確認">
    <g fill="none" stroke="#9ba6a8" stroke-width="2" stroke-dasharray="4 4">
      <path d="M18 65V35H58V15M58 35H98V15"/><path d="M110 65V35H150V15M150 35H190V15"/>
      <path d="M202 65V38M202 38L232 16M202 38L262 16"/>
    </g>
    <g fill="#fffdf8" stroke="#9ba6a8" stroke-width="2">
      <circle cx="58" cy="15" r="7"/><circle cx="98" cy="15" r="7"/>
      <circle cx="150" cy="15" r="7"/><circle cx="190" cy="15" r="7"/>
      <circle cx="232" cy="16" r="7"/><circle cx="262" cy="16" r="7"/>
    </g>
    <g fill="#5d6e72" font-size="9" font-family="monospace">
      <text x="20" y="80">F→P ?</text><text x="112" y="80">P→F ?</text><text x="205" y="80">sister ?</text>
    </g>
  </svg>`;
}
function stateGraph(c){
  const s=c.state_graph;
  return `<div class="state-layout">
    <div>
      <div class="state-matrix">
        <div class="axis-corner">first=focal<br>second=partner</div>
        <div class="axis-head">partner REF</div><div class="axis-head">partner ALT</div>
        <div class="axis-side">focal REF</div>
        <div class="state-cell"><b>${s.RR}</b><span>RR</span></div><div class="state-cell alt"><b>${s.RA}</b><span>RA</span></div>
        <div class="axis-side">focal ALT</div>
        <div class="state-cell"><b>${s.AR}</b><span>AR</span></div><div class="state-cell alt"><b>${s.AA}</b><span>AA</span></div>
      </div>
      <div class="xbox">X=${s.X}：partner 不可判讀；不進入 direct-edge fixed-state support。</div>
    </div>
    <div class="topology-box">${topologySvg()}
      <div class="topology-status">MUTATION ORDER：NOT IDENTIFIABLE</div>
      <p class="pane-note" style="margin-top:7px">${esc(c.endpoint_b.reason_zh)}</p>
      <p class="pane-note">compatible relation models = ${c.endpoint_b.n_compatible_relation_models}</p>
    </div>
  </div>`;
}
function wCards(c){
  return c.strict_w.map(w=>`<div class="w-card">
    <b>PS ${n(w.phase_set)} · HP${w.hp} · W k=${w.k}</b>
    <p>direct support=${n(w.direct_support)} · span=${n(w.span_bp)} bp · threshold≥3</p>
  </div>`).join('');
}
function caseCard(c,i){
  const m=c.statistics;
  const truthWarn=c.partner.truth==='FP'?'<span class="chip stop">PARTNER FP</span>':'<span class="chip ok">TP / TP</span>';
  const legacy=c.legacy_posthoc_topology;
  return `<article class="case" data-dataset="${esc(c.dataset)}">
    <div class="case-head">
      <div class="case-title"><div class="case-index">G1 LOCUS ${String(i+1).padStart(2,'0')} · ${esc(c.dataset)}</div>
        <h3>${c.focal.chrom}:${n(c.focal.pos)} ${c.focal.ref}&gt;${c.focal.alt} ↔ ${n(c.partner.pos)} ${c.partner.ref}&gt;${c.partner.alt}</h3>
        <div class="locus-sub">${truthWarn} <span class="chip info">CALLABILITY PASS</span> <span class="chip stop">ORDER UNRESOLVED</span></div>
      </div>
      <div class="case-verdict"><b>MG ↔ linked allele：確認</b><span>因果／clone／演化順序：未確認</span></div>
    </div>
    ${track(c)}
    <div class="case-grid">
      <section class="pane">
        <h4>A · METHYL GROUP × PARTNER ALLELE</h4>
        <p class="pane-note">每列為 [partner REF, ALT]；MG 是甲基群，不是 HP。</p>
        <div class="metric-strip">
          <div class="metric"><b>${qfmt(m.q_global_by)}</b><span>GLOBAL BY q</span></div>
          <div class="metric"><b>${m.cramers_v.toFixed(3)}</b><span>CRAMÉR'S V</span></div>
          <div class="metric"><b>${pct(m.delta_alt_fraction,1)}</b><span>Δ ALT FRACTION</span></div>
          <div class="metric"><b>${m.conditional_p.toFixed(3)}</b><span>CONDITIONAL p</span></div>
        </div>
        ${association(c)}
        <div class="legend"><span style="width:13px;height:9px;background:var(--cyan)"></span>partner REF <span style="width:13px;height:9px;background:var(--orange)"></span>partner ALT</div>
        <div class="statement">${esc(c.association_statement)}</div>
      </section>
      <section class="pane">
        <h4>B · 5mC READ × CpG PATTERN</h4>
        <p class="pane-note">全部 ${n(c.methylation.core_read_rows)} 條 stable focal-ALT core reads；顯示 ${c.methylation.display_cpg_positions.length}/${n(c.methylation.raw_cpg_columns)} CpGs。</p>
        <div class="heatmap-wrap"><canvas class="heatmap" data-case="${esc(c.id)}"></canvas><div class="heat-tooltip"></div></div>
        <div class="legend"><span>0</span><span class="legend-ramp"><i style="background:#17324d"></i><i style="background:#397d8a"></i><i style="background:#e8e3d8"></i><i style="background:#e9b949"></i><i style="background:#e76f51"></i></span><span>1 · 5mC probability</span><span style="margin-left:8px">灰=NA</span></div>
        <div class="group-summary">${c.methylation.groups.map((g,j)=>`<span class="group-pill" style="border-color:var(--group${j+1})">${g.group}: n=${g.n_reads}, all-CpG mean=${g.all_cpg_mean_5mc.toFixed(3)}</span>`).join('')}</div>
      </section>
      <section class="pane">
        <h4>C · OBSERVED MOLECULAR-STATE GRAPH</h4>
        <p class="pane-note">非演化樹；固定採 focal-first，避免基因組座標方向造成 RA/AR 顛倒。</p>
        ${stateGraph(c)}
      </section>
      <section class="pane">
        <h4>D · EXACT-PS×HP W LINKAGE</h4>
        <p class="pane-note">兩位點同一 W 且 direct edge 通過 read-support threshold；線是無向。</p>
        <div class="w-list">${wCards(c)}</div>
        <div class="claim-grid">
          <div class="claim"><b>可支持</b><br>${esc(c.claim.supported)}</div>
          <div class="claim no"><b>不可支持</b><br>${esc(c.claim.not_supported)}</div>
        </div>
        <details class="audit"><summary>展開 legacy/posthoc 與來源稽核欄位</summary><pre>${esc(JSON.stringify({legacy_posthoc_topology:legacy,cross_platform:c.cross_platform,endpoint_b_status:c.endpoint_b.status,matrix_sha256:c.methylation.source_sha256},null,2))}</pre></details>
      </section>
    </div>
  </article>`;
}
function color(v){
  if(v==null)return '#b9bebc';
  const stops=[[0,[23,50,77]],[.25,[57,125,138]],[.5,[232,227,216]],[.75,[233,185,73]],[1,[231,111,81]]];
  for(let i=1;i<stops.length;i++){if(v<=stops[i][0]){const [x0,c0]=stops[i-1],[x1,c1]=stops[i],t=(v-x0)/(x1-x0);return `rgb(${c0.map((c,j)=>Math.round(c+t*(c1[j]-c))).join(',')})`;}}
  return 'rgb(231,111,81)';
}
function drawHeatmap(canvas,c){
  const wrap=canvas.parentElement,tooltip=wrap.querySelector('.heat-tooltip');
  const cssW=Math.max(300,wrap.clientWidth),cssH=Math.max(235,Math.min(360,110+c.methylation.rows.length*1.35));
  const dpr=window.devicePixelRatio||1;canvas.width=cssW*dpr;canvas.height=cssH*dpr;canvas.style.height=cssH+'px';
  const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);ctx.clearRect(0,0,cssW,cssH);
  const rows=c.methylation.rows,cols=c.methylation.display_cpg_positions,left=62,right=8,top=24,bottom=26;
  const cw=(cssW-left-right)/cols.length,rh=(cssH-top-bottom)/rows.length;
  const groups=[...new Set(rows.map(r=>r.group))];
  rows.forEach((r,ri)=>{
    ctx.fillStyle=getComputedStyle(document.documentElement).getPropertyValue(`--group${groups.indexOf(r.group)+1}`).trim()||'#006d77';
    ctx.fillRect(0,top+ri*rh,8,Math.max(1,rh));
    r.values.forEach((v,ci)=>{ctx.fillStyle=color(v);ctx.fillRect(left+ci*cw,top+ri*rh,Math.ceil(cw+.2),Math.ceil(rh+.2));});
  });
  ctx.strokeStyle='#102a32';ctx.lineWidth=1;ctx.strokeRect(left,top,cssW-left-right,cssH-top-bottom);
  ctx.fillStyle='#5d6e72';ctx.font='10px monospace';ctx.fillText(n(cols[0]),left,13);
  const last=n(cols[cols.length-1]);ctx.fillText(last,cssW-right-ctx.measureText(last).width,13);
  let start=0;groups.forEach((g,gi)=>{const end=rows.findLastIndex(r=>r.group===g);const y=top+(start+end+1)*rh/2;ctx.fillStyle=getComputedStyle(document.documentElement).getPropertyValue(`--group${gi+1}`).trim();ctx.font='bold 9px monospace';ctx.fillText(g.replace('MG-',''),12,y);ctx.strokeStyle='#fff';ctx.beginPath();ctx.moveTo(left,top+(end+1)*rh);ctx.lineTo(cssW-right,top+(end+1)*rh);ctx.stroke();start=end+1;});
  ctx.fillStyle='#5d6e72';ctx.font='9px monospace';ctx.fillText(`${rows.length} reads · all rows`,left,cssH-7);
  canvas.onmousemove=e=>{const rect=canvas.getBoundingClientRect(),x=e.clientX-rect.left,y=e.clientY-rect.top,ci=Math.floor((x-left)/cw),ri=Math.floor((y-top)/rh);if(ci>=0&&ci<cols.length&&ri>=0&&ri<rows.length){const r=rows[ri],v=r.values[ci];tooltip.style.display='block';tooltip.style.left=Math.min(x+12,cssW-230)+'px';tooltip.style.top=Math.max(4,y-45)+'px';tooltip.textContent=`${r.group} · read ${r.read_name_short} · CpG ${cols[ci]} · 5mC ${v==null?'NA':v.toFixed(3)} · HP ${r.hp}`;}else tooltip.style.display='none';};
  canvas.onmouseleave=()=>tooltip.style.display='none';
}
function renderCases(){
  const list=document.getElementById('case-list');list.innerHTML=DATA.cases.map(caseCard).join('');
  requestAnimationFrame(()=>document.querySelectorAll('.heatmap').forEach(cv=>drawHeatmap(cv,DATA.cases.find(c=>c.id===cv.dataset.case))));
}
function renderFilters(){
  const datasets=['全部',...new Set(DATA.cases.map(c=>c.dataset))];
  const box=document.getElementById('filters');
  box.innerHTML=datasets.map((d,i)=>`<button class="filter ${i===0?'active':''}" data-filter="${d}">${d}${d==='全部'?' · 7':` · ${DATA.cases.filter(c=>c.dataset===d).length}`}</button>`).join('');
  box.onclick=e=>{if(!e.target.matches('.filter'))return;box.querySelectorAll('.filter').forEach(b=>b.classList.remove('active'));e.target.classList.add('active');const f=e.target.dataset.filter;document.querySelectorAll('.case').forEach(card=>card.hidden=f!=='全部'&&card.dataset.dataset!==f);};
}
function renderSources(){
  document.getElementById('sources').innerHTML=DATA.sources.map(s=>`<div class="source"><b>${esc(s.id)}</b> · ${esc(s.role)}<br>${esc(s.path)}<br>SHA-256 ${esc(s.sha256)}</div>`).join('');
}
renderFunnel();renderDatasetTable();renderCaseIndex();renderCases();renderFilters();renderSources();
document.getElementById('build-time').textContent='Built '+DATA.created_at_utc;
let rt;window.addEventListener('resize',()=>{clearTimeout(rt);rt=setTimeout(()=>document.querySelectorAll('.heatmap').forEach(cv=>drawHeatmap(cv,DATA.cases.find(c=>c.id===cv.dataset.case))),120)});
</script>
</body></html>
"""


def main() -> None:
    args = parse_args()
    args.output_html.parent.mkdir(parents=True, exist_ok=True)
    args.output_data.parent.mkdir(parents=True, exist_ok=True)
    args.output_receipt.parent.mkdir(parents=True, exist_ok=True)

    data = build_report_data()
    args.output_data.write_text(
        json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace(
        "</", "<\\/"
    )
    args.output_html.write_text(
        HTML_TEMPLATE.replace("__REPORT_DATA__", payload), encoding="utf-8"
    )
    receipt = {
        "schema_name": "intersubmod.g1_html_build_receipt",
        "schema_version": "1.0.0",
        "created_at_utc": data["created_at_utc"],
        "status": "PASS",
        "inputs": {
            "pair_results": str(PAIR_RESULTS),
            "summary": str(SUMMARY),
            "assignments": str(ASSIGNMENTS),
            "strict_root": str(STRICT_ROOT),
        },
        "outputs": {
            "html": {
                "path": str(args.output_html),
                "sha256": sha256_file(args.output_html),
                "size_bytes": args.output_html.stat().st_size,
            },
            "data": {
                "path": str(args.output_data),
                "sha256": sha256_file(args.output_data),
                "size_bytes": args.output_data.stat().st_size,
            },
        },
        "assertions": {
            "formal_g1_pairs": len(data["cases"]),
            "all_cases_have_methyl_matrix": all(
                case["methylation"]["core_read_rows"] > 0 for case in data["cases"]
            ),
            "all_cases_have_strict_w": all(case["strict_w"] for case in data["cases"]),
            "identified_mutation_orders": data["headline"]["identified_mutation_orders"],
            "whole_scope": data["scope"]["whole_genome"],
        },
    }
    args.output_receipt.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": "PASS",
                "html": str(args.output_html),
                "data": str(args.output_data),
                "receipt": str(args.output_receipt),
                "cases": len(data["cases"]),
                "core_reads": sum(
                    case["methylation"]["core_read_rows"] for case in data["cases"]
                ),
                "strict_w_containers": data["headline"]["strict_w_containers"],
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
