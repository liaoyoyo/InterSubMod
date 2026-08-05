#!/usr/bin/env python3
"""Build the canonical portable-report artifact for cross-sample topology matrices."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import sqlite3
import statistics
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SAMPLES = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
BIOLOGICAL_IDS = ("HCC1395", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954")
BIOLOGICAL_GROUPS = {
    "HCC1395": "breast",
    "HCC1937": "breast",
    "HCC1954": "breast",
    "H1437": "lung_adenocarcinoma",
    "H2009": "lung_adenocarcinoma",
    "COLO829": "melanoma",
}
SAMPLE_META = (
    {
        "dataset_order": 1,
        "dataset": "HCC1395",
        "biological_id": "HCC1395",
        "organ": "breast",
        "histology": "invasive ductal carcinoma",
        "technical_role": "canonical technical source",
    },
    {
        "dataset_order": 2,
        "dataset": "HCC1395_DORADO",
        "biological_id": "HCC1395",
        "organ": "breast",
        "histology": "invasive ductal carcinoma",
        "technical_role": "same biological ID; alternate technical source",
    },
    {
        "dataset_order": 3,
        "dataset": "COLO829",
        "biological_id": "COLO829",
        "organ": "skin / melanocytic",
        "histology": "cutaneous melanoma",
        "technical_role": "independent biological ID",
    },
    {
        "dataset_order": 4,
        "dataset": "H1437",
        "biological_id": "H1437",
        "organ": "lung",
        "histology": "NSCLC adenocarcinoma",
        "technical_role": "independent biological ID",
    },
    {
        "dataset_order": 5,
        "dataset": "H2009",
        "biological_id": "H2009",
        "organ": "lung",
        "histology": "adenocarcinoma",
        "technical_role": "independent biological ID",
    },
    {
        "dataset_order": 6,
        "dataset": "HCC1937",
        "biological_id": "HCC1937",
        "organ": "breast",
        "histology": "primary ductal carcinoma",
        "technical_role": "independent biological ID",
    },
    {
        "dataset_order": 7,
        "dataset": "HCC1954",
        "biological_id": "HCC1954",
        "organ": "breast",
        "histology": "primary invasive ductal carcinoma",
        "technical_role": "independent biological ID",
    },
)
PAIR_CLASS_LABELS = {
    "same_biological_id": "同 biological ID",
    "same_cancer_different_id": "同癌種、不同 ID",
    "cross_cancer": "跨癌種",
}
TITLE = "七資料集 sSNV／區域拓撲比對矩陣"
REPORT_DB_PATH = (
    "research/20260725_crosssample_topology_matrix_report/data/report.sqlite"
)


class BuildError(RuntimeError):
    """Raised when an artifact input or output contract fails."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise BuildError(message)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument("--cohort-similarity-json", required=True, type=Path)
    parser.add_argument("--artifact-json", required=True, type=Path)
    parser.add_argument("--sqlite-db", required=True, type=Path)
    return parser.parse_args()


def parse_scalar(value: str) -> Any:
    text = value.strip()
    if text == "":
        return None
    if text == "true":
        return True
    if text == "false":
        return False
    try:
        if all(character not in text.lower() for character in (".", "e")):
            return int(text)
        return float(text)
    except ValueError:
        return text


def read_tsv(path: Path) -> list[dict[str, Any]]:
    require(path.is_file(), f"missing TSV: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [
            {key: parse_scalar(value) for key, value in row.items()}
            for row in csv.DictReader(handle, delimiter="\t")
        ]
    require(rows, f"empty TSV: {path}")
    return rows


def read_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing JSON: {path}")
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    require(isinstance(payload, dict), f"expected JSON object: {path}")
    return payload


def read_matrix(path: Path, *, blank_diagonal: bool = True) -> list[dict[str, Any]]:
    rows = read_tsv(path)
    require(tuple(row["dataset"] for row in rows) == SAMPLES, f"matrix order drift: {path}")
    for row_order, row in enumerate(rows, 1):
        row["row_order"] = row_order
        if blank_diagonal:
            row[row["dataset"]] = None
    require(
        all(tuple(column for column in row if column not in {"dataset", "row_order"}) == SAMPLES for row in rows),
        f"matrix columns drift: {path}",
    )
    return rows


def biological_pair(left: str, right: str) -> tuple[str, str]:
    order = {sample: index for index, sample in enumerate(BIOLOGICAL_IDS)}
    return (left, right) if order[left] < order[right] else (right, left)


def sample_biological_id(sample: str) -> str:
    return "HCC1395" if sample == "HCC1395_DORADO" else sample


def sample_group(sample: str) -> str:
    return BIOLOGICAL_GROUPS[sample_biological_id(sample)]


def pair_class(left: str, right: str) -> str:
    left_id = sample_biological_id(left)
    right_id = sample_biological_id(right)
    if left_id == right_id:
        return "same_biological_id"
    if BIOLOGICAL_GROUPS[left_id] == BIOLOGICAL_GROUPS[right_id]:
        return "same_cancer_different_id"
    return "cross_cancer"


def build_composite_profile_rows(
    cohort: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, float]]:
    matrix = cohort["matrix"]
    require(tuple(matrix["sample_order"]) == SAMPLES, "cohort profile sample order drift")
    wide_rows: list[dict[str, Any]] = []
    for index, sample in enumerate(SAMPLES):
        row: dict[str, Any] = {"dataset": sample, "row_order": index + 1}
        for peer_index, peer in enumerate(SAMPLES):
            row[peer] = None if peer == sample else matrix["similarity"][index][peer_index]
        wide_rows.append(row)

    pair_values: dict[tuple[str, str], list[float]] = {}
    technical_rows: list[dict[str, Any]] = []
    for row in cohort["pairwise_profiles"]:
        left = row["sample_a"]
        right = row["sample_b"]
        similarity = float(row["profile_similarity"])
        technical_rows.append(
            {
                "pair": f"{left} × {right}",
                "dataset_a": left,
                "dataset_b": right,
                "pair_class": pair_class(left, right),
                "pair_class_label": PAIR_CLASS_LABELS[pair_class(left, right)],
                "profile_similarity": similarity,
                "rank_of_21": int(row["similarity_rank_of_21"]),
            }
        )
        left_id = sample_biological_id(left)
        right_id = sample_biological_id(right)
        if left_id != right_id:
            pair_values.setdefault(biological_pair(left_id, right_id), []).append(similarity)

    macro = {key: statistics.fmean(values) for key, values in pair_values.items()}
    require(len(macro) == 15, "composite profile biological pair count != 15")
    class_means: dict[str, float] = {}
    for class_name in ("same_cancer_different_id", "cross_cancer"):
        values = [
            value
            for (left, right), value in macro.items()
            if (
                "same_cancer_different_id"
                if BIOLOGICAL_GROUPS[left] == BIOLOGICAL_GROUPS[right]
                else "cross_cancer"
            )
            == class_name
        ]
        class_means[class_name] = statistics.fmean(values)
    target = next(
        row["profile_similarity"]
        for row in technical_rows
        if {row["dataset_a"], row["dataset_b"]} == {"HCC1395", "HCC1395_DORADO"}
    )
    class_means["same_biological_id"] = target
    return wide_rows, technical_rows, class_means


def exact_partition_test(
    biological_pair_values: Mapping[tuple[str, str], float],
) -> dict[str, Any]:
    require(len(biological_pair_values) == 15, "exact partition test needs 15 pairs")
    observed_groups = {
        "breast": {"HCC1395", "HCC1937", "HCC1954"},
        "lung": {"H1437", "H2009"},
        "melanoma": {"COLO829"},
    }

    def statistic(groups: Iterable[set[str]]) -> tuple[float, float, float]:
        groups = tuple(groups)
        membership = {
            biological_id: index
            for index, group in enumerate(groups)
            for biological_id in group
        }
        within: list[float] = []
        between: list[float] = []
        for pair, value in biological_pair_values.items():
            (within if membership[pair[0]] == membership[pair[1]] else between).append(value)
        require(len(within) == 4 and len(between) == 11, "3/2/1 denominator drift")
        within_mean = statistics.fmean(within)
        between_mean = statistics.fmean(between)
        return within_mean - between_mean, within_mean, between_mean

    observed, within_mean, between_mean = statistic(observed_groups.values())
    assignments: list[float] = []
    ids = set(BIOLOGICAL_IDS)
    for size_three in itertools.combinations(BIOLOGICAL_IDS, 3):
        remaining = ids.difference(size_three)
        for size_two in itertools.combinations(sorted(remaining), 2):
            size_one = remaining.difference(size_two)
            assignments.append(statistic((set(size_three), set(size_two), size_one))[0])
    require(len(assignments) == 60, "exact assignment count != 60")
    greater_or_equal = sum(value >= observed - 1e-15 for value in assignments)
    return {
        "observed_within_mean": within_mean,
        "observed_between_mean": between_mean,
        "observed_within_minus_between": observed,
        "permutation_n": 60,
        "greater_or_equal_n": greater_or_equal,
        "one_sided_exact_p": greater_or_equal / 60,
    }


def infer_sql_type(values: Iterable[Any]) -> str:
    non_null = [value for value in values if value is not None]
    if non_null and all(isinstance(value, (bool, int)) for value in non_null):
        return "INTEGER"
    if non_null and all(isinstance(value, (bool, int, float)) for value in non_null):
        return "REAL"
    return "TEXT"


def quote_identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


def create_table(
    connection: sqlite3.Connection, table: str, rows: Sequence[Mapping[str, Any]]
) -> None:
    require(rows, f"empty table: {table}")
    columns: list[str] = []
    for row in rows:
        for column in row:
            if column not in columns:
                columns.append(column)
    column_types = {
        column: infer_sql_type(row.get(column) for row in rows) for column in columns
    }
    connection.execute(
        f"CREATE TABLE {quote_identifier(table)} ("
        + ", ".join(
            f"{quote_identifier(column)} {column_types[column]}" for column in columns
        )
        + ")"
    )
    placeholders = ", ".join("?" for _ in columns)
    connection.executemany(
        f"INSERT INTO {quote_identifier(table)} VALUES ({placeholders})",
        [
            tuple(
                int(row.get(column)) if isinstance(row.get(column), bool) else row.get(column)
                for column in columns
            )
            for row in rows
        ],
    )


def execute_query(connection: sqlite3.Connection, sql: str) -> list[dict[str, Any]]:
    cursor = connection.execute(sql)
    columns = [description[0] for description in cursor.description]
    return [dict(zip(columns, row)) for row in cursor.fetchall()]


def source_record(
    source_id: str,
    label: str,
    sql: str,
    tables: Sequence[str],
    metric_definitions: Sequence[str],
    executed_at: str,
) -> dict[str, Any]:
    return {
        "id": source_id,
        "label": label,
        "path": REPORT_DB_PATH,
        "query": {
            "engine": "sqlite",
            "language": "sql",
            "sql": sql,
            "description": f"Executed and reviewed SQLite query for {label}.",
            "tables_used": list(tables),
            "filters": ["GRCh38 chr1–22", "7 technical datasets / 6 biological IDs"],
            "metric_definitions": list(metric_definitions),
            "executed_at": executed_at,
        },
    }


def heatmap_chart(
    chart_id: str,
    dataset: str,
    title: str,
    subtitle: str,
    source: Mapping[str, Any],
    *,
    value_format: str = "percent",
) -> dict[str, Any]:
    return {
        "id": chart_id,
        "dataset": dataset,
        "type": "heatmap",
        "title": title,
        "subtitle": subtitle,
        "intent": "relationship",
        "question": "七個 technical datasets 的兩兩相似度如何分布？",
        "rationale": "A common 7×7 ordering makes same-ID, same-cancer, and cross-cancer patterns directly comparable across evidence layers.",
        "comparisonContext": {
            "grain": "technical-dataset pair",
            "denominator": "metric-specific pair universe",
            "diagonal": "N/A to avoid self-comparison dominating the scale",
        },
        "encodings": {
            "x": {"field": "dataset", "type": "nominal", "label": "Row dataset"},
            "y": {
                "fields": list(SAMPLES),
                "labels": {"HCC1395_DORADO": "HCC1395 DORADO"},
            },
        },
        "valueFormat": value_format,
        "palette": {"kind": "sequential"},
        "sourceId": source["id"],
        "source": deepcopy(source),
    }


def bar_chart(
    chart_id: str,
    dataset: str,
    title: str,
    subtitle: str,
    source: Mapping[str, Any],
    x_field: str,
    y_field: str,
    color_field: str,
    *,
    value_format: str = "percent",
) -> dict[str, Any]:
    return {
        "id": chart_id,
        "dataset": dataset,
        "type": "bar",
        "title": title,
        "subtitle": subtitle,
        "intent": "compare",
        "question": "配對類型或證據層之間的相似度差異有多大？",
        "rationale": "Grouped bars preserve a common scale and keep pair classes visible.",
        "encodings": {
            "x": {"field": x_field, "type": "nominal", "label": x_field},
            "y": {"field": y_field, "type": "quantitative", "label": "Similarity"},
            "color": {
                "field": color_field,
                "type": "nominal",
                "label": color_field,
            },
        },
        "options": {"orientation": "vertical", "grouping": "grouped"},
        "valueFormat": value_format,
        "sourceId": source["id"],
        "source": deepcopy(source),
    }


def table_definition(
    table_id: str,
    dataset: str,
    title: str,
    subtitle: str,
    source: Mapping[str, Any],
    columns: Sequence[Mapping[str, Any]],
    sort_field: str,
    direction: str = "asc",
) -> dict[str, Any]:
    return {
        "id": table_id,
        "dataset": dataset,
        "title": title,
        "subtitle": subtitle,
        "density": "compact",
        "layout": "full",
        "defaultSort": {"field": sort_field, "direction": direction},
        "columns": list(columns),
        "sourceId": source["id"],
        "source": deepcopy(source),
    }


def main() -> int:
    args = parse_args()
    require(args.analysis_dir.is_dir(), f"missing analysis directory: {args.analysis_dir}")
    require(not args.artifact_json.exists(), f"artifact already exists: {args.artifact_json}")
    require(not args.sqlite_db.exists(), f"SQLite DB already exists: {args.sqlite_db}")

    analysis = read_json(args.analysis_dir / "analysis_summary.json")
    receipt = read_json(args.analysis_dir / "validation_receipt.json")
    cohort = read_json(args.cohort_similarity_json)
    require(receipt["checks"]["all_pass"] is True, "matrix receipt is not all-pass")
    require(all(receipt["checks"].values()), "one or more matrix checks failed")
    require(all(cohort["checks"].values()), "one or more cohort-profile checks failed")

    target = analysis["hcc1395_target"]
    pair_rows_raw = read_tsv(args.analysis_dir / "technical_pair_metrics.tsv")
    permutation_rows_raw = read_tsv(
        args.analysis_dir / "biological_id_exact_permutation_tests.tsv"
    )
    biological_rows = read_tsv(args.analysis_dir / "biological_pair_aggregates.tsv")
    composite_matrix, _, composite_class_means = build_composite_profile_rows(cohort)

    matrices = {
        "candidate_matrix": read_matrix(
            args.analysis_dir / "matrices/candidate_sites_jaccard.tsv"
        ),
        "active_matrix": read_matrix(
            args.analysis_dir / "matrices/active_sites_jaccard.tsv"
        ),
        "w_matrix": read_matrix(args.analysis_dir / "matrices/w_sites_jaccard.tsv"),
        "edge_matrix": read_matrix(
            args.analysis_dir / "matrices/primary_edges_jaccard.tsv"
        ),
        "component_matrix": read_matrix(
            args.analysis_dir / "matrices/exact_components_jaccard.tsv"
        ),
        "vaf_matrix": read_matrix(
            args.analysis_dir / "matrices/vaf_latest_js_similarity.tsv"
        ),
        "coarse_topology_matrix": read_matrix(
            args.analysis_dir / "matrices/coarse_topology_js_similarity.tsv"
        ),
        "exact_coordinate_n_matrix": read_matrix(
            args.analysis_dir / "matrices/exact_coordinate_one_to_one_n.tsv"
        ),
        "composite_profile_matrix": composite_matrix,
    }

    pair_rows: list[dict[str, Any]] = []
    pair_order = sorted(
        pair_rows_raw,
        key=lambda row: (-float(row["candidate_sites_jaccard"]), row["dataset_a"], row["dataset_b"]),
    )
    composite_lookup = {
        frozenset((row["sample_a"], row["sample_b"])): row["profile_similarity"]
        for row in cohort["pairwise_profiles"]
    }
    for order, row in enumerate(pair_order, 1):
        key = frozenset((row["dataset_a"], row["dataset_b"]))
        pair_rows.append(
            {
                "pair_order": order,
                "pair": f"{row['dataset_a']} × {row['dataset_b']}",
                "pair_class": PAIR_CLASS_LABELS[row["pair_class"]],
                "candidate_jaccard": row["candidate_sites_jaccard"],
                "active_jaccard": row["active_sites_jaccard"],
                "w_jaccard": row["w_sites_jaccard"],
                "primary_edge_jaccard": row["primary_edges_jaccard"],
                "exact_component_jaccard": row["exact_components_jaccard"],
                "vaf_js_similarity": row["vaf_latest_js_similarity"],
                "coarse_topology_js_similarity": row["coarse_topology_js_similarity"],
                "composite_profile_similarity": composite_lookup[key],
                "exact_coordinate_1to1_n": row["exact_coordinate_one_to_one_n"],
                "k_ge_2_shape_equal_rate": row["both_unique_k_ge_2_shape_equal_rate"],
                "ancestor_equal_rate": row["both_unique_ancestor_equal_rate"],
            }
        )

    global_ladder = [
        {"layer_order": 1, "layer": "Candidate sites", "similarity": target["candidate_sites_jaccard"]},
        {"layer_order": 2, "layer": "Active sites", "similarity": target["active_sites_jaccard"]},
        {"layer_order": 3, "layer": "W sites", "similarity": target["w_sites_jaccard"]},
        {"layer_order": 4, "layer": "Primary edges", "similarity": target["primary_edges_jaccard"]},
        {"layer_order": 5, "layer": "Exact components", "similarity": target["exact_components_jaccard"]},
        {"layer_order": 6, "layer": "VAF profile", "similarity": target["vaf_latest_js_similarity"]},
        {"layer_order": 7, "layer": "Coarse topology", "similarity": target["coarse_topology_js_similarity"]},
        {
            "layer_order": 8,
            "layer": "Composite topology profile",
            "similarity": cohort["hcc1395_dorado_chromosome_block_bootstrap"]["point_similarity"],
        },
    ]
    conditional_topology = [
        {
            "metric_order": 1,
            "metric": "Topology set compatible",
            "rate": target["topology_signature_set_compatible_rate"],
            "numerator": target["topology_signature_set_compatible_n"],
            "denominator": target["exact_coordinate_one_to_one_n"],
        },
        {
            "metric_order": 2,
            "metric": "Unique shape equal",
            "rate": target["both_unique_shape_equal_rate"],
            "numerator": target["both_unique_shape_equal_n"],
            "denominator": target["both_unique_n"],
        },
        {
            "metric_order": 3,
            "metric": "k≥2 unique shape equal",
            "rate": target["both_unique_k_ge_2_shape_equal_rate"],
            "numerator": target["both_unique_k_ge_2_shape_equal_n"],
            "denominator": target["both_unique_k_ge_2_n"],
        },
        {
            "metric_order": 4,
            "metric": "Ancestor relation equal",
            "rate": target["both_unique_ancestor_equal_rate"],
            "numerator": target["both_unique_ancestor_equal_n"],
            "denominator": target["both_unique_ancestor_evaluable_n"],
        },
        {
            "metric_order": 5,
            "metric": "DORADO edge contained in HCC",
            "rate": target["shared_candidate_primary_edges_dor_to_hcc"],
            "numerator": 11110,
            "denominator": 11305,
        },
    ]

    metric_labels = {
        "candidate_sites_jaccard": "Candidate-site Jaccard",
        "active_sites_jaccard": "Active-site Jaccard",
        "w_sites_jaccard": "W-site Jaccard",
        "primary_edges_jaccard": "Primary-edge Jaccard",
        "exact_components_jaccard": "Exact-component Jaccard",
        "coarse_topology_similarity": "Coarse topology TVD similarity",
        "coarse_topology_js_similarity": "Coarse topology JS similarity",
        "vaf_latest_js_similarity": "Latest VAF JS similarity",
        "vaf_baseline_js_similarity": "Baseline VAF JS similarity",
    }
    permutation_rows: list[dict[str, Any]] = []
    for order, row in enumerate(permutation_rows_raw, 1):
        permutation_rows.append(
            {
                "metric_order": order,
                "metric": metric_labels[row["metric"]],
                "within_mean": row["observed_within_mean"],
                "between_mean": row["observed_between_mean"],
                "within_minus_between": row["observed_within_minus_between"],
                "permutation_n": row["permutation_n"],
                "one_sided_exact_p": row["one_sided_exact_p"],
                "interpretation": (
                    "exploratory group signal"
                    if float(row["one_sided_exact_p"]) <= 0.05
                    else "no group-level support at 0.05"
                ),
            }
        )

    profile_macro_values: dict[tuple[str, str], list[float]] = {}
    for row in cohort["pairwise_profiles"]:
        left = sample_biological_id(row["sample_a"])
        right = sample_biological_id(row["sample_b"])
        if left == right:
            continue
        profile_macro_values.setdefault(biological_pair(left, right), []).append(
            float(row["profile_similarity"])
        )
    profile_test = exact_partition_test(
        {pair: statistics.fmean(values) for pair, values in profile_macro_values.items()}
    )
    permutation_rows.append(
        {
            "metric_order": len(permutation_rows) + 1,
            "metric": "Composite exact-PS profile similarity",
            "within_mean": profile_test["observed_within_mean"],
            "between_mean": profile_test["observed_between_mean"],
            "within_minus_between": profile_test["observed_within_minus_between"],
            "permutation_n": 60,
            "one_sided_exact_p": profile_test["one_sided_exact_p"],
            "interpretation": "exploratory group signal",
        }
    )

    labeled_class_rows: list[dict[str, Any]] = []
    profile_class_rows: list[dict[str, Any]] = []
    target_class_values = {
        "candidate_sites_jaccard": target["candidate_sites_jaccard"],
        "active_sites_jaccard": target["active_sites_jaccard"],
        "w_sites_jaccard": target["w_sites_jaccard"],
        "primary_edges_jaccard": target["primary_edges_jaccard"],
        "coarse_topology_js_similarity": target["coarse_topology_js_similarity"],
        "vaf_latest_js_similarity": target["vaf_latest_js_similarity"],
        "composite_profile_similarity": composite_class_means["same_biological_id"],
    }
    permutation_lookup = {row["metric"]: row for row in permutation_rows_raw}
    for class_name in ("same_biological_id", "same_cancer_different_id", "cross_cancer"):
        class_label = PAIR_CLASS_LABELS[class_name]
        for metric, layer in (
            ("candidate_sites_jaccard", "Candidate sites"),
            ("active_sites_jaccard", "Active sites"),
            ("w_sites_jaccard", "W sites"),
            ("primary_edges_jaccard", "Primary edges"),
        ):
            if class_name == "same_biological_id":
                value = target_class_values[metric]
            elif class_name == "same_cancer_different_id":
                value = permutation_lookup[metric]["observed_within_mean"]
            else:
                value = permutation_lookup[metric]["observed_between_mean"]
            labeled_class_rows.append(
                {"layer": layer, "pair_class": class_label, "similarity": value}
            )
        for metric, layer in (
            ("vaf_latest_js_similarity", "VAF profile"),
            ("coarse_topology_js_similarity", "Coarse topology"),
            ("composite_profile_similarity", "Composite topology"),
        ):
            if metric == "composite_profile_similarity":
                value = composite_class_means[class_name]
            elif class_name == "same_biological_id":
                value = target_class_values[metric]
            elif class_name == "same_cancer_different_id":
                value = permutation_lookup[metric]["observed_within_mean"]
            else:
                value = permutation_lookup[metric]["observed_between_mean"]
            profile_class_rows.append(
                {"layer": layer, "pair_class": class_label, "similarity": value}
            )

    bootstrap = cohort["hcc1395_dorado_chromosome_block_bootstrap"]
    robustness_rows = [
        {
            "order": 1,
            "check": "Point profile similarity",
            "estimate": bootstrap["point_similarity"],
            "lower": None,
            "upper": None,
            "note": "rank 1/21 at point estimate",
        },
        {
            "order": 2,
            "check": "22-autosome bootstrap similarity",
            "estimate": bootstrap["similarity"]["median"],
            "lower": bootstrap["similarity"]["p2_5"],
            "upper": bootstrap["similarity"]["p97_5"],
            "note": "2,000 chromosome-block replicates",
        },
        {
            "order": 3,
            "check": "Rank-1 rate",
            "estimate": bootstrap["rank_of_21"]["rank_1_rate"],
            "lower": None,
            "upper": None,
            "note": "top-3 rate = 0.986",
        },
        {
            "order": 4,
            "check": "Gap versus best other pair",
            "estimate": bootstrap["gap_vs_best_other_pair"]["median"],
            "lower": bootstrap["gap_vs_best_other_pair"]["p2_5"],
            "upper": bootstrap["gap_vs_best_other_pair"]["p97_5"],
            "note": "interval crosses zero; exact rank is chromosome-sensitive",
        },
    ]

    args.sqlite_db.parent.mkdir(parents=True, exist_ok=False)
    executed_at = datetime.now(timezone.utc).isoformat()
    datasets: dict[str, list[dict[str, Any]]] = {}
    sources: list[dict[str, Any]] = []
    connection = sqlite3.connect(args.sqlite_db)
    try:
        table_rows: dict[str, list[dict[str, Any]]] = {
            **matrices,
            "sample_metadata": list(SAMPLE_META),
            "global_ladder": global_ladder,
            "conditional_topology": conditional_topology,
            "labeled_class_summary": labeled_class_rows,
            "profile_class_summary": profile_class_rows,
            "permutation_tests": permutation_rows,
            "robustness": robustness_rows,
            "pairwise_metrics": pair_rows,
        }
        for table, rows in table_rows.items():
            create_table(connection, table, rows)
        connection.commit()

        query_specs: dict[str, tuple[str, list[str]]] = {}
        matrix_columns = ", ".join(
            [quote_identifier("dataset")]
            + [quote_identifier(sample) for sample in SAMPLES]
        )
        for table in matrices:
            query_specs[table] = (
                f"SELECT {matrix_columns} FROM {quote_identifier(table)} ORDER BY row_order;",
                [
                    "Off-diagonal cells are pairwise values; diagonal cells are N/A.",
                    "Values retain the metric-specific denominator and claim boundary.",
                ],
            )
        query_specs.update(
            {
                "sample_metadata": (
                    "SELECT dataset_order, dataset, biological_id, organ, histology, technical_role "
                    "FROM sample_metadata ORDER BY dataset_order;",
                    ["HCC1395_DORADO is a technical source of HCC1395, not an independent biological ID."],
                ),
                "global_ladder": (
                    "SELECT layer_order, layer, similarity FROM global_ladder ORDER BY layer_order;",
                    ["Set Jaccard and profile similarities are distinct measures on a shared 0–1 display scale."],
                ),
                "conditional_topology": (
                    "SELECT metric_order, metric, rate, numerator, denominator "
                    "FROM conditional_topology ORDER BY metric_order;",
                    ["Rates are conditional on exact shared/evaluable universes and do not equal global coverage."],
                ),
                "labeled_class_summary": (
                    "SELECT layer, pair_class, similarity FROM labeled_class_summary "
                    "ORDER BY layer, pair_class;",
                    ["Same-cancer and cross-cancer values use 6-biological-ID macro averages."],
                ),
                "profile_class_summary": (
                    "SELECT layer, pair_class, similarity FROM profile_class_summary "
                    "ORDER BY layer, pair_class;",
                    ["VAF and topology profiles compare distributions, not locus identity or cellular clones."],
                ),
                "permutation_tests": (
                    "SELECT metric_order, metric, within_mean, between_mean, within_minus_between, "
                    "permutation_n, one_sided_exact_p, interpretation "
                    "FROM permutation_tests ORDER BY metric_order;",
                    ["Exact one-sided label permutation over all 60 group assignments with sizes 3/2/1."],
                ),
                "robustness": (
                    "SELECT \"order\", \"check\", estimate, lower, upper, note "
                    "FROM robustness ORDER BY \"order\";",
                    ["Autosomes are resampled as 22 blocks with replacement for the composite profile."],
                ),
                "pairwise_metrics": (
                    "SELECT pair_order, pair, pair_class, candidate_jaccard, active_jaccard, "
                    "w_jaccard, primary_edge_jaccard, exact_component_jaccard, "
                    "vaf_js_similarity, coarse_topology_js_similarity, "
                    "composite_profile_similarity, exact_coordinate_1to1_n, "
                    "k_ge_2_shape_equal_rate, ancestor_equal_rate "
                    "FROM pairwise_metrics ORDER BY pair_order;",
                    ["Exactly 21 upper-triangle technical-dataset pairs; no duplicated symmetric cells."],
                ),
            }
        )
        source_by_dataset: dict[str, dict[str, Any]] = {}
        for dataset, (sql, definitions) in query_specs.items():
            rows = execute_query(connection, sql)
            require(rows, f"query returned no rows: {dataset}")
            datasets[dataset] = rows
            source = source_record(
                f"source_{dataset}",
                dataset.replace("_", " "),
                sql,
                [dataset],
                definitions,
                executed_at,
            )
            source_by_dataset[dataset] = source
            sources.append(source)
    finally:
        connection.close()

    charts = [
        bar_chart(
            "global_ladder_chart",
            "global_ladder",
            "HCC1395 technical pair 的分層相似度",
            "底層位點與 profile 高度相似；global edge/component Jaccard 較低，不能概括成每條邊一致。",
            source_by_dataset["global_ladder"],
            "layer",
            "similarity",
            "layer",
        ),
        heatmap_chart(
            "candidate_heatmap",
            "candidate_matrix",
            "Candidate-site Jaccard matrix",
            "只有 HCC1395 × DORADO 具有大規模 exact allele-key overlap；不同 biological IDs 接近零。",
            source_by_dataset["candidate_matrix"],
        ),
        heatmap_chart(
            "active_heatmap",
            "active_matrix",
            "Active-site Jaccard matrix",
            "進入 strict graph 的 active loci 仍保留強烈 biological-ID specificity。",
            source_by_dataset["active_matrix"],
        ),
        heatmap_chart(
            "w_heatmap",
            "w_matrix",
            "W-site Jaccard matrix",
            "區域可辨識性下降使 HCC pair Jaccard 降至 0.428，但仍遠高於任何不同 ID。",
            source_by_dataset["w_matrix"],
        ),
        heatmap_chart(
            "edge_heatmap",
            "edge_matrix",
            "Primary undirected-edge Jaccard matrix",
            "HCC pair 為 0.198；這是全域 set overlap，不等於 shared-opportunity 的 conditional containment。",
            source_by_dataset["edge_matrix"],
        ),
        heatmap_chart(
            "component_heatmap",
            "component_matrix",
            "Exact-component Jaccard matrix",
            "HCC pair 仍排名第一；不同 ID 沒有相同 exact components。",
            source_by_dataset["component_matrix"],
        ),
        heatmap_chart(
            "vaf_heatmap",
            "vaf_matrix",
            "Truth-confirmed marginal raw-VAF JS-similarity matrix",
            "HCC pair 排名第一，但同癌種 group-level permutation p=0.30；頻譜相似不能代替 clone identity。",
            source_by_dataset["vaf_matrix"],
        ),
        heatmap_chart(
            "coarse_topology_heatmap",
            "coarse_topology_matrix",
            "Coarse topology-composition JS-similarity matrix",
            "同癌種不同 ID 呈較高 profile similarity；這是 shape composition，不是相同 locus tree。",
            source_by_dataset["coarse_topology_matrix"],
        ),
        heatmap_chart(
            "composite_profile_heatmap",
            "composite_profile_matrix",
            "Exact-PS 五維 composite-profile similarity matrix",
            "HCC pair point rank 1/21；same-cancer profile 亦靠近，顯示粗 profile 不是專一 clone 指紋。",
            source_by_dataset["composite_profile_matrix"],
        ),
        bar_chart(
            "conditional_topology_chart",
            "conditional_topology",
            "HCC exact-coordinate／shared-opportunity conditional topology",
            "k≥2 unique shape 89.51%、ancestor relation 79.72%；這些只適用於可 1:1 對齊的局部區域。",
            source_by_dataset["conditional_topology"],
            "metric",
            "rate",
            "metric",
        ),
        bar_chart(
            "labeled_class_chart",
            "labeled_class_summary",
            "Exact-labeled 位點／邊相似度",
            "same-ID 遠高於 same-cancer 與 cross-cancer；癌種沒有產生相同變異或相同 edge set。",
            source_by_dataset["labeled_class_summary"],
            "layer",
            "similarity",
            "pair_class",
        ),
        bar_chart(
            "profile_class_chart",
            "profile_class_summary",
            "VAF 與拓撲 profile 相似度",
            "same-cancer 的 topology profile 較接近，但 VAF group-level evidence 較弱。",
            source_by_dataset["profile_class_summary"],
            "layer",
            "similarity",
            "pair_class",
        ),
    ]

    tables = [
        table_definition(
            "sample_table",
            "sample_metadata",
            "7 technical datasets 與 6 biological IDs",
            "Organ、histology 與取樣部位分開理解；DORADO 不算第二個乳癌 biological replicate。",
            source_by_dataset["sample_metadata"],
            [
                {"field": "dataset_order", "label": "#", "format": "number"},
                {"field": "dataset", "label": "Technical dataset", "type": "text"},
                {"field": "biological_id", "label": "Biological ID", "type": "text"},
                {"field": "organ", "label": "Organ", "type": "text"},
                {"field": "histology", "label": "Histology", "type": "text"},
                {"field": "technical_role", "label": "Role", "type": "text"},
            ],
            "dataset_order",
        ),
        table_definition(
            "permutation_table",
            "permutation_tests",
            "6-biological-ID exact label-permutation tests",
            "60 個 3/2/1 group assignments；p 值為探索性，4 個 within-group pairs 並不足以一般化。",
            source_by_dataset["permutation_tests"],
            [
                {"field": "metric_order", "label": "#", "format": "number"},
                {"field": "metric", "label": "Metric", "type": "text"},
                {"field": "within_mean", "label": "Within mean", "format": "number"},
                {"field": "between_mean", "label": "Between mean", "format": "number"},
                {"field": "within_minus_between", "label": "Δ within−between", "format": "number"},
                {"field": "permutation_n", "label": "Permutations", "format": "number"},
                {"field": "one_sided_exact_p", "label": "One-sided p", "format": "number"},
                {"field": "interpretation", "label": "Interpretation", "type": "text"},
            ],
            "metric_order",
        ),
        table_definition(
            "robustness_table",
            "robustness",
            "Composite profile 的 chromosome-block robustness",
            "HCC point rank 1，但 rank-1 bootstrap rate 66.7%，且與次佳差的 interval 跨零。",
            source_by_dataset["robustness"],
            [
                {"field": "order", "label": "#", "format": "number"},
                {"field": "check", "label": "Check", "type": "text"},
                {"field": "estimate", "label": "Estimate", "format": "number"},
                {"field": "lower", "label": "2.5%", "format": "number"},
                {"field": "upper", "label": "97.5%", "format": "number"},
                {"field": "note", "label": "Note", "type": "text"},
            ],
            "order",
        ),
        table_definition(
            "pairwise_table",
            "pairwise_metrics",
            "21 組 technical-dataset pair 的逐步對帳",
            "Exact topology rates 對不同 biological IDs 為 N/A，因 exact-coordinate 1:1 evaluable n=0；不可把 N/A 當 0 agreement。",
            source_by_dataset["pairwise_metrics"],
            [
                {"field": "pair_order", "label": "#", "format": "number"},
                {"field": "pair", "label": "Pair", "type": "text"},
                {"field": "pair_class", "label": "Class", "type": "text"},
                {"field": "candidate_jaccard", "label": "Candidate J", "format": "number"},
                {"field": "active_jaccard", "label": "Active J", "format": "number"},
                {"field": "w_jaccard", "label": "W J", "format": "number"},
                {"field": "primary_edge_jaccard", "label": "Edge J", "format": "number"},
                {"field": "exact_component_jaccard", "label": "Component J", "format": "number"},
                {"field": "vaf_js_similarity", "label": "VAF JS sim", "format": "number"},
                {"field": "coarse_topology_js_similarity", "label": "Coarse topo sim", "format": "number"},
                {"field": "composite_profile_similarity", "label": "Composite sim", "format": "number"},
                {"field": "exact_coordinate_1to1_n", "label": "Exact 1:1 n", "format": "number"},
                {"field": "k_ge_2_shape_equal_rate", "label": "k≥2 shape", "format": "number"},
                {"field": "ancestor_equal_rate", "label": "Ancestor", "format": "number"},
            ],
            "pair_order",
        ),
    ]

    blocks: list[dict[str, Any]] = [
        {
            "id": "title",
            "type": "markdown",
            "body": (
                f"# {TITLE}\n\n"
                "Task Type B comprehensive validation · GRCh38 chr1–22 · "
                "7 technical datasets · 6 biological IDs · 21 technical pairs · "
                "15 biological-ID-deduplicated pairs"
            ),
        },
        {
            "id": "technical_summary",
            "type": "markdown",
            "body": (
                "## 技術摘要\n\n"
                "**結論：方法可清楚重現 HCC1395 biological ID；但「相似」必須分層。** "
                "HCC1395 × HCC1395_DORADO 的 candidate/active/W/edge Jaccard 分別為 "
                "92.76%／78.59%／42.75%／19.80%，全部排名 1/21。"
                "在 2,825 個 fail-closed 1:1 exact-coordinate units 中，k≥2 unique shape "
                "相等 89.51%，ancestor relation 相等 79.72%；shared-candidate universe 中 "
                "DORADO→HCC edge containment 為 98.28%。因此可說底層位點高度一致、"
                "可辨識局部拓撲高度重現，但不能說 global 每條 edge 或唯一融合樹已確認。\n\n"
                "不同 biological ID 但同癌種的 **coarse topology profile** 較相似"
                "（exact 60-label permutation p=1/60），但 candidate/edge p=0.10、"
                "VAF p=0.30，且只有 4 個去重 same-cancer pairs。這是探索性"
                " cancer-type-associated profile signal，不是相同 clone tree。"
            ),
        },
        {
            "id": "scope_heading",
            "type": "markdown",
            "body": (
                "## Cohort、癌種與獨立性\n\n"
                "HCC1395_DORADO 與 HCC1395 是同一 biological ID 的不同技術來源；"
                "不能當兩個獨立乳癌樣本。去重後的同癌種 pairs 為乳癌 3 組、肺腺癌 1 組；"
                "COLO829 沒有第二個 melanoma 對照。"
            ),
        },
        {"id": "sample_table_block", "type": "table", "tableId": "sample_table", "layout": "full"},
        {
            "id": "ladder_heading",
            "type": "markdown",
            "body": (
                "## HCC1395 same-ID 證據階梯\n\n"
                "Set Jaccard 回答「整體集合是否相同」；profile similarity 回答"
                "「分布輪廓是否相似」。兩者不能互換。global edge Jaccard 低，"
                "但 shared-opportunity containment 高，表示解析量與可觀測母群不同。"
            ),
        },
        {"id": "global_ladder_block", "type": "chart", "chartId": "global_ladder_chart"},
        {
            "id": "site_heading",
            "type": "markdown",
            "body": (
                "## 位點底層矩陣\n\n"
                "Candidate 與 active site 使用 exact allele-key set Jaccard。"
                "HCC pair 為唯一明顯高相似格；同癌種不同 ID 不共享相同突變背景。"
            ),
        },
        {"id": "candidate_block", "type": "chart", "chartId": "candidate_heatmap"},
        {"id": "active_block", "type": "chart", "chartId": "active_heatmap"},
        {
            "id": "region_heading",
            "type": "markdown",
            "body": (
                "## sSNV 區域、edge 與 component 矩陣\n\n"
                "W、primary edge 與 exact component 逐步加入 phase／partial-read／"
                "edge-identifiability 限制，因此 HCC pair 的全域 Jaccard 下降；"
                "它仍在每層排名第一。Primary edge 是 phase-invariant 無向 read linkage，"
                "不是演化 parent→child edge。"
            ),
        },
        {"id": "w_block", "type": "chart", "chartId": "w_heatmap"},
        {"id": "edge_block", "type": "chart", "chartId": "edge_heatmap"},
        {"id": "component_block", "type": "chart", "chartId": "component_heatmap"},
        {
            "id": "vaf_heading",
            "type": "markdown",
            "body": (
                "## VAF 分布矩陣\n\n"
                "此矩陣比較 truth-confirmed marginal raw-VAF 50-bin distributions，"
                "位點不必相同。HCC pair 是最近鄰，但其他 cell line 也可有相似頻譜；"
                "VAF profile 不能單獨識別 clone。"
            ),
        },
        {"id": "vaf_block", "type": "chart", "chartId": "vaf_heatmap"},
        {
            "id": "topology_profile_heading",
            "type": "markdown",
            "body": (
                "## 拓撲 profile 矩陣\n\n"
                "Coarse matrix 比較 Single/Sister/Direct/Sister+direct/Unresolved 的 composition；"
                "composite matrix 再加入 status、resolution、active-k 與 HP-symmetric profile。"
                "HCC pair point estimate 均排名第一，但 same-cancer pairs 也較靠近，"
                "所以 profile similarity 不是 biological-ID 專屬指紋。"
            ),
        },
        {"id": "coarse_topology_block", "type": "chart", "chartId": "coarse_topology_heatmap"},
        {"id": "composite_profile_block", "type": "chart", "chartId": "composite_profile_heatmap"},
        {
            "id": "exact_topology_heading",
            "type": "markdown",
            "body": (
                "## Exact-coordinate 局部拓撲\n\n"
                "只有 HCC same-ID pair 有可觀的 exact-coordinate 1:1 母群；"
                "其餘 20 pairs 的 evaluable n 都是 0，因此不能以空白矩陣判為不一致。"
                "下圖只呈現 HCC pair 的條件式局部拓撲，且 k=1 不納入 k≥2 shape／ancestor 主張。"
            ),
        },
        {"id": "conditional_topology_block", "type": "chart", "chartId": "conditional_topology_chart"},
        {
            "id": "same_cancer_heading",
            "type": "markdown",
            "body": (
                "## 同 ID、同癌種與跨癌種\n\n"
                "Exact-labeled 層顯示同 ID 的壓倒性 specificity；distribution-profile 層則出現"
                " same-ID > same-cancer > cross-cancer 的描述性梯度。這兩個結果並不衝突："
                "同癌種可能共享複雜度／形狀分布，卻不共享相同突變、edge 或 clone ancestry。"
            ),
        },
        {"id": "labeled_class_block", "type": "chart", "chartId": "labeled_class_chart"},
        {"id": "profile_class_block", "type": "chart", "chartId": "profile_class_chart"},
        {
            "id": "permutation_heading",
            "type": "markdown",
            "body": (
                "## 同癌種 exact label-permutation\n\n"
                "HCC1395 的兩技術來源先等權 macro-average，再以 6 biological IDs 枚舉"
                "全部 60 種 3/2/1 分組。最小可能 p 值是 1/60；本檢定為 post-hoc、"
                "sample-size 極小，不能當作癌種效應的 confirmatory proof。"
            ),
        },
        {"id": "permutation_table_block", "type": "table", "tableId": "permutation_table", "layout": "full"},
        {
            "id": "robustness_heading",
            "type": "markdown",
            "body": (
                "## Robustness 與排名不確定度\n\n"
                "Composite profile 的 point rank 是 1/21，但 22-autosome block bootstrap 中"
                " rank-1 rate 只有 66.7%、top-3 rate 98.6%，且與次佳 pair 的差值 interval 跨零。"
                "因此可說 HCC pair 屬最相似群，但不能說粗 profile 永遠唯一第一。"
            ),
        },
        {"id": "robustness_table_block", "type": "table", "tableId": "robustness_table", "layout": "full"},
        {
            "id": "pairwise_heading",
            "type": "markdown",
            "body": (
                "## 21 組配對逐步明細\n\n"
                "表格依 candidate Jaccard 由高到低排序；N/A 表示沒有可評估共同母群，"
                "不是 0 agreement。"
            ),
        },
        {"id": "pairwise_table_block", "type": "table", "tableId": "pairwise_table", "layout": "full"},
        {
            "id": "methodology",
            "type": "markdown",
            "body": (
                "## 方法與 metric 定義\n\n"
                "Candidate/active/W/edge/component 採 symmetric Jaccard。VAF 與 topology "
                "JS similarity 採 `1 − sqrt(JSD_nats / ln 2)`。Coarse topology 只比較"
                "五類 shape composition；composite profile 對 status、resolution、morphology、"
                "active-k、HP-symmetric 五維等權。Exact-coordinate topology 先排除任一側"
                "重複 coordinate signature，再於 1:1 units 比較完整 shape set、k≥2 unique shape "
                "與 recurrence-free unique-tree ancestor relation。"
            ),
        },
        {
            "id": "limitations",
            "type": "markdown",
            "body": (
                "## 限制與主張上限\n\n"
                "只有一組 same-ID positive pair；same-cancer 去重後只有 4 pairs，且癌種與"
                "平台、coverage、phase fragmentation、read length、passage／aliquot 可能混雜。"
                "本輪沒有單細胞或多時間點 clone truth，也沒有跨區域 private-node contraction、"
                "partial-order consensus 與 conflict-aware fusion。因此可驗證的是位點、局部"
                "read-linkage、條件式數學拓撲與 dataset-level profile；不能確認 cellular clone K、"
                "全域 biological ancestry 或唯一 fusion tree。"
            ),
        },
        {
            "id": "next_steps",
            "type": "markdown",
            "body": (
                "## 下一步\n\n"
                "1. 對 HCC pair 做 matched-depth、PS/opportunity 與 read-length ablation，量化"
                " phase／region／partial-read 各自貢獻。\n"
                "2. 以 shared loci 建三態 ancestor matrix（支持／反對／不可辨識），再做"
                " bootstrap consensus，而非硬選單棵樹。\n"
                "3. 新增至少第二組 same-ID technical replicate，並擴大每癌種 biological IDs。\n"
                "4. 最後才以 CN/CCF、SEQC2-compatible context 或單細胞資料驗證 global clone tree。"
            ),
        },
        {
            "id": "further_questions",
            "type": "markdown",
            "body": (
                "## 後續可回答的問題\n\n"
                "粗 topology cancer-group 訊號是否在 matched coverage／PS fragmentation 後仍存在？"
                "同一 biological ID 的第二組 replicate 是否也在 exact-labeled 層排名第一？"
                "局部 79.72% ancestor agreement 中，衝突集中於 k=3+、CN/LOH 或短 phase blocks 嗎？"
            ),
        },
    ]

    manifest_sources: list[dict[str, Any]] = [
        {
            "id": "matrix_receipt",
            "label": "Cross-sample matrix validation receipt v2",
            "path": "research/20260725_crosssample_topology_matrix_report/results_v2/validation_receipt.json",
        },
        {
            "id": "cohort_similarity_source",
            "label": "Exact-PS cohort profile and HCC chromosome-block bootstrap",
            "path": "output/synthesis/observation_workspaces/20260725_exact_ps_cohort_similarity/all7_v1/cohort_similarity.v1.json",
        },
        {
            "id": "seqc2_reference",
            "label": "SEQC2 HCC1395 somatic mutation benchmark reference",
            "href": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8532138/",
        },
        {
            "id": "atcc_hcc1395",
            "label": "ATCC HCC1395",
            "href": "https://www.atcc.org/products/crl-2324",
        },
        {
            "id": "atcc_colo829",
            "label": "ATCC COLO829",
            "href": "https://www.atcc.org/products/crl-1974",
        },
        {
            "id": "atcc_h1437",
            "label": "ATCC H1437",
            "href": "https://www.atcc.org/products/crl-5872",
        },
        {
            "id": "atcc_h2009",
            "label": "ATCC H2009",
            "href": "https://www.atcc.org/products/crl-5911",
        },
        {
            "id": "atcc_hcc1937",
            "label": "ATCC HCC1937",
            "href": "https://www.atcc.org/products/crl-2336",
        },
        {
            "id": "atcc_hcc1954",
            "label": "ATCC HCC1954",
            "href": "https://www.atcc.org/products/crl-2338",
        },
        *sources,
    ]
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": TITLE,
            "description": (
                "Seven-dataset chr1–22 sSNV, regional linkage, VAF, and topology "
                "comparison matrices with biological-ID-deduplicated cancer-group tests."
            ),
            "generatedAt": executed_at,
            "blocks": blocks,
            "charts": charts,
            "tables": tables,
            "sources": manifest_sources,
        },
        "snapshot": {
            "version": 1,
            "status": "ready",
            "generatedAt": executed_at,
            "datasets": datasets,
        },
        "sources": manifest_sources,
        "package_info": {
            "root": ".",
            "manifestPath": "artifact.json",
            "snapshotPath": "artifact.json",
        },
    }
    args.artifact_json.parent.mkdir(parents=True, exist_ok=True)
    with args.artifact_json.open("x", encoding="utf-8") as handle:
        json.dump(artifact, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    require(len(datasets) <= 50, "snapshot dataset count exceeds 50")
    require(all(len(rows) <= 2000 for rows in datasets.values()), "snapshot row bound exceeded")
    print(
        json.dumps(
            {
                "status": "PASS",
                "artifact": str(args.artifact_json.resolve()),
                "sqlite": str(args.sqlite_db.resolve()),
                "datasets": len(datasets),
                "snapshot_rows": sum(len(rows) for rows in datasets.values()),
                "charts": len(charts),
                "tables": len(tables),
                "source_queries": len(sources),
                "profile_permutation_p": profile_test["one_sided_exact_p"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BuildError as exc:
        print(f"ERROR: {exc}")
        raise SystemExit(2)
