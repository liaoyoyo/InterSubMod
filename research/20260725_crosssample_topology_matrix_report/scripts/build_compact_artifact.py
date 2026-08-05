#!/usr/bin/env python3
"""Compact a validated report artifact into a browser-friendly visual set."""

from __future__ import annotations

import argparse
import json
import sqlite3
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


OLD_DB_PATH = "research/20260725_crosssample_topology_matrix_report/data/report.sqlite"
NEW_DB_PATH = (
    "research/20260725_crosssample_topology_matrix_report/data/report_compact.sqlite"
)
KEEP_CHARTS = {
    "candidate_heatmap",
    "composite_profile_heatmap",
    "conditional_topology_chart",
    "profile_class_chart",
}
KEEP_TABLES = {
    "sample_table",
    "permutation_table",
    "robustness_table",
    "pairwise_table",
}


class CompactError(RuntimeError):
    """Raised when the compact artifact contract fails."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CompactError(message)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-artifact", required=True, type=Path)
    parser.add_argument("--input-sqlite", required=True, type=Path)
    parser.add_argument("--output-artifact", required=True, type=Path)
    parser.add_argument("--output-sqlite", required=True, type=Path)
    return parser.parse_args()


def replace_source_path(value: Any) -> None:
    if isinstance(value, dict):
        if value.get("path") == OLD_DB_PATH:
            value["path"] = NEW_DB_PATH
        for child in value.values():
            replace_source_path(child)
    elif isinstance(value, list):
        for child in value:
            replace_source_path(child)


def main() -> int:
    args = parse_args()
    require(args.input_artifact.is_file(), f"missing artifact: {args.input_artifact}")
    require(args.input_sqlite.is_file(), f"missing SQLite DB: {args.input_sqlite}")
    require(not args.output_artifact.exists(), f"output exists: {args.output_artifact}")
    require(not args.output_sqlite.exists(), f"output exists: {args.output_sqlite}")

    artifact = json.loads(args.input_artifact.read_text(encoding="utf-8"))
    pair_rows = artifact["snapshot"]["datasets"]["pairwise_metrics"]
    require(len(pair_rows) == 21, "pairwise_metrics row count != 21")
    pair_fields = [f"p{index:02d}" for index in range(1, 22)]
    pair_labels = {
        field: pair_rows[index]["pair"] for index, field in enumerate(pair_fields)
    }
    step_specs = (
        ("Candidate sites", "candidate_jaccard"),
        ("Active sites", "active_jaccard"),
        ("W sites", "w_jaccard"),
        ("Primary edges", "primary_edge_jaccard"),
        ("Exact components", "exact_component_jaccard"),
        ("VAF profile", "vaf_js_similarity"),
        ("Coarse topology", "coarse_topology_js_similarity"),
        ("Composite topology", "composite_profile_similarity"),
        ("k≥2 unique shape", "k_ge_2_shape_equal_rate"),
        ("Ancestor relation", "ancestor_equal_rate"),
    )
    step_rows: list[dict[str, Any]] = []
    for order, (layer, value_field) in enumerate(step_specs, 1):
        row: dict[str, Any] = {"row_order": order, "layer": layer}
        for pair_field, pair_row in zip(pair_fields, pair_rows):
            row[pair_field] = pair_row[value_field]
        step_rows.append(row)

    args.output_sqlite.parent.mkdir(parents=True, exist_ok=True)
    source_connection = sqlite3.connect(args.input_sqlite)
    target_connection = sqlite3.connect(args.output_sqlite)
    try:
        source_connection.backup(target_connection)
        target_connection.execute(
            "CREATE TABLE step_pair_matrix ("
            + '"row_order" INTEGER, "layer" TEXT, '
            + ", ".join(f'"{field}" REAL' for field in pair_fields)
            + ")"
        )
        target_connection.executemany(
            "INSERT INTO step_pair_matrix VALUES ("
            + ", ".join("?" for _ in range(2 + len(pair_fields)))
            + ")",
            [
                tuple([row["row_order"], row["layer"]] + [row[field] for field in pair_fields])
                for row in step_rows
            ],
        )
        target_connection.commit()
        cursor = target_connection.execute(
            "SELECT "
            + ", ".join(['"layer"'] + [f'"{field}"' for field in pair_fields])
            + " FROM step_pair_matrix ORDER BY row_order"
        )
        columns = [description[0] for description in cursor.description]
        reviewed_rows = [dict(zip(columns, row)) for row in cursor.fetchall()]
    finally:
        source_connection.close()
        target_connection.close()
    require(len(reviewed_rows) == 10, "step-pair query row count != 10")

    compact = deepcopy(artifact)
    generated_at = datetime.now(timezone.utc).isoformat()
    compact["manifest"]["generatedAt"] = generated_at
    compact["snapshot"]["generatedAt"] = generated_at
    compact["snapshot"]["datasets"]["step_pair_matrix"] = reviewed_rows
    replace_source_path(compact)

    query_sql = (
        "SELECT "
        + ", ".join(['"layer"'] + [f'"{field}"' for field in pair_fields])
        + " FROM step_pair_matrix ORDER BY row_order;"
    )
    step_source = {
        "id": "source_step_pair_matrix",
        "label": "21 technical pairs × 10 evidence steps",
        "path": NEW_DB_PATH,
        "query": {
            "engine": "sqlite",
            "language": "sql",
            "sql": query_sql,
            "description": "Executed and reviewed compact pair-by-step similarity matrix.",
            "tables_used": ["step_pair_matrix"],
            "filters": [
                "GRCh38 chr1–22",
                "21 upper-triangle technical-dataset pairs",
            ],
            "metric_definitions": [
                "Rows retain metric-specific denominators; values share a 0–1 display scale.",
                "N/A means no exact-coordinate evaluable universe, not zero agreement.",
            ],
            "executed_at": generated_at,
        },
    }
    step_chart = {
        "id": "step_pair_heatmap",
        "dataset": "step_pair_matrix",
        "type": "heatmap",
        "title": "21 technical pairs × 10 evidence steps",
        "subtitle": (
            "HCC1395 same-ID pair is the only pair with high exact-labeled overlap; "
            "same-cancer structure appears mainly in distribution-level topology profiles."
        ),
        "intent": "relationship",
        "question": "每一組 technical pair 在各證據步驟的相似度如何？",
        "rationale": (
            "One compact matrix preserves every pair-by-step value while avoiding "
            "repeated 7×7 chart rendering."
        ),
        "comparisonContext": {
            "grain": "evidence step × technical-dataset pair",
            "denominator": "metric-specific",
            "unit": "0–1 similarity / agreement",
        },
        "encodings": {
            "x": {"field": "layer", "type": "nominal", "label": "Evidence step"},
            "y": {"fields": pair_fields, "labels": pair_labels},
        },
        "valueFormat": "percent",
        "palette": {"kind": "sequential"},
        "sourceId": step_source["id"],
        "source": deepcopy(step_source),
    }

    chart_by_id = {
        chart["id"]: chart
        for chart in compact["manifest"]["charts"]
        if chart["id"] in KEEP_CHARTS
    }
    require(set(chart_by_id) == KEEP_CHARTS, "compact chart selection incomplete")
    compact["manifest"]["charts"] = [
        step_chart,
        chart_by_id["candidate_heatmap"],
        chart_by_id["composite_profile_heatmap"],
        chart_by_id["conditional_topology_chart"],
        chart_by_id["profile_class_chart"],
    ]
    compact["manifest"]["tables"] = [
        table
        for table in compact["manifest"]["tables"]
        if table["id"] in KEEP_TABLES
    ]
    require(
        {table["id"] for table in compact["manifest"]["tables"]} == KEEP_TABLES,
        "compact table selection incomplete",
    )

    original_blocks = {block["id"]: block for block in compact["manifest"]["blocks"]}
    original_blocks["ladder_heading"]["body"] = (
        "## 全部配對 × 全部步驟\n\n"
        "總矩陣同時列出 21 組 technical pairs 與 candidate→active→W→edge→component→"
        "VAF→topology 的完整階梯。各列使用自己的 denominator；最後兩列只有 HCC pair "
        "有 exact-coordinate evaluable 母群。"
    )
    compact["manifest"]["blocks"] = [
        original_blocks["title"],
        original_blocks["technical_summary"],
        original_blocks["scope_heading"],
        original_blocks["sample_table_block"],
        original_blocks["ladder_heading"],
        {"id": "step_pair_block", "type": "chart", "chartId": "step_pair_heatmap"},
        original_blocks["site_heading"],
        {"id": "candidate_block", "type": "chart", "chartId": "candidate_heatmap"},
        original_blocks["topology_profile_heading"],
        {
            "id": "composite_profile_block",
            "type": "chart",
            "chartId": "composite_profile_heatmap",
        },
        original_blocks["exact_topology_heading"],
        {
            "id": "conditional_topology_block",
            "type": "chart",
            "chartId": "conditional_topology_chart",
        },
        original_blocks["same_cancer_heading"],
        {"id": "profile_class_block", "type": "chart", "chartId": "profile_class_chart"},
        original_blocks["permutation_heading"],
        original_blocks["permutation_table_block"],
        original_blocks["robustness_heading"],
        original_blocks["robustness_table_block"],
        original_blocks["pairwise_heading"],
        original_blocks["pairwise_table_block"],
        original_blocks["methodology"],
        original_blocks["limitations"],
        original_blocks["next_steps"],
        original_blocks["further_questions"],
    ]
    compact["manifest"]["sources"].append(step_source)
    compact["sources"].append(step_source)
    compact["package_info"] = {
        "root": ".",
        "manifestPath": args.output_artifact.name,
        "snapshotPath": args.output_artifact.name,
    }

    with args.output_artifact.open("x", encoding="utf-8") as handle:
        json.dump(compact, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print(
        json.dumps(
            {
                "status": "PASS",
                "artifact": str(args.output_artifact.resolve()),
                "sqlite": str(args.output_sqlite.resolve()),
                "charts": len(compact["manifest"]["charts"]),
                "tables": len(compact["manifest"]["tables"]),
                "blocks": len(compact["manifest"]["blocks"]),
                "step_matrix_rows": len(reviewed_rows),
                "step_matrix_pairs": len(pair_fields),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CompactError as exc:
        print(f"ERROR: {exc}")
        raise SystemExit(2)
