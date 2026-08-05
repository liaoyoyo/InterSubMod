#!/usr/bin/env python3
"""Create the canonical Data Analytics artifact for the strict-region HTML report."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sqlite3
from typing import Any, Mapping, Sequence
from uuid import uuid4


TITLE = "Exact PS × HP 嚴格 read-linkage 區域：7 資料集全量報告"
RAW_SOURCE_ID = "all7_strict_region_data"
CANONICAL_DATASETS = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
AUTOSOMES = tuple(f"chr{number}" for number in range(1, 23))
PRIMARY_THRESHOLD = 3
ARTIFACT_RECEIPT_NAME = "artifact_receipt.json"
ARTIFACT_RECEIPT_SIDECAR_NAME = "artifact_receipt.json.sha256"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: JSON root must be an object")
    return value


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def verify_file_identity(spec: Mapping[str, Any], path: Path) -> None:
    observed = file_identity(path)
    if any(spec.get(key) != observed[key] for key in observed):
        raise ValueError(f"{path}: identity mismatch against report-data receipt")


def fmt_int(value: int) -> str:
    return f"{value:,}"


def pct_fraction(value: float | int | None) -> float | None:
    return None if value is None else float(value) / 100.0


def field(field: str, label: str, *, kind: str = "number") -> dict[str, str]:
    result = {"field": field, "label": label}
    if kind == "text":
        result["type"] = "text"
    return result


def source_id(dataset: str) -> str:
    return f"source_{dataset}"


def chart(
    *,
    chart_id: str,
    title: str,
    subtitle: str,
    dataset: str,
    chart_type: str,
    x_field: str,
    x_type: str,
    x_label: str,
    y_field: str,
    y_label: str,
    color_field: str | None = None,
    color_label: str | None = None,
    value_format: str | None = None,
    grouping: str | None = None,
    orientation: str = "vertical",
    intent: str,
    question: str,
    rationale: str,
    tooltip: Sequence[Mapping[str, str]] = (),
) -> dict[str, Any]:
    encodings: dict[str, Any] = {
        "x": {"field": x_field, "type": x_type, "label": x_label},
        "y": {"field": y_field, "type": "quantitative", "label": y_label},
    }
    if color_field is not None:
        encodings["color"] = {
            "field": color_field,
            "type": "nominal",
            "label": color_label or color_field,
        }
    if tooltip:
        encodings["tooltip"] = list(tooltip)
    result: dict[str, Any] = {
        "id": chart_id,
        "title": title,
        "subtitle": subtitle,
        "type": chart_type,
        "dataset": dataset,
        "sourceId": source_id(dataset),
        "encodings": encodings,
        "intent": intent,
        "question": question,
        "rationale": rationale,
        "options": {"orientation": orientation},
    }
    if grouping is not None:
        result["options"]["grouping"] = grouping
    if value_format is not None:
        result["valueFormat"] = value_format
    return result


def method_svg(total_w: int, singleton: int) -> str:
    return f"""
<section aria-label="Strict read-linkage region method" style="padding:8px 0 2px">
  <svg viewBox="0 0 1240 350" role="img" aria-labelledby="flow-title flow-desc"
       style="width:100%;height:auto;display:block;font-family:Inter,'Noto Sans TC',sans-serif">
    <title id="flow-title">Exact PS × HP strict read-linkage region 流程</title>
    <desc id="flow-desc">候選 sSNV 經 exact PS 與 HP 分容器，以至少三個 distinct molecules
      的 fixed R/A endpoint observation 建邊，再取最大連通分量；單點排除，兩點以上定義為 W。</desc>
    <defs>
      <marker id="arrow" markerWidth="10" markerHeight="10" refX="8" refY="3"
              orient="auto" markerUnits="strokeWidth">
        <path d="M0,0 L0,6 L9,3 z" fill="#46556b"/>
      </marker>
    </defs>
    <rect x="20" y="55" width="190" height="96" rx="16" fill="#eef4f8" stroke="#35556f" stroke-width="2"/>
    <text x="115" y="88" text-anchor="middle" font-size="20" font-weight="700" fill="#17324a">候選 sSNV</text>
    <text x="115" y="117" text-anchor="middle" font-size="15" fill="#405467">chr1–22 · PASS</text>
    <text x="115" y="140" text-anchor="middle" font-size="14" fill="#657687">物理位點 S</text>

    <line x1="210" y1="103" x2="255" y2="103" stroke="#46556b" stroke-width="3" marker-end="url(#arrow)"/>
    <rect x="265" y="55" width="205" height="96" rx="16" fill="#eef4f8" stroke="#35556f" stroke-width="2"/>
    <text x="368" y="88" text-anchor="middle" font-size="20" font-weight="700" fill="#17324a">合法容器</text>
    <text x="368" y="117" text-anchor="middle" font-size="15" fill="#405467">dataset × chromosome</text>
    <text x="368" y="140" text-anchor="middle" font-size="15" fill="#405467">× exact PS × HP1 / HP2</text>

    <line x1="470" y1="103" x2="515" y2="103" stroke="#46556b" stroke-width="3" marker-end="url(#arrow)"/>
    <rect x="525" y="55" width="235" height="96" rx="16" fill="#fff6df" stroke="#a66a09" stroke-width="2"/>
    <text x="643" y="86" text-anchor="middle" font-size="20" font-weight="700" fill="#573b0b">直接 read 邊</text>
    <text x="643" y="115" text-anchor="middle" font-size="15" fill="#5e513b">兩端皆 fixed R/A</text>
    <text x="643" y="140" text-anchor="middle" font-size="15" fill="#5e513b">≥3 distinct molecules</text>

    <line x1="760" y1="103" x2="805" y2="103" stroke="#46556b" stroke-width="3" marker-end="url(#arrow)"/>
    <rect x="815" y="55" width="205" height="96" rx="16" fill="#e8f5ef" stroke="#28715a" stroke-width="2"/>
    <text x="918" y="86" text-anchor="middle" font-size="20" font-weight="700" fill="#164f3d">連通分量</text>
    <text x="918" y="115" text-anchor="middle" font-size="15" fill="#375f53">maximal component</text>
    <text x="918" y="140" text-anchor="middle" font-size="15" fill="#375f53">距離只記錄、不切割</text>

    <line x1="1020" y1="103" x2="1060" y2="103" stroke="#46556b" stroke-width="3"/>
    <line x1="1060" y1="103" x2="1060" y2="52" stroke="#46556b" stroke-width="3"/>
    <line x1="1060" y1="103" x2="1060" y2="217" stroke="#46556b" stroke-width="3"/>
    <line x1="1060" y1="52" x2="1092" y2="52" stroke="#46556b" stroke-width="3" marker-end="url(#arrow)"/>
    <line x1="1060" y1="217" x2="1092" y2="217" stroke="#46556b" stroke-width="3" marker-end="url(#arrow)"/>

    <rect x="1100" y="15" width="120" height="75" rx="14" fill="#f1f3f5" stroke="#78838d" stroke-width="2"/>
    <text x="1160" y="45" text-anchor="middle" font-size="18" font-weight="700" fill="#3f4a54">k=1</text>
    <text x="1160" y="70" text-anchor="middle" font-size="14" fill="#596771">abstain</text>

    <rect x="1100" y="178" width="120" height="78" rx="14" fill="#dff2e8" stroke="#28715a" stroke-width="2"/>
    <text x="1160" y="208" text-anchor="middle" font-size="18" font-weight="700" fill="#164f3d">k≥2</text>
    <text x="1160" y="234" text-anchor="middle" font-size="14" fill="#375f53">W · 建樹候選</text>

    <rect x="20" y="284" width="1200" height="45" rx="12" fill="#f7f8fa" stroke="#d9dee3"/>
    <circle cx="50" cy="307" r="6" fill="#28715a"/>
    <text x="67" y="313" font-size="15" fill="#3e4d5a">全資料集 W={fmt_int(total_w)}</text>
    <circle cx="255" cy="307" r="6" fill="#78838d"/>
    <text x="272" y="313" font-size="15" fill="#3e4d5a">k=1 audit components={fmt_int(singleton)}</text>
    <circle cx="595" cy="307" r="6" fill="#a66a09"/>
    <text x="612" y="313" font-size="15" fill="#3e4d5a">edge 的座標距離不參與納入；另列 50 kb QC</text>
  </svg>
</section>
""".strip()


def quote_identifier(value: str) -> str:
    if not value or "\x00" in value:
        raise ValueError(f"invalid SQLite identifier: {value!r}")
    return '"' + value.replace('"', '""') + '"'


def sqlite_column_type(values: Sequence[Any]) -> str:
    nonmissing = [value for value in values if value is not None]
    if not nonmissing:
        return "TEXT"
    if all(isinstance(value, (bool, int)) and not isinstance(value, float) for value in nonmissing):
        return "INTEGER"
    if all(isinstance(value, (bool, int, float)) for value in nonmissing):
        return "REAL"
    return "TEXT"


def write_text_exclusive(path: Path, text: str) -> None:
    with path.open("x", encoding="utf-8") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())


def write_json_exclusive(path: Path, value: Mapping[str, Any]) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(
            value,
            handle,
            ensure_ascii=False,
            sort_keys=True,
            indent=2,
            allow_nan=False,
        )
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def artifact_release_checks(artifact: Mapping[str, Any]) -> dict[str, bool]:
    """Validate the structural release scope independently of the renderer."""
    manifest = artifact.get("manifest")
    snapshot = artifact.get("snapshot")
    if not isinstance(manifest, Mapping) or not isinstance(snapshot, Mapping):
        raise ValueError("artifact manifest/snapshot contract is absent")
    snapshot_datasets = snapshot.get("datasets")
    if not isinstance(snapshot_datasets, Mapping):
        raise ValueError("artifact snapshot datasets are absent")

    dataset_rows = snapshot_datasets.get("dataset_summary")
    chromosome_rows = snapshot_datasets.get("chromosome_detail")
    if not isinstance(dataset_rows, list) or not isinstance(chromosome_rows, list):
        raise ValueError("artifact release scope datasets are absent")
    dataset_names = [
        row.get("dataset") if isinstance(row, Mapping) else None
        for row in dataset_rows
    ]
    chromosome_pairs = [
        (
            row.get("dataset") if isinstance(row, Mapping) else None,
            row.get("chrom") if isinstance(row, Mapping) else None,
        )
        for row in chromosome_rows
    ]
    expected_pairs = {
        (dataset, chrom)
        for dataset in CANONICAL_DATASETS
        for chrom in AUTOSOMES
    }

    charts = manifest.get("charts")
    tables = manifest.get("tables")
    blocks = manifest.get("blocks")
    if not all(isinstance(value, list) for value in (charts, tables, blocks)):
        raise ValueError("artifact chart/table/block arrays are absent")
    chart_ids = [
        item.get("id") if isinstance(item, Mapping) else None for item in charts
    ]
    table_ids = [
        item.get("id") if isinstance(item, Mapping) else None for item in tables
    ]
    visible_chart_ids = [
        block.get("chartId")
        for block in blocks
        if isinstance(block, Mapping) and block.get("type") == "chart"
    ]
    visible_table_ids = [
        block.get("tableId")
        for block in blocks
        if isinstance(block, Mapping) and block.get("type") == "table"
    ]
    checks = {
        "surface_is_report": (
            artifact.get("surface") == "report"
            and manifest.get("surface") == "report"
        ),
        "snapshot_status_ready": snapshot.get("status") == "ready",
        "dataset_scope_exact_7": (
            len(dataset_names) == len(CANONICAL_DATASETS)
            and set(dataset_names) == set(CANONICAL_DATASETS)
        ),
        "chromosome_scope_exact_7x22": (
            len(chromosome_pairs) == len(expected_pairs)
            and set(chromosome_pairs) == expected_pairs
        ),
        "chart_ids_nonempty_unique": (
            bool(chart_ids)
            and all(isinstance(value, str) and value for value in chart_ids)
            and len(set(chart_ids)) == len(chart_ids)
        ),
        "table_ids_nonempty_unique": (
            bool(table_ids)
            and all(isinstance(value, str) and value for value in table_ids)
            and len(set(table_ids)) == len(table_ids)
        ),
        "visible_chart_ids_match_manifest": (
            sorted(visible_chart_ids) == sorted(chart_ids)
            and len(visible_chart_ids) == len(chart_ids)
        ),
        "visible_table_ids_match_manifest": (
            sorted(visible_table_ids) == sorted(table_ids)
            and len(visible_table_ids) == len(table_ids)
        ),
    }
    if not all(checks.values()):
        raise ValueError(f"artifact release checks failed: {checks}")
    return checks


def validate_report_data_package(data: Mapping[str, Any], data_path: Path) -> None:
    """Fail closed on scope, row uniqueness, checks, and the sibling receipt."""
    if (
        data.get("schema_name") != "intersubmod.strict_region_all7_report_data"
        or data.get("schema_version") != "1.1.0"
        or data.get("all_pass") is not True
    ):
        raise ValueError("report data schema/version/PASS contract mismatch")
    scope = data.get("scope")
    if not isinstance(scope, Mapping) or (
        scope.get("task_type") != "B_comprehensive_validation"
        or scope.get("datasets") != list(CANONICAL_DATASETS)
        or scope.get("chromosomes") != list(AUTOSOMES)
        or scope.get("primary_minimum_distinct_molecules") != PRIMARY_THRESHOLD
    ):
        raise ValueError("report data full-scope contract mismatch")
    checks = data.get("checks")
    if (
        not isinstance(checks, Mapping)
        or not checks
        or any(value is not True for value in checks.values())
    ):
        raise ValueError("report data contains a missing or failed check")
    if checks.get("all_154_extraction_receipt_identities_verified") is not True:
        raise ValueError("report data lacks the verified 7×22 extraction receipt grid")

    datasets = data.get("datasets")
    chromosomes = data.get("chromosomes")
    thresholds = data.get("threshold_sensitivity")
    k_bands = data.get("k_bands")
    span_bands = data.get("span_bands")
    exact_k = data.get("exact_k")
    topology_completion = data.get("topology_completion")
    topology_summary = data.get("topology_completion_summary")
    if not all(
        isinstance(rows, list)
        for rows in (
            datasets,
            chromosomes,
            thresholds,
            k_bands,
            span_bands,
            exact_k,
            topology_completion,
        )
    ):
        raise ValueError("one or more required report datasets are absent")
    dataset_names = [row.get("dataset") for row in datasets]
    if dataset_names != list(CANONICAL_DATASETS) or any(
        row.get("all_pass") is not True for row in datasets
    ):
        raise ValueError("dataset summary does not contain exactly one canonical PASS row")
    expected_pairs = {
        (dataset, chrom)
        for dataset in CANONICAL_DATASETS
        for chrom in AUTOSOMES
    }
    observed_pairs = {
        (row.get("dataset"), row.get("chrom"))
        for row in chromosomes
    }
    if (
        len(chromosomes) != len(expected_pairs)
        or observed_pairs != expected_pairs
        or any(row.get("all_pass") is not True for row in chromosomes)
    ):
        raise ValueError("chromosome detail is not the exact canonical 7×22 PASS grid")
    expected_threshold_pairs = {
        (dataset, threshold)
        for dataset in CANONICAL_DATASETS
        for threshold in (1, 2, 3, 5)
    }
    observed_threshold_pairs = {
        (row.get("dataset"), row.get("minimum_distinct_molecules"))
        for row in thresholds
    }
    if (
        len(thresholds) != len(expected_threshold_pairs)
        or observed_threshold_pairs != expected_threshold_pairs
    ):
        raise ValueError("threshold sensitivity is not the canonical 7×{1,2,3,5} grid")
    for row in thresholds:
        threshold = row.get("minimum_distinct_molecules")
        expected_basis = (
            "independent_row_recomputation_and_receipt_reconciliation"
            if threshold == 3
            else "chromosome_receipt_reaggregation_and_summary_reconciliation"
        )
        if row.get("validation_basis") != expected_basis:
            raise ValueError(
                "threshold sensitivity validation_basis does not match its evidence path"
            )
    for label, rows, bands in (
        ("k", k_bands, {"2", "3", "4", "5", "6–8", "9–12", ">12"}),
        (
            "span",
            span_bands,
            {"≤10 kb", "10–25 kb", "25–50 kb", "50–100 kb", ">100 kb"},
        ),
    ):
        band_field = f"{label}_band"
        expected = {
            (dataset, band)
            for dataset in CANONICAL_DATASETS
            for band in bands
        }
        observed = {(row.get("dataset"), row.get(band_field)) for row in rows}
        if len(rows) != len(expected) or observed != expected:
            raise ValueError(f"{label}-band data is not the complete canonical grid")
    exact_k_pairs = {
        (row.get("dataset"), row.get("k"))
        for row in exact_k
    }
    if (
        len(exact_k) != len(exact_k_pairs)
        or any(
            dataset not in CANONICAL_DATASETS
            or not isinstance(k, int)
            or k < 2
            for dataset, k in exact_k_pairs
        )
        or {
            row.get("dataset")
            for row in exact_k
        }
        != set(CANONICAL_DATASETS)
    ):
        raise ValueError("exact-k data is incomplete, duplicated, or invalid")

    if (
        [row.get("dataset") for row in topology_completion]
        != list(CANONICAL_DATASETS)
        or not isinstance(topology_summary, Mapping)
        or topology_summary.get("strict_linkage_complete_datasets") != 7
        or topology_summary.get(
            "strict_directed_topology_production_complete_datasets"
        )
        != 0
        or topology_summary.get("cellular_clone_count_validated_datasets") != 0
        or any(
            row.get("strict_linkage_status") != "COMPLETE_22_OF_22"
            or row.get("strict_directed_topology_status")
            not in {"NOT_RUN", "NOT_RUN_LATEST_V4_PARTITION_ONLY"}
            or row.get("cellular_clone_count_status") != "NOT_VALIDATED"
            or row.get("exact_cellular_parent_child_status") != "NOT_VALIDATED"
            or row.get("cross_hp_or_technical_fusion_status") != "NOT_VALIDATED"
            for row in topology_completion
        )
    ):
        raise ValueError("topology completion boundary is incomplete or overstated")

    inputs = data.get("inputs")
    if not isinstance(inputs, Mapping) or set(inputs) != set(CANONICAL_DATASETS):
        raise ValueError("report data input provenance is incomplete")
    extraction_grid_count = 0
    for dataset, provenance in inputs.items():
        if (
            not isinstance(provenance, Mapping)
            or not isinstance(provenance.get("summary"), Mapping)
            or not isinstance(provenance.get("chromosome_receipts"), Mapping)
            or set(provenance["chromosome_receipts"]) != set(AUTOSOMES)
            or not isinstance(provenance.get("extraction_receipts"), Mapping)
            or set(provenance["extraction_receipts"]) != set(AUTOSOMES)
        ):
            raise ValueError(f"{dataset}: input provenance is incomplete")
        for chrom, receipt_identity in provenance["extraction_receipts"].items():
            if (
                not isinstance(receipt_identity, Mapping)
                or not isinstance(receipt_identity.get("path"), str)
                or not Path(receipt_identity["path"]).is_absolute()
                or not isinstance(receipt_identity.get("size_bytes"), int)
                or isinstance(receipt_identity.get("size_bytes"), bool)
                or receipt_identity["size_bytes"] < 0
                or not isinstance(receipt_identity.get("sha256"), str)
                or len(receipt_identity["sha256"]) != 64
                or any(
                    char not in "0123456789abcdef"
                    for char in receipt_identity["sha256"]
                )
            ):
                raise ValueError(
                    f"{dataset}/{chrom}: extraction receipt identity is invalid"
                )
            extraction_grid_count += 1
    if extraction_grid_count != len(CANONICAL_DATASETS) * len(AUTOSOMES):
        raise ValueError("report data extraction receipt grid is not exactly 7×22")

    supporting_inputs = data.get("supporting_inputs")
    if not isinstance(supporting_inputs, Mapping):
        raise ValueError("report data supporting-input provenance is absent")
    for key in (
        "topology_completion_audit",
        "topology_completion_audit_receipt",
        "topology_completion_audit_receipt_sidecar",
    ):
        spec = supporting_inputs.get(key)
        if (
            not isinstance(spec, Mapping)
            or not isinstance(spec.get("path"), str)
            or not Path(spec["path"]).is_absolute()
        ):
            raise ValueError(f"report data supporting identity is invalid: {key}")
        verify_file_identity(spec, Path(spec["path"]))

    receipt_path = data_path.parent / "receipt.json"
    sidecar_path = data_path.parent / "receipt.json.sha256"
    receipt = read_json(receipt_path)
    expected_sidecar = f"{sha256_path(receipt_path)}  receipt.json"
    if sidecar_path.read_text(encoding="ascii").strip() != expected_sidecar:
        raise ValueError(f"{sidecar_path}: report-data receipt checksum mismatch")
    if (
        receipt.get("schema_name")
        != "intersubmod.strict_region_all7_report_data_receipt"
        or receipt.get("schema_version") != "1.1.0"
        or receipt.get("all_pass") is not True
        or receipt.get("checks") != checks
    ):
        raise ValueError(f"{receipt_path}: report-data receipt contract mismatch")
    output_spec = receipt.get("outputs", {}).get(data_path.name)
    if not isinstance(output_spec, Mapping):
        raise ValueError(f"{receipt_path}: data JSON identity is absent")
    verify_file_identity(output_spec, data_path)


def stage_and_query_snapshot(
    *,
    output_root: Path,
    raw_datasets: Mapping[str, Sequence[Mapping[str, Any]]],
    order_by: Mapping[str, Sequence[str]],
    generated_at: str,
) -> tuple[dict[str, list[dict[str, Any]]], list[dict[str, Any]]]:
    """Persist reviewed rows via an auditable, atomically published SQLite stage."""
    data_dir = output_root / "data"
    query_dir = output_root / "queries"
    data_dir.mkdir(parents=True, exist_ok=True)
    db_path = data_dir / "report.sqlite"
    if db_path.exists() or query_dir.exists():
        raise ValueError(
            f"SQLite/query staging outputs already exist: {db_path}, {query_dir}"
        )

    stage_id = uuid4().hex
    pending_db = output_root / f".report.sqlite.pending-{stage_id}"
    pending_queries = output_root / f".queries.pending-{stage_id}"
    pending_queries.mkdir(parents=False, exist_ok=False)
    queried: dict[str, list[dict[str, Any]]] = {}
    query_specs: list[tuple[str, str]] = []
    sources: list[dict[str, Any]] = []
    connection: sqlite3.Connection | None = None
    published_db = False
    published_queries = False
    query_executed_at = utc_now()
    try:
        connection = sqlite3.connect(pending_db)
        for dataset, input_rows in raw_datasets.items():
            rows = [dict(row) for row in input_rows]
            if not rows:
                raise ValueError(f"snapshot dataset must not be empty: {dataset}")
            columns = list(rows[0])
            if not columns or any(set(row) != set(columns) for row in rows):
                raise ValueError(f"snapshot dataset has inconsistent columns: {dataset}")
            column_types = {
                column: sqlite_column_type([row[column] for row in rows])
                for column in columns
            }
            column_ddl = ", ".join(
                f"{quote_identifier(column)} {column_types[column]}"
                for column in columns
            )
            connection.execute(
                f"CREATE TABLE {quote_identifier(dataset)} ({column_ddl})"
            )
            placeholders = ", ".join("?" for _ in columns)
            connection.executemany(
                (
                    f"INSERT INTO {quote_identifier(dataset)} "
                    f"({', '.join(quote_identifier(column) for column in columns)}) "
                    f"VALUES ({placeholders})"
                ),
                [
                    tuple(
                        int(value) if isinstance(value, bool) else value
                        for value in (row[column] for column in columns)
                    )
                    for row in rows
                ],
            )
            requested_order = tuple(order_by.get(dataset, ()))
            if any(column not in columns for column in requested_order):
                raise ValueError(f"{dataset}: ORDER BY references an unknown column")
            order_sql = (
                " ORDER BY "
                + ", ".join(quote_identifier(column) for column in requested_order)
                if requested_order
                else ""
            )
            sql = (
                f"SELECT {', '.join(quote_identifier(column) for column in columns)} "
                f"FROM {quote_identifier(dataset)}{order_sql};"
            )
            cursor = connection.execute(sql)
            observed_columns = [description[0] for description in cursor.description]
            queried[dataset] = [
                dict(zip(observed_columns, values, strict=True))
                for values in cursor.fetchall()
            ]
            if len(queried[dataset]) != len(rows):
                raise ValueError(f"{dataset}: SQLite query did not conserve row count")
            query_path = pending_queries / f"{dataset}.sql"
            write_text_exclusive(query_path, f"{sql}\n")
            query_specs.append((dataset, sql))
        connection.commit()
        connection.close()
        connection = None
        database_identity = file_identity(pending_db)
        for dataset, sql in query_specs:
            sources.append(
                {
                    "id": source_id(dataset),
                    "label": f"Executed SQLite snapshot query: {dataset}",
                    "path": f"queries/{dataset}.sql",
                    "query": {
                        "engine": "sqlite",
                        "description": (
                            f"Read reviewed table {dataset} from data/report.sqlite "
                            f"(sha256={database_identity['sha256']}, "
                            f"bytes={database_identity['size_bytes']})."
                        ),
                        "sql": sql,
                        "tables_used": [dataset],
                        "executed_at": query_executed_at,
                        "source_data_created_at": generated_at,
                        "metric_definitions": [
                            "All counts retain the unit and denominator encoded in their field names.",
                            "Primary strict edge support is at least 3 distinct canonical molecules.",
                            "Distance fields are diagnostics and never region inclusion rules.",
                        ],
                    },
                }
            )
        os.replace(pending_db, db_path)
        published_db = True
        os.replace(pending_queries, query_dir)
        published_queries = True
    except Exception as error:
        if connection is not None:
            connection.close()
        archive_dir = (
            output_root
            / "archive"
            / "failed_artifact_staging"
            / stage_id
        )
        archive_dir.mkdir(parents=True, exist_ok=False)
        candidates = (
            (db_path if published_db else pending_db, "report.sqlite"),
            (query_dir if published_queries else pending_queries, "queries"),
        )
        for candidate, archive_name in candidates:
            if candidate.exists():
                os.replace(candidate, archive_dir / archive_name)
        raise ValueError(
            f"SQLite artifact staging failed; preserved evidence at {archive_dir}"
        ) from error
    return queried, sources


def build_artifact(data: Mapping[str, Any], data_path: Path, output_path: Path) -> dict[str, Any]:
    validate_report_data_package(data, data_path)

    relative_source = os.path.relpath(data_path.resolve(), output_path.parent.resolve())
    if relative_source.startswith("..") or os.path.isabs(relative_source):
        raise ValueError("report data must be located inside the artifact root")

    datasets = [dict(row) for row in data["datasets"]]
    chromosomes = [dict(row) for row in data["chromosomes"]]
    k_bands = [dict(row) for row in data["k_bands"]]
    span_bands = [dict(row) for row in data["span_bands"]]
    exact_k = [dict(row) for row in data["exact_k"]]
    threshold_rows = [dict(row) for row in data["threshold_sensitivity"]]
    topology_completion = [dict(row) for row in data["topology_completion"]]
    W_by_dataset = {row["dataset"]: row["tree_eligible_W"] for row in datasets}
    k_band_order = {
        "2": 1,
        "3": 2,
        "4": 3,
        "5": 4,
        "6–8": 5,
        "9–12": 6,
        ">12": 7,
    }
    span_band_order = {
        "≤10 kb": 1,
        "10–25 kb": 2,
        "25–50 kb": 3,
        "50–100 kb": 4,
        ">100 kb": 5,
    }
    for rows in (k_bands, span_bands):
        for row in rows:
            denominator = W_by_dataset[row["dataset"]]
            row["W_denominator"] = denominator
            row["fraction_of_W"] = row["W"] / denominator if denominator else None
    for row in k_bands:
        row["k_band_order"] = k_band_order[row["k_band"]]
    for row in span_bands:
        row["span_band_order"] = span_band_order[row["span_band"]]

    total_s = sum(row["candidate_loci_S"] for row in datasets)
    total_active = sum(row["active_unique_loci"] for row in datasets)
    total_linked_loci = sum(row["unique_loci_in_any_W"] for row in datasets)
    total_components = sum(row["all_components"] for row in datasets)
    total_singletons = sum(row["singleton_k1_components"] for row in datasets)
    total_w = sum(row["tree_eligible_W"] for row in datasets)
    total_edges = sum(row["retained_direct_edges"] for row in datasets)
    total_w_span50 = sum(row["W_span_gt_50kb"] for row in datasets)
    total_w_gap50 = sum(row["W_adjacent_gap_gt_50kb"] for row in datasets)
    total_long_edges = sum(row["direct_edges_gt_50kb"] for row in datasets)
    total_w_with_long_edge = sum(
        row["W_with_direct_edge_gt_50kb"] for row in datasets
    )
    total_w_changed_by_cut = sum(
        row["W_partition_changed_if_50kb_cut"] for row in datasets
    )
    total_w_lost_by_cut = sum(row["W_fully_lost_if_50kb_cut"] for row in datasets)
    total_counterfactual_w = sum(
        row["counterfactual_W_after_50kb_cut"] for row in datasets
    )
    total_counterfactual_w_delta = sum(
        row["counterfactual_W_delta_after_50kb_cut"] for row in datasets
    )
    total_memberships_lost_by_cut = sum(
        row["W_memberships_lost_to_singletons_if_50kb_cut"] for row in datasets
    )

    headline = [
        {
            "technical_datasets": 7,
            "biological_cell_lines": 6,
            "autosome_dataset_rows": 154,
            "dataset_summed_candidate_locus_records": total_s,
            "dataset_summed_active_locus_records": total_active,
            "dataset_summed_locus_records_in_W": total_linked_loci,
            "all_components": total_components,
            "singleton_k1_components": total_singletons,
            "tree_eligible_W": total_w,
            "retained_direct_edges": total_edges,
            "W_span_gt_50kb": total_w_span50,
            "W_adjacent_gap_gt_50kb": total_w_gap50,
            "direct_edges_gt_50kb": total_long_edges,
            "W_partition_changed_if_50kb_cut": total_w_changed_by_cut,
            "counterfactual_W_after_50kb_cut": total_counterfactual_w,
        }
    ]
    locus_stage = [
        {"dataset": row["dataset"], "stage": stage, "loci": row[field_name]}
        for row in datasets
        for stage, field_name in (
            ("候選 S", "candidate_loci_S"),
            ("Active 去重", "active_unique_loci"),
            ("進入 ≥1 W", "unique_loci_in_any_W"),
        )
    ]
    molecule_row_outcome = [
        {
            "dataset": row["dataset"],
            "outcome": outcome,
            "molecule_rows": row[field_name],
            "canonical_molecule_rows": row["canonical_molecule_rows"],
            "fraction_of_canonical_rows": (
                row[field_name] / row["canonical_molecule_rows"]
                if row["canonical_molecule_rows"]
                else None
            ),
        }
        for row in datasets
        for outcome, field_name in (
            ("Primary exact PS×HP", "primary_known_PS_molecule_rows"),
            ("排除：missing PS", "excluded_missing_PS_rows"),
            ("排除：nonprimary HP", "excluded_nonprimary_HP_rows"),
        )
    ]
    component_outcome = [
        {
            "dataset": row["dataset"],
            "outcome": outcome,
            "components": row[field_name],
            "all_components": row["all_components"],
            "fraction_of_components": (
                row[field_name] / row["all_components"]
                if row["all_components"]
                else None
            ),
        }
        for row in datasets
        for outcome, field_name in (
            ("k=1 不建樹", "singleton_k1_components"),
            ("k≥2 W", "tree_eligible_W"),
        )
    ]
    distance_w_qc = [
        {
            "dataset": row["dataset"],
            "criterion": criterion,
            "fraction_of_W": pct_fraction(row[field_name]),
            "regions": row[count_name],
            "W_denominator": row["tree_eligible_W"],
        }
        for row in datasets
        for criterion, field_name, count_name in (
            ("總 span >50 kb", "W_span_gt_50kb_pct", "W_span_gt_50kb"),
            (
                "最大相鄰 gap >50 kb",
                "W_adjacent_gap_gt_50kb_pct",
                "W_adjacent_gap_gt_50kb",
            ),
        )
    ]
    direct_long_edge = [
        {
            "dataset": row["dataset"],
            "fraction_of_edges": pct_fraction(row["direct_edges_gt_50kb_pct"]),
            "direct_edges_gt_50kb": row["direct_edges_gt_50kb"],
            "retained_direct_edges": row["retained_direct_edges"],
            "max_direct_edge_distance_bp": row["max_direct_edge_distance_bp"],
        }
        for row in datasets
    ]
    distance_cut_impact = [
        {
            "dataset": row["dataset"],
            "criterion": criterion,
            "fraction_of_W": (
                count / row["tree_eligible_W"]
                if row["tree_eligible_W"]
                else None
            ),
            "regions": count,
            "W_denominator": row["tree_eligible_W"],
        }
        for row in datasets
        for criterion, count in (
            (
                "含 >50 kb 直接邊",
                row["W_with_direct_edge_gt_50kb"],
            ),
            (
                "硬切後 partition 改變",
                row["W_partition_changed_if_50kb_cut"],
            ),
            (
                "硬切後完全失去 W",
                row["W_fully_lost_if_50kb_cut"],
            ),
        )
    ]
    long_edge_support = [
        {
            "dataset": row["dataset"],
            "support_class": support_class,
            "direct_edges_gt_50kb": row[field_name],
        }
        for row in datasets
        for support_class, field_name in (
            ("support=3", "direct_edges_gt_50kb_support_3"),
            ("support=4", "direct_edges_gt_50kb_support_4"),
            ("support≥5", "direct_edges_gt_50kb_support_ge_5"),
        )
    ]
    long_edge_state = [
        {
            "dataset": row["dataset"],
            "state_class": state_class,
            "direct_edges_gt_50kb": row[field_name],
        }
        for row in datasets
        for state_class, field_name in (
            ("僅 RR callability", "direct_edges_gt_50kb_RR_only"),
            (
                "含 ALT 資訊（RA/AR/AA）",
                "direct_edges_gt_50kb_alt_informative",
            ),
        )
    ]
    definitions = [
        {
            "order": 1,
            "term": "S",
            "unit": "physical locus",
            "definition": (
                "chr1–22 的 upstream autosomal biallelic LongPhase-S PASS sSNV "
                "target；只在各 technical dataset 內去重。"
            ),
        },
        {
            "order": 2,
            "term": "exact PS",
            "unit": "phase-block identifier",
            "definition": (
                "同一 dataset/染色體內的局部 phase 座標系；不同 PS 的 HP "
                "orientation 未知，PS 也不是 clone boundary 或 phasing 正確性的保證。"
            ),
        },
        {
            "order": 3,
            "term": "Primary HP family",
            "unit": "read-side haplotype",
            "definition": (
                "raw 1/1-1/1-2→HP1；2/2-1/2-2→HP2。HP3、HP4、unphased "
                "與 missing PS 不建立 primary graph。"
            ),
        },
        {
            "order": 4,
            "term": "Container H",
            "unit": "dataset×chromosome×PS×HP",
            "definition": (
                "dataset × chromosome × exact nonmissing PS × primary HP family；"
                "這四欄全相同才可共同連邊。"
            ),
        },
        {
            "order": 5,
            "term": "Canonical molecule",
            "unit": "distinct read molecule",
            "definition": (
                "SHA256(dataset, RG-or-'.', QNAME)；只納入 mapped primary alignment，"
                "排除 secondary/supplementary/QC-fail/duplicate，MAPQ≥20。"
            ),
        },
        {
            "order": 6,
            "term": "Fixed R/A call",
            "unit": "molecule×endpoint",
            "definition": (
                "BQ≥20 且鹼基等於 REF(R) 或 ALT(A)。X/L/O/D/S 代表缺失、"
                "低品質、其他鹼基、deletion、ref-skip，均不替 endpoint edge 投票。"
            ),
        },
        {
            "order": 7,
            "term": "Node",
            "unit": "locus membership",
            "definition": "容器中至少由一個 canonical molecule fixed-observe 為 R 或 A 的 sSNV。",
        },
        {
            "order": 8,
            "term": "Direct edge",
            "unit": "endpoint pair",
            "definition": (
                "至少 3 個 distinct canonical molecules 同時 fixed-observe 兩個 "
                "endpoint；每 molecule 對每 pair 最多一票。"
            ),
        },
        {
            "order": 9,
            "term": "W",
            "unit": "read-chain candidate region",
            "definition": (
                "合格直接邊的 maximal connected component 且 k≥2；只表示可共同分析，"
                "不表示一條 full-span read、同一細胞或真實 clone。"
            ),
        },
        {
            "order": 10,
            "term": "k",
            "unit": "sSNV memberships per W",
            "definition": "一個 W 內的節點數；不是 clone 數，也不是樹候選數 T。",
        },
        {
            "order": 11,
            "term": "k=1 abstain",
            "unit": "audit component",
            "definition": (
                "在指定 threshold 下沒有合格 edge 的單點；保留稽核但排除建樹，"
                "不等於該 physical locus 全域沒有 read linkage。"
            ),
        },
        {
            "order": 12,
            "term": "Membership",
            "unit": "locus×PS×HP",
            "definition": (
                "同一 physical locus 可在不同 exact PS/HP 容器重複出現，"
                "因此不可與 unique loci 混為同一分母。"
            ),
        },
        {
            "order": 13,
            "term": "Span",
            "unit": "bp per W",
            "definition": "max(position)−min(position) 的顯示 envelope；不決定成員資格。",
        },
        {
            "order": 14,
            "term": "Max adjacent gap",
            "unit": "bp per W",
            "definition": "W 內排序相鄰 membership 的最大座標差；僅為距離 QC。",
        },
        {
            "order": 15,
            "term": "Direct-edge distance",
            "unit": "bp per edge",
            "definition": "同一合格 read-supported endpoint pair 的座標距離；不作 hard cutoff。",
        },
        {
            "order": 16,
            "term": "RR-only edge",
            "unit": "edge evidence class",
            "definition": (
                "只證明兩端在 reference-state molecules 上共同 callable/linkable，"
                "不能解讀成 somatic ALT 共現。"
            ),
        },
        {
            "order": 17,
            "term": "ALT-informative edge",
            "unit": "edge evidence class",
            "definition": "RA、AR 或 AA 至少一票；仍不是 clone truth 或祖先—後代方向。",
        },
        {
            "order": 18,
            "term": "Tree candidate",
            "unit": "T",
            "definition": (
                "後續 Boolean-hypercube/Steiner 階段的候選 mutation-state tree；"
                "本報告尚未統計或宣稱唯一拓撲。"
            ),
        },
    ]
    qa_rows = [
        {"check": check, "passed": "PASS" if passed is True else "FAIL"}
        for check, passed in data["checks"].items()
    ]
    for row in chromosomes:
        chrom_number = int(row["chrom"].removeprefix("chr"))
        row["dataset_chr_order"] = f"{row['dataset']} · chr{chrom_number:02d}"
        row["all_pass_display"] = "PASS" if row["all_pass"] is True else "FAIL"
    for row in threshold_rows:
        row["dataset_threshold_order"] = (
            f"{row['dataset']} · {row['minimum_distinct_molecules']:02d}"
        )
    for row in exact_k:
        row["dataset_k_order"] = f"{row['dataset']} · k={row['k']:04d}"

    generated_at = data.get("created_at_utc") or utc_now()
    raw_source = {
        "id": RAW_SOURCE_ID,
        "label": "7 datasets × chr1–22 strict-region independently recomputed package",
        "path": relative_source,
    }
    raw_snapshot = {
        "headline": headline,
        "dataset_summary": datasets,
        "locus_stage": locus_stage,
        "molecule_row_outcome": molecule_row_outcome,
        "component_outcome": component_outcome,
        "k_bands": k_bands,
        "exact_k": exact_k,
        "span_bands": span_bands,
        "distance_w_qc": distance_w_qc,
        "direct_long_edge": direct_long_edge,
        "distance_cut_impact": distance_cut_impact,
        "long_edge_support": long_edge_support,
        "long_edge_state": long_edge_state,
        "threshold_sensitivity": threshold_rows,
        "chromosome_detail": chromosomes,
        "definitions": definitions,
        "qa_checks": qa_rows,
        "topology_completion": topology_completion,
    }
    order_by = {
        "headline": ("technical_datasets",),
        "dataset_summary": ("dataset",),
        "locus_stage": ("dataset", "stage"),
        "molecule_row_outcome": ("dataset", "outcome"),
        "component_outcome": ("dataset", "outcome"),
        "k_bands": ("dataset", "k_band_order"),
        "exact_k": ("dataset_k_order",),
        "span_bands": ("dataset", "span_band_order"),
        "distance_w_qc": ("dataset", "criterion"),
        "direct_long_edge": ("dataset",),
        "distance_cut_impact": ("dataset", "criterion"),
        "long_edge_support": ("dataset", "support_class"),
        "long_edge_state": ("dataset", "state_class"),
        "threshold_sensitivity": (
            "dataset",
            "minimum_distinct_molecules",
        ),
        "chromosome_detail": ("dataset_chr_order",),
        "definitions": ("order",),
        "qa_checks": ("check",),
        "topology_completion": ("dataset",),
    }
    queried_snapshot, query_sources = stage_and_query_snapshot(
        output_root=output_path.parent,
        raw_datasets=raw_snapshot,
        order_by=order_by,
        generated_at=generated_at,
    )
    sources = [raw_source, *query_sources]

    charts = [
        chart(
            chart_id="molecule_row_outcome_chart",
            title="各資料集 canonical molecule rows 的 PS×HP routing",
            subtitle=(
                "三類比例以各資料集全部 canonical molecule rows 為共同分母，"
                "合計 100%；missing PS 與 nonprimary HP 不進 primary graph。"
            ),
            dataset="molecule_row_outcome",
            chart_type="bar",
            x_field="dataset",
            x_type="nominal",
            x_label="技術資料集",
            y_field="fraction_of_canonical_rows",
            y_label="Molecule-row 比例",
            color_field="outcome",
            color_label="Routing 結果",
            value_format="percent",
            grouping="grouped",
            intent="composition",
            question="每個資料集有多少 canonical molecule rows 真正進入 exact PS×primary HP graph？",
            rationale="Grouped percentage bars preserve the common canonical-molecule-row denominator.",
            tooltip=(
                {"field": "molecule_rows", "type": "quantitative", "label": "Molecule rows"},
                {
                    "field": "canonical_molecule_rows",
                    "type": "quantitative",
                    "label": "All canonical molecule rows",
                },
            ),
        ),
        chart(
            chart_id="locus_stage_chart",
            title="各資料集的物理 sSNV 位點保留情況",
            subtitle="三個階段都以 unique physical loci 為單位，不把 region 數混入漏斗。",
            dataset="locus_stage",
            chart_type="bar",
            x_field="dataset",
            x_type="nominal",
            x_label="技術資料集",
            y_field="loci",
            y_label="去重物理位點",
            color_field="stage",
            color_label="位點階段",
            grouping="grouped",
            intent="compare",
            question="各資料集由候選位點到具 read-linkage 的 physical loci 保留多少？",
            rationale="Grouped bars preserve the common physical-locus denominator across stages.",
        ),
        chart(
            chart_id="component_outcome_chart",
            title="各資料集的 graph component 分類比例",
            subtitle=(
                "兩類比例以各資料集全部 components 為共同分母；"
                "每個 component 只能屬於 k=1 abstain 或 k≥2 W，合計 100%。"
            ),
            dataset="component_outcome",
            chart_type="bar",
            x_field="dataset",
            x_type="nominal",
            x_label="技術資料集",
            y_field="fraction_of_components",
            y_label="Component 比例",
            color_field="outcome",
            color_label="分類結果",
            value_format="percent",
            grouping="grouped",
            intent="composition",
            question="各資料集的 graph components 有多少比例具備多位點 read linkage？",
            rationale="Grouped percentage bars expose the conserved component partition.",
            tooltip=(
                {"field": "components", "type": "quantitative", "label": "Components"},
                {
                    "field": "all_components",
                    "type": "quantitative",
                    "label": "All components",
                },
            ),
        ),
        chart(
            chart_id="k_distribution_chart",
            title="W 內 sSNV 數（k）分布",
            subtitle=(
                "各 k 分層以該資料集全部 W 為共同分母且合計 100%；"
                "k 是一個 W 內的 sSNV membership 數，不代表 clone 數。"
            ),
            dataset="k_bands",
            chart_type="bar",
            x_field="dataset",
            x_type="nominal",
            x_label="技術資料集",
            y_field="fraction_of_W",
            y_label="W 比例",
            color_field="k_band",
            color_label="k 分層",
            value_format="percent",
            grouping="grouped",
            intent="composition",
            question="各資料集的 W 主要由幾個 sSNV 構成？",
            rationale="Grouped shares make the k composition comparable across unequal W totals.",
            tooltip=(
                {"field": "W", "type": "quantitative", "label": "W"},
                {"field": "W_denominator", "type": "quantitative", "label": "All W"},
            ),
        ),
        chart(
            chart_id="span_distribution_chart",
            title="W 的座標跨度分布",
            subtitle=(
                "五個互斥 span 分層以該資料集全部 W 為共同分母且合計 100%；"
                "span 是座標 envelope，不參與建立或切割 W。"
            ),
            dataset="span_bands",
            chart_type="bar",
            x_field="dataset",
            x_type="nominal",
            x_label="技術資料集",
            y_field="fraction_of_W",
            y_label="W 比例",
            color_field="span_band",
            color_label="座標跨度",
            value_format="percent",
            grouping="grouped",
            intent="composition",
            question="各資料集的 W coordinate span 如何分布？",
            rationale=(
                "Ordered grouped shares preserve the complete W denominator and reveal "
                "long chained regions without treating 50 kb as a boundary."
            ),
            tooltip=(
                {"field": "W", "type": "quantitative", "label": "W"},
                {"field": "W_denominator", "type": "quantitative", "label": "All W"},
            ),
        ),
        chart(
            chart_id="distance_w_qc_chart",
            title="50 kb 的 region-level 距離 QC",
            subtitle="兩個百分比的共同分母都是各資料集的全部 W；total span 與單一 adjacent gap 不可混稱。",
            dataset="distance_w_qc",
            chart_type="bar",
            x_field="dataset",
            x_type="nominal",
            x_label="技術資料集",
            y_field="fraction_of_W",
            y_label="占全部 W 的比例",
            color_field="criterion",
            color_label="距離條件",
            value_format="percent",
            grouping="grouped",
            intent="compare",
            question="跨 50 kb 的 W 是由整體 chain 延伸，還是含單一大座標 gap？",
            rationale="Both measures share W as denominator, so grouped percentages are valid.",
            tooltip=(
                {"field": "regions", "type": "quantitative", "label": "Regions"},
                {"field": "W_denominator", "type": "quantitative", "label": "All W"},
            ),
        ),
        chart(
            chart_id="direct_long_edge_chart",
            title="端點距離 >50 kb 的直接 read edges",
            subtitle="分母是各資料集全部 retained direct edges；這是 edge-level QC，不能與 W 百分比相加。",
            dataset="direct_long_edge",
            chart_type="bar",
            x_field="dataset",
            x_type="nominal",
            x_label="技術資料集",
            y_field="fraction_of_edges",
            y_label="占全部 retained edges 的比例",
            value_format="percent",
            grouping="single",
            intent="compare",
            question="實際有多少合格直接 read edge 的端點距離超過 50 kb？",
            rationale="A separate chart preserves the edge-level denominator.",
            tooltip=(
                {
                    "field": "direct_edges_gt_50kb",
                    "type": "quantitative",
                    "label": "Edges >50 kb",
                },
                {
                    "field": "retained_direct_edges",
                    "type": "quantitative",
                    "label": "All retained edges",
                },
                {
                    "field": "max_direct_edge_distance_bp",
                    "type": "quantitative",
                    "label": "Maximum distance (bp)",
                },
            ),
        ),
        chart(
            chart_id="distance_cut_impact_chart",
            title="50 kb 硬切的反事實影響",
            subtitle=(
                "逐一移除 retained direct edges>50 kb 後重建 components；"
                "三項共同分母為原始 W。"
            ),
            dataset="distance_cut_impact",
            chart_type="bar",
            x_field="dataset",
            x_type="nominal",
            x_label="技術資料集",
            y_field="fraction_of_W",
            y_label="占原始 W 的比例",
            color_field="criterion",
            color_label="反事實結果",
            value_format="percent",
            grouping="grouped",
            intent="compare",
            question="若恢復 50 kb hard cutoff，多少原始 W 的 read-linkage partition 會改變？",
            rationale=(
                "A direct removal counterfactual tests the practical consequence of the "
                "distance rule instead of inferring it from span alone."
            ),
            tooltip=(
                {
                    "field": "regions",
                    "type": "quantitative",
                    "label": "Regions meeting criterion",
                },
                {
                    "field": "W_denominator",
                    "type": "quantitative",
                    "label": "All original W",
                },
            ),
        ),
        chart(
            chart_id="long_edge_support_chart",
            title=">50 kb 直接邊的 molecule-support 分層",
            subtitle="只分解已跨過 50 kb 的直接邊；support 是 distinct canonical molecule 數。",
            dataset="long_edge_support",
            chart_type="bar",
            x_field="dataset",
            x_type="nominal",
            x_label="技術資料集",
            y_field="direct_edges_gt_50kb",
            y_label=">50 kb 直接邊",
            color_field="support_class",
            color_label="Molecule 支持數",
            grouping="grouped",
            intent="composition",
            question="跨 50 kb 的直接邊主要落在哪一個 molecule-support 層級？",
            rationale="Grouped counts expose whether long edges concentrate near the operational threshold.",
        ),
        chart(
            chart_id="long_edge_state_chart",
            title=">50 kb 直接邊的 endpoint-state 證據",
            subtitle="ALT-informative 表示 RA、AR 或 AA 至少一票；RR-only 只支持 callable linkage，不是 somatic ALT 共現。",
            dataset="long_edge_state",
            chart_type="bar",
            x_field="dataset",
            x_type="nominal",
            x_label="技術資料集",
            y_field="direct_edges_gt_50kb",
            y_label=">50 kb 直接邊",
            color_field="state_class",
            color_label="Endpoint-state 證據",
            grouping="grouped",
            intent="composition",
            question="跨 50 kb 的直接邊有多少含 somatic ALT pattern，而非只有 RR？",
            rationale="Separating RR-only from ALT-informative edges prevents callability from being overread as mutation co-occurrence.",
        ),
        chart(
            chart_id="threshold_sensitivity_chart",
            title="不同 molecule-support 門檻下的可建樹 W",
            subtitle=(
                "Threshold=3 是 primary setting；edge 數隨 threshold 單調不增，"
                "但 W 可因大 component 裂成多個 k≥2 component，數量不保證單調。"
                "Threshold=3 由 row-level 獨立重算；1/2/5 由目前 22 份 chromosome "
                "receipts 重加總並與 dataset summary 對帳。"
            ),
            dataset="threshold_sensitivity",
            chart_type="line",
            x_field="minimum_distinct_molecules",
            x_type="quantitative",
            x_label="最少 distinct molecules",
            y_field="multisite_components",
            y_label="可建樹多位點 components",
            color_field="dataset",
            color_label="技術資料集",
            intent="trend",
            question="W 數量對 read-support threshold 有多敏感？",
            rationale="A line per dataset exposes threshold dependence without implying calibration truth.",
        ),
    ]

    tables = [
        {
            "id": "dataset_funnel_audit",
            "title": "各資料集的 locus、membership 與 component 對帳",
            "subtitle": "各欄位保留自身單位與分母；HCC1395 與 HCC1395_DORADO 是兩個 technical datasets。",
            "dataset": "dataset_summary",
            "sourceId": source_id("dataset_summary"),
            "defaultSort": {"field": "tree_eligible_W", "direction": "desc"},
            "density": "compact",
            "layout": "full",
            "columns": [
                field("dataset", "Dataset", kind="text"),
                field("canonical_molecule_rows", "Canonical molecule rows"),
                field("primary_known_PS_molecule_rows", "Primary exact-PS×HP rows"),
                field("excluded_missing_PS_rows", "Excluded missing-PS rows"),
                field("excluded_nonprimary_HP_rows", "Excluded nonprimary-HP rows"),
                field("candidate_loci_S", "S loci"),
                field("active_unique_loci", "Active loci"),
                field("active_unique_pct_of_S", "Active / S (%)"),
                field("unique_loci_in_any_W", "Loci in ≥1 W"),
                field("unique_loci_in_W_pct_of_active", "In W / active (%)"),
                field("exact_chromosome_PS", "Exact chr×PS"),
                field("exact_PS_HP_containers", "PS×HP containers"),
                field("active_node_memberships", "Active memberships"),
                field("singleton_memberships", "Singleton memberships"),
                field("W_memberships", "W memberships"),
                field("all_components", "All components"),
                field("singleton_k1_components", "k=1 abstain"),
                field("singleton_pct_of_all_components", "k=1 / all (%)"),
                field("tree_eligible_W", "W (k≥2)"),
                field("W_pct_of_all_components", "W / all (%)"),
                field("HP1_W", "HP1 W"),
                field("HP2_W", "HP2 W"),
            ],
        },
        {
            "id": "dataset_distance_audit",
            "title": "各資料集的 k、span 與 direct-edge 明細",
            "subtitle": "Long-region、large-gap 與 direct long-edge 三種指標分開列出。",
            "dataset": "dataset_summary",
            "sourceId": source_id("dataset_summary"),
            "defaultSort": {"field": "W_span_gt_50kb", "direction": "desc"},
            "density": "compact",
            "layout": "full",
            "columns": [
                field("dataset", "Dataset", kind="text"),
                field("tree_eligible_W", "All W"),
                field("mean_W_k", "Mean k"),
                field("median_W_k", "Median k"),
                field("p95_W_k", "P95 k"),
                field("max_W_k", "Max k"),
                field("W_k_gt12", "W with k>12"),
                field("median_W_span_bp", "Median span (bp)"),
                field("p95_W_span_bp", "P95 span (bp)"),
                field("max_W_span_bp", "Max span (bp)"),
                field("W_span_gt_50kb", "W span>50 kb"),
                field("W_span_gt_50kb_pct", "Span>50 / W (%)"),
                field("W_adjacent_gap_gt_50kb", "W max-gap>50 kb"),
                field("W_adjacent_gap_gt_50kb_pct", "Max-gap>50 / W (%)"),
                field(
                    "W_span_gt_50kb_without_adjacent_gap_gt_50kb",
                    "Long-span chained W",
                ),
                field("retained_direct_edges", "Retained edges"),
                field("direct_edges_gt_50kb", "Direct edges>50 kb"),
                field("direct_edges_gt_50kb_pct", "Long edges / all (%)"),
                field("direct_edges_gt_50kb_support_3", "Long edge support=3"),
                field("direct_edges_gt_50kb_support_4", "Long edge support=4"),
                field("direct_edges_gt_50kb_support_ge_5", "Long edge support≥5"),
                field("direct_edges_gt_50kb_RR_only", "Long edge RR-only"),
                field(
                    "direct_edges_gt_50kb_alt_informative",
                    "Long edge ALT-informative",
                ),
                field("W_with_direct_edge_gt_50kb", "W containing long edge"),
                field(
                    "W_partition_changed_if_50kb_cut",
                    "W partition changed if cut",
                ),
                field(
                    "W_partition_changed_if_50kb_cut_pct",
                    "Changed / original W (%)",
                ),
                field("W_fully_lost_if_50kb_cut", "W fully lost if cut"),
                field(
                    "W_reduced_to_one_W_if_50kb_cut",
                    "Changed but remains one W",
                ),
                field(
                    "W_split_to_multiple_W_if_50kb_cut",
                    "Split into multiple W",
                ),
                field(
                    "counterfactual_W_after_50kb_cut",
                    "Counterfactual W after cut",
                ),
                field(
                    "counterfactual_W_delta_after_50kb_cut",
                    "Counterfactual ΔW",
                ),
                field(
                    "W_memberships_lost_to_singletons_if_50kb_cut",
                    "W memberships lost to k=1",
                ),
                field("max_direct_edge_distance_bp", "Max direct distance (bp)"),
            ],
        },
        {
            "id": "exact_k_audit",
            "title": "各資料集的完整 exact-k 分布",
            "subtitle": "每列是一個實際出現的 dataset × k；W 與 pct_of_W 的分母均為該資料集全部 W。",
            "dataset": "exact_k",
            "sourceId": source_id("exact_k"),
            "defaultSort": {"field": "dataset_k_order", "direction": "asc"},
            "density": "compact",
            "columns": [
                field("dataset_k_order", "Dataset · k", kind="text"),
                field("dataset", "Dataset", kind="text"),
                field("k", "sSNV memberships (k)"),
                field("W", "W count"),
                field("pct_of_W", "占該資料集 W (%)"),
            ],
        },
        {
            "id": "threshold_audit",
            "title": "Read-support threshold 敏感度明細",
            "subtitle": (
                "每列是一個 dataset × minimum distinct molecules 設定；"
                "threshold=3 另由原始列重算，其餘門檻由 22 份 chromosome receipts 重新加總並與 summary 對帳。"
            ),
            "dataset": "threshold_sensitivity",
            "sourceId": source_id("threshold_sensitivity"),
            "defaultSort": {"field": "dataset_threshold_order", "direction": "asc"},
            "density": "compact",
            "columns": [
                field("dataset_threshold_order", "Dataset · threshold", kind="text"),
                field("dataset", "Dataset", kind="text"),
                field("minimum_distinct_molecules", "Minimum molecules"),
                field("all_components", "All components"),
                field("singleton_k1_components", "k=1 abstain"),
                field("multisite_components", "k≥2 components"),
                field("retained_direct_edges", "Retained direct edges"),
                field("validation_basis", "Validation basis", kind="text"),
            ],
        },
        {
            "id": "chromosome_audit",
            "title": "全部 154 列 dataset × chromosome 明細",
            "subtitle": "完整列出 7 technical datasets × chr1–22；dataset totals 可由本表加總重算。",
            "dataset": "chromosome_detail",
            "sourceId": source_id("chromosome_detail"),
            "defaultSort": {"field": "dataset_chr_order", "direction": "asc"},
            "density": "compact",
            "layout": "full",
            "columns": [
                field("dataset_chr_order", "Dataset · chromosome", kind="text"),
                field("candidate_loci_S", "S loci"),
                field("active_unique_loci", "Active loci"),
                field("unique_loci_in_any_W", "Loci in ≥1 W"),
                field("active_node_memberships", "Memberships"),
                field("all_components", "All components"),
                field("singleton_k1_components", "k=1 abstain"),
                field("tree_eligible_W", "W"),
                field("W_memberships", "W memberships"),
                field("retained_direct_edges", "Direct edges"),
                field("direct_edges_gt_50kb", "Edges>50 kb"),
                field("W_span_gt_50kb", "W span>50 kb"),
                field("W_adjacent_gap_gt_50kb", "W max-gap>50 kb"),
                field(
                    "W_partition_changed_if_50kb_cut",
                    "W changed if 50 kb cut",
                ),
                field(
                    "counterfactual_W_after_50kb_cut",
                    "W after 50 kb cut",
                ),
                field(
                    "W_memberships_lost_to_singletons_if_50kb_cut",
                    "Memberships lost if cut",
                ),
                field("W_k_gt12", "W k>12"),
                field("median_W_k", "Median k"),
                field("median_W_span_bp", "Median span (bp)"),
                field("max_W_span_bp", "Max span (bp)"),
                field("max_direct_edge_distance_bp", "Max edge distance (bp)"),
                field("all_pass_display", "Validation", kind="text"),
            ],
        },
        {
            "id": "definition_table",
            "title": "名詞、單位與證據邊界",
            "subtitle": "同名數字若單位不同，不可直接相除或串成同一漏斗。",
            "dataset": "definitions",
            "sourceId": source_id("definitions"),
            "defaultSort": {"field": "order", "direction": "asc"},
            "columns": [
                field("order", "#"),
                field("term", "Term", kind="text"),
                field("unit", "Unit", kind="text"),
                field("definition", "Definition", kind="text"),
            ],
        },
        {
            "id": "qa_table",
            "title": "資料品質與守恆檢查",
            "subtitle": "所有檢查都必須為 true，artifact 才能標記為 ready。",
            "dataset": "qa_checks",
            "sourceId": source_id("qa_checks"),
            "defaultSort": {"field": "check", "direction": "asc"},
            "columns": [
                field("check", "Check", kind="text"),
                field("passed", "Status", kind="text"),
            ],
        },
        {
            "id": "topology_completion_table",
            "title": "各資料集最新完成層級",
            "subtitle": (
                "L1 strict W 已完成 7/7；新版 production strict directed topology "
                "為 0/7，clone／精確 parent–child／融合樹均為 0/7。"
            ),
            "dataset": "topology_completion",
            "sourceId": source_id("topology_completion"),
            "defaultSort": {"field": "dataset", "direction": "asc"},
            "columns": [
                field("dataset", "Dataset", kind="text"),
                field("strict_W", "Strict W"),
                field("strict_linkage_status", "L1 read-linkage", kind="text"),
                field(
                    "exploratory_exact_ps_topology_status",
                    "Earlier exploratory topology",
                    kind="text",
                ),
                field(
                    "strict_directed_topology_status",
                    "L3 production strict topology",
                    kind="text",
                ),
                field(
                    "cellular_clone_count_status",
                    "Clone count",
                    kind="text",
                ),
                field(
                    "exact_cellular_parent_child_status",
                    "Exact cellular parent→child",
                    kind="text",
                ),
                field(
                    "cross_hp_or_technical_fusion_status",
                    "Cross-HP / technical fusion",
                    kind="text",
                ),
            ],
        },
    ]

    summary_body = (
        "## 結論\n\n"
        f"- **50 kb 不應作唯一、無條件的 hard cutoff。** 正式 W 完全由 exact PS×HP "
        f"內、通過目前 filters 的直接 read evidence 建立；"
        f"7 個 technical datasets 合計有 **{fmt_int(total_w)} 個 W**，其中 "
        f"**{fmt_int(total_w_span50)} 個**總 span>50 kb，但只有 "
        f"**{fmt_int(total_w_gap50)} 個**含單一 adjacent gap>50 kb。\n"
        f"- 實際保留的 **{fmt_int(total_edges)} 條 direct edges** 中，"
        f"**{fmt_int(total_long_edges)} 條**端點距離>50 kb；這些證據若被全域 cutoff 丟棄，"
        f"涉及 **{fmt_int(total_w_with_long_edge)} 個 W**。反事實重建顯示 "
        f"**{fmt_int(total_w_changed_by_cut)} 個 W** 的 linkage partition 會改變，"
        f"**{fmt_int(total_w_lost_by_cut)} 個**完全失去 W 資格；跨資料集 W 加總由 "
        f"**{fmt_int(total_w)}** 變為 **{fmt_int(total_counterfactual_w)}**"
        f"（Δ={total_counterfactual_w_delta:+,}），並有 "
        f"**{fmt_int(total_memberships_lost_by_cut)} 個 memberships** 退回 k=1。\n"
        "- 本報告是 **7 個資料產品、6 個 biological cell lines**；"
        "HCC1395 與 HCC1395_DORADO 不可當成兩個獨立生物樣本。\n"
        "- 完成層級：**L1 strict read-linkage 為 7/7；新版 production strict directed "
        "topology 為 0/7；clone 數、精確 cellular parent–child edge 與完整融合樹均為 "
        "0/7。** HCC1395 曾有較早的 exploratory exact-PS topology 技術觀察，但其 "
        "upstream receipt 為 PARTIAL、`validation_evidence_eligible=false`，且不是目前 "
        "85,621 個 strict W 的 production topology 結果。"
    )

    blocks: list[dict[str, Any]] = [
        {
            "id": "title",
            "type": "markdown",
            "body": (
                f"# {TITLE}\n\n"
                "Task Type B comprehensive validation · 7 technical datasets × chr1–22 · "
                "primary edge support ≥3 distinct canonical molecules · "
                "跨資料集 headline 為 per-dataset records 加總，並非跨樣本去重 loci"
            ),
        },
        {
            "id": "summary_text",
            "type": "markdown",
            "body": summary_body,
            "sourceId": source_id("dataset_summary"),
        },
        {
            "id": "headline_metrics",
            "type": "metric-strip",
            "cardIds": [
                "scope_card",
                "locus_card",
                "component_card",
                "edge_card",
                "distance_card",
                "counterfactual_card",
            ],
        },
        {
            "id": "completion_boundary_heading",
            "type": "markdown",
            "body": (
                "## 最新完成邊界：不要把 W 當成已完成的演化樹\n\n"
                "**已完成的是 L1：** 7/7 datasets 的 exact PS×HP strict read-linkage "
                "region。**尚未完成的是 L3/L4：** production strict directed topology "
                "0/7；clone 數、精確 cellular parent→child、跨 HP 或跨 technical "
                "融合樹 0/7。舊版 7-dataset / 6-cell-line candidate-tree census 使用 "
                "50 kb/legacy "
                "grouping 且含 mixed PS，只能作歷史參考，不能與新版 W 合併或相減。"
            ),
        },
        {
            "id": "topology_completion_table_block",
            "type": "table",
            "tableId": "topology_completion_table",
            "layout": "full",
        },
        {
            "id": "definition_heading",
            "type": "markdown",
            "body": (
                "## 區域怎麼定義\n\n"
                "`W = maximal connected component(G_H)`，其中 "
                "`H = dataset × chromosome × exact PS × primary HP family`；"
                "邊只由 fixed R/A endpoint 的 distinct-molecule support 決定。"
                "座標距離不進入 edge 或 component 的目標函數。"
            ),
        },
        {"id": "method_flow", "type": "html", "body": method_svg(total_w, total_singletons)},
        {
            "id": "partial_read_rule",
            "type": "markdown",
            "body": (
                "## Partial reads 與 read chain\n\n"
                "對排序位點 `s1<s2<s3`，pattern `R-X-A` 只替 `(s1,s3)` 投一票 `RA`；"
                "不得補出 `(s1,s2)` 或 `(s2,s3)`。若 `(s1,s2)` 與 `(s2,s3)` "
                "由各自的 molecule 集合達到 threshold，三點仍可形成同一 W；"
                "這不表示有一條 read 或同一細胞橫跨整個 W。`RA/AR` 依 genomic site "
                "order 編碼，不代表祖先—後代方向。"
            ),
        },
        {"id": "definition_table_block", "type": "table", "tableId": "definition_table", "layout": "full"},
        {
            "id": "molecule_row_heading",
            "type": "markdown",
            "body": (
                "## Molecule-row 容器過濾\n\n"
                "每個 canonical molecule row 只能落入 primary exact-PS×HP、"
                "missing PS 排除或 nonprimary HP 排除之一；三類在各資料集內守恆。"
            ),
        },
        {
            "id": "molecule_row_chart_block",
            "type": "chart",
            "chartId": "molecule_row_outcome_chart",
        },
        {
            "id": "locus_heading",
            "type": "markdown",
            "body": (
                "## Locus、membership 與 component 分母分開\n\n"
                "Physical loci、PS×HP memberships 與 graph components 是三種不同統計單位。"
                "因此不畫 `S loci → W regions` 的假漏斗；每條圖只比較同一單位。"
            ),
        },
        {"id": "locus_chart_block", "type": "chart", "chartId": "locus_stage_chart"},
        {"id": "component_chart_block", "type": "chart", "chartId": "component_outcome_chart"},
        {"id": "dataset_funnel_table_block", "type": "table", "tableId": "dataset_funnel_audit", "layout": "full"},
        {
            "id": "k_heading",
            "type": "markdown",
            "body": (
                "## W 的 k 與座標尺度\n\n"
                "k 反映每個 read-linked region 含多少 sSNV memberships。"
                "Span 只描述從最左到最右位點的 envelope；長 span 常由多段 read-supported chain 串接。"
            ),
        },
        {"id": "k_chart_block", "type": "chart", "chartId": "k_distribution_chart"},
        {
            "id": "exact_k_table_block",
            "type": "table",
            "tableId": "exact_k_audit",
            "layout": "full",
        },
        {"id": "span_chart_block", "type": "chart", "chartId": "span_distribution_chart"},
        {
            "id": "distance_heading",
            "type": "markdown",
            "body": (
                "## 為何 50 kb 降為 QC\n\n"
                "同一 read 對兩端的 fixed calls 是直接觀測，read alignment 本身已自然限制可連距離；"
                "額外的全域 50 kb 門檻會丟棄通過目前 filters 的觀測支持。另一方面，"
                "這不代表所有長距離 edge 都是真實生物連結；它們仍可能受 repeat、"
                "misalignment、CNV/LOH 或低支持影響，所以必須保留分層 QC。"
                "本報告也實際移除所有 >50 kb direct edges 後重建 components，"
                "直接量化硬切造成的 partition、W 與 membership 改變。"
            ),
        },
        {"id": "distance_w_chart_block", "type": "chart", "chartId": "distance_w_qc_chart"},
        {"id": "direct_edge_chart_block", "type": "chart", "chartId": "direct_long_edge_chart"},
        {
            "id": "distance_cut_impact_chart_block",
            "type": "chart",
            "chartId": "distance_cut_impact_chart",
        },
        {"id": "long_edge_support_chart_block", "type": "chart", "chartId": "long_edge_support_chart"},
        {"id": "long_edge_state_chart_block", "type": "chart", "chartId": "long_edge_state_chart"},
        {"id": "dataset_distance_table_block", "type": "table", "tableId": "dataset_distance_audit", "layout": "full"},
        {
            "id": "threshold_heading",
            "type": "markdown",
            "body": (
                "## Read-support threshold 敏感度\n\n"
                "Primary threshold=3 是可稽核的 operational setting，不是已由 held-out error model "
                "校準的生物真值。Threshold=3 由 component、membership、edge rows 獨立重算；"
                "1/2/5 則由目前 22 份 chromosome receipts 重加總並與 dataset summary 對帳。"
                "提高 threshold 時 retained "
                "edges 必然不增，但一個 W 可能裂成兩個仍為 k≥2 的 W，所以 W 數不必然單調下降。"
            ),
        },
        {"id": "threshold_chart_block", "type": "chart", "chartId": "threshold_sensitivity_chart"},
        {"id": "threshold_table_block", "type": "table", "tableId": "threshold_audit"},
        {
            "id": "chromosome_heading",
            "type": "markdown",
            "body": (
                "## 全部 154 個 dataset × chromosome 明細\n\n"
                "下表完整列出 chr1–22；每個 dataset 的 total 均可由表中 22 列重算。"
            ),
        },
        {"id": "chromosome_table_block", "type": "table", "tableId": "chromosome_audit", "layout": "full"},
        {
            "id": "qa_heading",
            "type": "markdown",
            "body": (
                "## 驗證、限制與可主張範圍\n\n"
                "可主張：exact-PS/HP-local、multi-molecule-supported、fixed-endpoint read-linkage region；"
                "distance 不再決定 region。\n\n"
                "不可主張：唯一 topology、cellular clone 數、subclone truth、祖先—後代方向、"
                "跨 HP 同細胞配對、跨 PS 的 HP 對應、每個 W 都存在一條 full-span molecule，"
                "或已由 VAF/CCF/CN/甲基證實。**DORADO 與其餘 6 組目前都沒有可用的 "
                "production strict directed-topology PASS receipt。** RR-only edge 也不等於 "
                "somatic ALT 共現；endpoint edge 是無向 linkage，不是 evolutionary "
                "parent→child。"
            ),
        },
        {"id": "qa_table_block", "type": "table", "tableId": "qa_table"},
        {
            "id": "source_heading",
            "type": "markdown",
            "body": (
                "## 資料來源與 freshness\n\n"
                "- Quantitative source：all-7 strict-region independently recomputed package；"
                "所有 chromosome receipt、TSV identity 與守恆檢查均 fail-closed。\n"
                "- Completion audit：在 declared observation/research search roots 盤點 "
                "v4 strict run receipts；找到的最新 HCC1395 v4 run 只到 partition，"
                "sample topology identity 為 null，未找到任何 7×22 production topology "
                "receipt。HCC1395 早期 technical topology receipt 明示不具 cohort-final "
                "資格；legacy all-7 census 的 claim scope 也只到 candidate mutation-state "
                "trees。\n"
                "- Builder provenance：HCC1395、HCC1395_DORADO、COLO829、H1437、H2009、"
                "HCC1937 的 132 份 chromosome artifacts 使用已封存的 pre-hotfix builder "
                "`7260a763…`；HCC1954 的 22 份使用 post-hotfix builder `912721f9…`。"
                "前 6 組 production container TSV 共 156,316 列均未命中 missing-PS aliases，"
                "且 HCC1937 chr21 以新版重建後 component、membership、edge、container "
                "解壓內容 byte-identical；因此本報告採 data-specific no-trigger equivalence，"
                "但不宣稱 154 份 artifacts 為同一 builder 版本。\n"
                "- Contextual KB：`kb-05-tools-longphase-s`（runtime-fact, verified, "
                "last_verified 2026-07-11）與 `kb-03-file-formats-vcf-longphase` "
                "（reference-summary, verified, last_verified 2026-07-11）。"
                "LongPhase-S block N50 僅作背景，不能替代本次 exact-PS 實測。\n"
                "- 圖表由 canonical artifact renderer 產生並抽取為 inline SVG；表格保留語意 fallback。"
            ),
        },
    ]

    manifest = {
        "version": 1,
        "surface": "report",
        "title": TITLE,
        "description": (
            "Exact PS×HP strict fixed-endpoint read-linkage regions across seven "
            "technical datasets and autosomes, with 50 kb retained as QC only."
        ),
        "generatedAt": generated_at,
        "sources": sources,
        "cards": [
            {
                "id": "scope_card",
                "description": "Full Task Type B scope and independence boundary.",
                "dataset": "headline",
                "sourceId": source_id("headline"),
                "metrics": [
                    {"label": "Technical datasets", "field": "technical_datasets", "format": "number"},
                    {"label": "Biological cell lines", "field": "biological_cell_lines", "format": "number"},
                    {"label": "Dataset×chromosome rows", "field": "autosome_dataset_rows", "format": "number"},
                ],
            },
            {
                "id": "locus_card",
                "description": (
                    "Per-dataset unique-locus records summed across seven technical "
                    "datasets; not cross-dataset globally deduplicated loci."
                ),
                "dataset": "headline",
                "sourceId": source_id("headline"),
                "metrics": [
                    {
                        "label": "Dataset-summed S records",
                        "field": "dataset_summed_candidate_locus_records",
                        "format": "number",
                    },
                    {
                        "label": "Active records",
                        "field": "dataset_summed_active_locus_records",
                        "format": "number",
                    },
                    {
                        "label": "Records in ≥1 W",
                        "field": "dataset_summed_locus_records_in_W",
                        "format": "number",
                    },
                ],
            },
            {
                "id": "component_card",
                "description": "Conserved graph-component partition.",
                "dataset": "headline",
                "sourceId": source_id("headline"),
                "metrics": [
                    {"label": "All components", "field": "all_components", "format": "number"},
                    {"label": "k=1 abstain", "field": "singleton_k1_components", "format": "number"},
                    {"label": "Tree-eligible W", "field": "tree_eligible_W", "format": "number"},
                ],
            },
            {
                "id": "edge_card",
                "description": "Retained direct endpoint-edge evidence.",
                "dataset": "headline",
                "sourceId": source_id("headline"),
                "metrics": [
                    {"label": "Retained direct edges", "field": "retained_direct_edges", "format": "number"},
                    {"label": "Direct edges>50 kb", "field": "direct_edges_gt_50kb", "format": "number"},
                ],
            },
            {
                "id": "distance_card",
                "description": "Observed 50 kb region diagnostics; both metrics use W as denominator.",
                "dataset": "headline",
                "sourceId": source_id("headline"),
                "metrics": [
                    {"label": "W span>50 kb", "field": "W_span_gt_50kb", "format": "number"},
                    {"label": "W max-gap>50 kb", "field": "W_adjacent_gap_gt_50kb", "format": "number"},
                ],
            },
            {
                "id": "counterfactual_card",
                "description": "Counterfactual result after removing only direct edges longer than 50 kb.",
                "dataset": "headline",
                "sourceId": source_id("headline"),
                "metrics": [
                    {
                        "label": "W changed if cut",
                        "field": "W_partition_changed_if_50kb_cut",
                        "format": "number",
                    },
                    {
                        "label": "W after cut",
                        "field": "counterfactual_W_after_50kb_cut",
                        "format": "number",
                    },
                ],
            },
        ],
        "charts": charts,
        "tables": tables,
        "blocks": blocks,
    }
    snapshot = {
        "version": 1,
        "generatedAt": generated_at,
        "status": "ready",
        "datasets": queried_snapshot,
    }
    return {
        "surface": "report",
        "manifest": manifest,
        "snapshot": snapshot,
        "sources": sources,
        "package_info": {
            "root": ".",
            "manifestPath": output_path.name,
            "snapshotPath": output_path.name,
        },
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    data_path = args.data.resolve(strict=True)
    output_path = args.output.resolve()
    if output_path.exists():
        raise ValueError(f"output must not exist: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path = output_path.parent / ARTIFACT_RECEIPT_NAME
    sidecar_path = output_path.parent / ARTIFACT_RECEIPT_SIDECAR_NAME
    if receipt_path.exists() or sidecar_path.exists():
        raise ValueError(
            "artifact receipt outputs must not exist: "
            f"{receipt_path}, {sidecar_path}"
        )
    stage_id = uuid4().hex
    candidate_path = output_path.parent / f".{output_path.name}.pending-{stage_id}"
    candidate_receipt_path = (
        output_path.parent / f".{ARTIFACT_RECEIPT_NAME}.pending-{stage_id}"
    )
    candidate_sidecar_path = (
        output_path.parent / f".{ARTIFACT_RECEIPT_SIDECAR_NAME}.pending-{stage_id}"
    )
    artifact: dict[str, Any] | None = None
    release_checks: dict[str, bool] | None = None
    try:
        artifact = build_artifact(read_json(data_path), data_path, output_path)
        release_checks = artifact_release_checks(artifact)
        with candidate_path.open("x", encoding="utf-8") as handle:
            json.dump(
                artifact,
                handle,
                ensure_ascii=False,
                sort_keys=True,
                indent=2,
                allow_nan=False,
            )
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(candidate_path, output_path)

        query_paths = sorted((output_path.parent / "queries").glob("*.sql"))
        expected_query_names = {
            f"{dataset_id}.sql"
            for dataset_id in artifact["snapshot"]["datasets"]
        }
        if (
            not query_paths
            or {path.name for path in query_paths} != expected_query_names
        ):
            raise ValueError(
                "published SQLite query files do not exactly match snapshot datasets"
            )
        required_outputs = [
            output_path,
            output_path.parent / "data" / "report.sqlite",
            *query_paths,
        ]
        if any(not path.is_file() for path in required_outputs):
            raise ValueError("one or more required artifact outputs are absent")
        receipt = {
            "schema_name": "intersubmod.strict_region_all7_artifact_receipt",
            "schema_version": "1.0.0",
            "created_at_utc": utc_now(),
            "all_pass": True,
            "scope": {
                "task_type": "B_comprehensive_validation",
                "datasets": list(CANONICAL_DATASETS),
                "chromosomes": list(AUTOSOMES),
            },
            "checks": release_checks,
            "inputs": {
                "data": file_identity(data_path),
                "data_receipt": file_identity(data_path.parent / "receipt.json"),
                "data_receipt_sidecar": file_identity(
                    data_path.parent / "receipt.json.sha256"
                ),
            },
            "outputs": {
                str(path.relative_to(output_path.parent)): file_identity(path)
                for path in required_outputs
            },
        }
        write_json_exclusive(candidate_receipt_path, receipt)
        write_text_exclusive(
            candidate_sidecar_path,
            (
                f"{sha256_path(candidate_receipt_path)}  "
                f"{ARTIFACT_RECEIPT_NAME}\n"
            ),
        )
        os.replace(candidate_receipt_path, receipt_path)
        os.replace(candidate_sidecar_path, sidecar_path)
    except Exception as error:
        candidates = (
            (candidate_path, candidate_path.name),
            (candidate_receipt_path, candidate_receipt_path.name),
            (candidate_sidecar_path, candidate_sidecar_path.name),
            (receipt_path, receipt_path.name),
            (sidecar_path, sidecar_path.name),
            (output_path, output_path.name),
            (output_path.parent / "data" / "report.sqlite", "report.sqlite"),
            (output_path.parent / "queries", "queries"),
        )
        if any(path.exists() for path, _ in candidates):
            archive_dir = (
                output_path.parent
                / "archive"
                / "failed_artifact_publish"
                / stage_id
            )
            archive_dir.mkdir(parents=True, exist_ok=False)
            for path, archive_name in candidates:
                if path.exists():
                    os.replace(path, archive_dir / archive_name)
            raise ValueError(
                f"artifact publication failed; preserved evidence at {archive_dir}"
            ) from error
        raise
    assert artifact is not None
    print(
        json.dumps(
            {
                "output": str(output_path),
                "charts": len(artifact["manifest"]["charts"]),
                "tables": len(artifact["manifest"]["tables"]),
                "blocks": len(artifact["manifest"]["blocks"]),
                "snapshot_datasets": len(artifact["snapshot"]["datasets"]),
                "receipt": str(receipt_path),
                "release_checks": release_checks,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
