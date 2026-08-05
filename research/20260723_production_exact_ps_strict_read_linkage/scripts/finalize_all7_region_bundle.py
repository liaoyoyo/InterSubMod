#!/usr/bin/env python3
"""Fail closed and publish READY.json last for the complete all-7 report bundle."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence
from uuid import uuid4


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
DATA_OUTPUT_NAMES = frozenset(
    {
        "all7_report_data.json",
        "dataset_summary.tsv",
        "chromosome_summary.tsv",
        "k_band_distribution.tsv",
        "span_band_distribution.tsv",
        "exact_k_distribution.tsv",
        "threshold_sensitivity.tsv",
        "topology_completion_status.tsv",
    }
)
EXPECTED_QA_EVIDENCE = frozenset(
    {
        "desktop_overview.png",
        "mobile_overview.png",
        "desktop_method_flow.png",
        "desktop_distance_qc.png",
        "desktop_distance_cut_impact.png",
        "nojs_desktop_overview.png",
        "nojs_mobile_overview.png",
    }
)
EXPECTED_QA_MODES = frozenset(
    {
        ("javascript_enabled", "desktop", 1440, 1100),
        ("javascript_enabled", "mobile", 390, 844),
        ("javascript_disabled", "desktop", 1440, 1100),
        ("javascript_disabled", "mobile", 390, 844),
    }
)


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
    path = path.resolve(strict=True)
    if not path.is_file():
        raise ValueError(f"expected a regular file: {path}")
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def verify_identity(
    spec: Mapping[str, Any],
    path: Path,
    *,
    label: str,
) -> None:
    observed = file_identity(path)
    declared_size = (
        spec.get("size_bytes")
        if "size_bytes" in spec
        else spec.get("bytes")
    )
    if (
        spec.get("path") != observed["path"]
        or declared_size != observed["size_bytes"]
        or spec.get("sha256") != observed["sha256"]
    ):
        raise ValueError(f"{label}: file identity mismatch for {path}")


def verify_sidecar(sidecar_path: Path, target_path: Path) -> None:
    expected = f"{sha256_path(target_path)}  {target_path.name}"
    if sidecar_path.read_text(encoding="ascii").strip() != expected:
        raise ValueError(f"{sidecar_path}: checksum sidecar mismatch")


def require_all_true(value: Any, *, label: str) -> Mapping[str, Any]:
    if (
        not isinstance(value, Mapping)
        or not value
        or any(item is not True for item in value.values())
    ):
        raise ValueError(f"{label}: expected a nonempty all-true check map")
    return value


def require_scope(scope: Any, *, label: str) -> None:
    if (
        not isinstance(scope, Mapping)
        or scope.get("task_type") != "B_comprehensive_validation"
        or scope.get("datasets") != list(CANONICAL_DATASETS)
        or scope.get("chromosomes") != list(AUTOSOMES)
    ):
        raise ValueError(f"{label}: canonical full-scope contract mismatch")


def require_extraction_receipt_grid(data: Mapping[str, Any]) -> None:
    checks = data.get("checks")
    if (
        not isinstance(checks, Mapping)
        or checks.get("all_154_extraction_receipt_identities_verified") is not True
    ):
        raise ValueError("report data lacks the verified 7×22 extraction receipt grid")
    inputs = data.get("inputs")
    if not isinstance(inputs, Mapping) or set(inputs) != set(CANONICAL_DATASETS):
        raise ValueError("report data input provenance is incomplete")
    count = 0
    for dataset in CANONICAL_DATASETS:
        provenance = inputs.get(dataset)
        receipts = (
            provenance.get("extraction_receipts")
            if isinstance(provenance, Mapping)
            else None
        )
        if not isinstance(receipts, Mapping) or set(receipts) != set(AUTOSOMES):
            raise ValueError(
                f"{dataset}: extraction receipt provenance is not an exact chr1–22 grid"
            )
        for chrom, spec in receipts.items():
            if (
                not isinstance(spec, Mapping)
                or not isinstance(spec.get("path"), str)
                or not Path(spec["path"]).is_absolute()
                or not isinstance(spec.get("size_bytes"), int)
                or isinstance(spec.get("size_bytes"), bool)
                or spec["size_bytes"] < 0
                or not isinstance(spec.get("sha256"), str)
                or len(spec["sha256"]) != 64
                or any(char not in "0123456789abcdef" for char in spec["sha256"])
            ):
                raise ValueError(
                    f"{dataset}/{chrom}: extraction receipt identity is invalid"
                )
            count += 1
    if count != len(CANONICAL_DATASETS) * len(AUTOSOMES):
        raise ValueError("report data extraction receipt grid is not exactly 7×22")


def artifact_visual_ids(
    artifact: Mapping[str, Any],
) -> tuple[list[str], list[str]]:
    manifest = artifact.get("manifest")
    if not isinstance(manifest, Mapping):
        raise ValueError("artifact manifest is absent")
    charts = manifest.get("charts")
    tables = manifest.get("tables")
    blocks = manifest.get("blocks")
    if not all(isinstance(value, list) for value in (charts, tables, blocks)):
        raise ValueError("artifact chart/table/block arrays are absent")

    def ids(rows: Sequence[Any], key: str, label: str) -> list[str]:
        values = [
            row.get(key) if isinstance(row, Mapping) else None for row in rows
        ]
        if (
            any(not isinstance(value, str) or not value for value in values)
            or len(values) != len(set(values))
        ):
            raise ValueError(f"artifact {label} IDs are missing or duplicated")
        return sorted(values)

    chart_ids = ids(charts, "id", "chart")
    table_ids = ids(tables, "id", "table")
    visible_chart_ids = ids(
        [
            block
            for block in blocks
            if isinstance(block, Mapping) and block.get("type") == "chart"
        ],
        "chartId",
        "visible chart",
    )
    visible_table_ids = ids(
        [
            block
            for block in blocks
            if isinstance(block, Mapping) and block.get("type") == "table"
        ],
        "tableId",
        "visible table",
    )
    if chart_ids != visible_chart_ids or table_ids != visible_table_ids:
        raise ValueError("artifact visible chart/table IDs do not match definitions")
    return chart_ids, table_ids


def validate_data_package(
    *,
    data_path: Path,
    receipt_path: Path,
    sidecar_path: Path,
) -> dict[str, Any]:
    data = read_json(data_path)
    if (
        data.get("schema_name") != "intersubmod.strict_region_all7_report_data"
        or data.get("schema_version") != "1.1.0"
        or data.get("all_pass") is not True
    ):
        raise ValueError("report data schema/version/PASS contract mismatch")
    require_scope(data.get("scope"), label="report data")
    data_checks = require_all_true(data.get("checks"), label="report data")
    require_extraction_receipt_grid(data)

    dataset_rows = data.get("datasets")
    chromosome_rows = data.get("chromosomes")
    if not isinstance(dataset_rows, list) or not isinstance(chromosome_rows, list):
        raise ValueError("report data scope rows are absent")
    if (
        [row.get("dataset") for row in dataset_rows] != list(CANONICAL_DATASETS)
        or any(row.get("all_pass") is not True for row in dataset_rows)
    ):
        raise ValueError("report data does not contain the exact seven PASS datasets")
    observed_pairs = [
        (row.get("dataset"), row.get("chrom")) for row in chromosome_rows
    ]
    expected_pairs = {
        (dataset, chrom)
        for dataset in CANONICAL_DATASETS
        for chrom in AUTOSOMES
    }
    if (
        len(observed_pairs) != len(expected_pairs)
        or set(observed_pairs) != expected_pairs
        or any(row.get("all_pass") is not True for row in chromosome_rows)
    ):
        raise ValueError("report data is not the exact 7×22 PASS grid")

    verify_sidecar(sidecar_path, receipt_path)
    receipt = read_json(receipt_path)
    if (
        receipt.get("schema_name")
        != "intersubmod.strict_region_all7_report_data_receipt"
        or receipt.get("schema_version") != "1.1.0"
        or receipt.get("all_pass") is not True
        or receipt.get("checks") != data_checks
    ):
        raise ValueError("report-data receipt contract mismatch")
    outputs = receipt.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != DATA_OUTPUT_NAMES:
        raise ValueError("report-data receipt output set is incomplete or unexpected")
    for name, spec in outputs.items():
        if not isinstance(spec, Mapping):
            raise ValueError(f"report-data output identity is invalid: {name}")
        verify_identity(
            spec,
            data_path.parent / name,
            label=f"report-data output {name}",
        )
    return data


def validate_artifact_package(
    *,
    artifact_path: Path,
    receipt_path: Path,
    sidecar_path: Path,
    data_path: Path,
    data_receipt_path: Path,
    data_sidecar_path: Path,
) -> tuple[dict[str, Any], list[str], list[str]]:
    artifact = read_json(artifact_path)
    manifest = artifact.get("manifest")
    snapshot = artifact.get("snapshot")
    if (
        artifact.get("surface") != "report"
        or not isinstance(manifest, Mapping)
        or manifest.get("surface") != "report"
        or not isinstance(snapshot, Mapping)
        or snapshot.get("status") != "ready"
    ):
        raise ValueError("artifact ready-report contract mismatch")
    snapshot_datasets = snapshot.get("datasets")
    if not isinstance(snapshot_datasets, Mapping):
        raise ValueError("artifact snapshot datasets are absent")
    dataset_rows = snapshot_datasets.get("dataset_summary")
    chromosome_rows = snapshot_datasets.get("chromosome_detail")
    if not isinstance(dataset_rows, list) or not isinstance(chromosome_rows, list):
        raise ValueError("artifact full-scope datasets are absent")
    artifact_dataset_names = [row.get("dataset") for row in dataset_rows]
    if (
        len(artifact_dataset_names) != len(CANONICAL_DATASETS)
        or set(artifact_dataset_names) != set(CANONICAL_DATASETS)
    ):
        raise ValueError("artifact dataset_summary is not the exact canonical seven")
    observed_pairs = [
        (row.get("dataset"), row.get("chrom")) for row in chromosome_rows
    ]
    expected_pairs = {
        (dataset, chrom)
        for dataset in CANONICAL_DATASETS
        for chrom in AUTOSOMES
    }
    if len(observed_pairs) != len(expected_pairs) or set(observed_pairs) != expected_pairs:
        raise ValueError("artifact chromosome_detail is not the exact 7×22 grid")
    chart_ids, table_ids = artifact_visual_ids(artifact)

    verify_sidecar(sidecar_path, receipt_path)
    receipt = read_json(receipt_path)
    if (
        receipt.get("schema_name")
        != "intersubmod.strict_region_all7_artifact_receipt"
        or receipt.get("schema_version") != "1.0.0"
        or receipt.get("all_pass") is not True
    ):
        raise ValueError("artifact receipt schema/version/PASS contract mismatch")
    require_scope(receipt.get("scope"), label="artifact receipt")
    require_all_true(receipt.get("checks"), label="artifact receipt")
    inputs = receipt.get("inputs")
    if not isinstance(inputs, Mapping):
        raise ValueError("artifact receipt inputs are absent")
    for key, path in {
        "data": data_path,
        "data_receipt": data_receipt_path,
        "data_receipt_sidecar": data_sidecar_path,
    }.items():
        spec = inputs.get(key)
        if not isinstance(spec, Mapping):
            raise ValueError(f"artifact receipt input identity is absent: {key}")
        verify_identity(spec, path, label=f"artifact receipt input {key}")

    expected_query_paths = {
        f"queries/{dataset_id}.sql" for dataset_id in snapshot_datasets
    }
    expected_outputs = {
        artifact_path.name,
        "data/report.sqlite",
        *expected_query_paths,
    }
    outputs = receipt.get("outputs")
    if not isinstance(outputs, Mapping) or set(outputs) != expected_outputs:
        raise ValueError("artifact receipt output set is incomplete or unexpected")
    for relative_path, spec in outputs.items():
        if not isinstance(spec, Mapping):
            raise ValueError(f"artifact output identity is invalid: {relative_path}")
        verify_identity(
            spec,
            artifact_path.parent / relative_path,
            label=f"artifact output {relative_path}",
        )
    return artifact, chart_ids, table_ids


def validate_delivery(
    *,
    receipt_path: Path,
    artifact_path: Path,
    html_path: Path,
    chart_ids: Sequence[str],
) -> dict[str, Any]:
    receipt = read_json(receipt_path)
    if receipt.get("version") != 1 or receipt.get("ok") is not True:
        raise ValueError("delivery receipt version/PASS contract mismatch")
    stages = receipt.get("stages")
    required_stages = {
        "canonicalValidation",
        "package",
        "staticSvgExtraction",
        "officialBrowserVerification",
        "atomicOutputRename",
    }
    if (
        not isinstance(stages, Mapping)
        or set(stages) != required_stages
        or any(value != "passed" for value in stages.values())
    ):
        raise ValueError("delivery receipt stage map is incomplete or failed")
    for key, path in {"input": artifact_path, "output": html_path}.items():
        spec = receipt.get(key)
        if not isinstance(spec, Mapping):
            raise ValueError(f"delivery receipt identity is absent: {key}")
        verify_identity(spec, path, label=f"delivery {key}")
    browser = receipt.get("browser")
    if not isinstance(browser, Mapping):
        raise ValueError("delivery browser identity is absent")
    browser_path = browser.get("path")
    if not isinstance(browser_path, str):
        raise ValueError("delivery browser path is absent")
    verify_identity(browser, Path(browser_path), label="delivery browser")
    static_svg = receipt.get("staticSvg")
    if (
        not isinstance(static_svg, Mapping)
        or static_svg.get("chartIds") != sorted(chart_ids)
        or static_svg.get("count") != len(chart_ids)
        or static_svg.get("themeVariants") != len(chart_ids) * 2
        or static_svg.get("embeddedThemeVariants") != len(chart_ids) * 2
    ):
        raise ValueError("delivery static-SVG ID/count contract mismatch")
    verification = receipt.get("verification")
    if (
        not isinstance(verification, Mapping)
        or verification.get("viewports") != [1440, 390]
        or verification.get("sourceDialog") != "passed"
        or verification.get("counts", {}).get("charts") != len(chart_ids)
    ):
        raise ValueError("delivery official browser verification contract mismatch")
    return receipt


def validate_visual_qa(
    *,
    receipt_path: Path,
    artifact_path: Path,
    html_path: Path,
    chart_ids: Sequence[str],
    table_ids: Sequence[str],
) -> dict[str, Any]:
    receipt = read_json(receipt_path)
    if (
        receipt.get("schema_name") != "intersubmod.strict_region_html_visual_qa"
        or receipt.get("schema_version") != "2.0.0"
        or receipt.get("all_pass") is not True
    ):
        raise ValueError("visual-QA receipt schema/version/PASS contract mismatch")
    inputs = receipt.get("inputs")
    if not isinstance(inputs, Mapping):
        raise ValueError("visual-QA inputs are absent")
    for key, path in {"artifact": artifact_path, "html": html_path}.items():
        spec = inputs.get(key)
        if not isinstance(spec, Mapping):
            raise ValueError(f"visual-QA input identity is absent: {key}")
        verify_identity(spec, path, label=f"visual-QA input {key}")
    expected = receipt.get("expected")
    if (
        not isinstance(expected, Mapping)
        or expected.get("chart_ids") != sorted(chart_ids)
        or expected.get("table_ids") != sorted(table_ids)
        or expected.get("visible_charts") != len(chart_ids)
        or expected.get("visible_tables") != len(table_ids)
        or expected.get("static_svg_variants") != len(chart_ids) * 2
        or expected.get("modes")
        != ["javascript_enabled", "javascript_disabled"]
    ):
        raise ValueError("visual-QA expected ID/count contract mismatch")
    browser = receipt.get("browser")
    executable = browser.get("executable") if isinstance(browser, Mapping) else None
    if not isinstance(executable, Mapping) or not isinstance(
        executable.get("path"), str
    ):
        raise ValueError("visual-QA browser executable identity is absent")
    verify_identity(
        executable,
        Path(executable["path"]),
        label="visual-QA browser executable",
    )
    viewports = receipt.get("viewports")
    if not isinstance(viewports, list):
        raise ValueError("visual-QA viewport results are absent")
    observed_modes = set()
    for result in viewports:
        if not isinstance(result, Mapping):
            raise ValueError("visual-QA viewport result is invalid")
        viewport = result.get("viewport")
        checks = result.get("checks")
        metrics = result.get("metrics")
        if (
            not isinstance(viewport, Mapping)
            or not isinstance(checks, Mapping)
            or not isinstance(metrics, Mapping)
            or result.get("all_pass") is not True
            or not checks
            or any(value is not True for value in checks.values())
            or metrics.get("horizontalOverflowPx", 2) > 1
        ):
            raise ValueError("visual-QA viewport result contains a failed check")
        observed_modes.add(
            (
                result.get("mode"),
                viewport.get("label"),
                viewport.get("width"),
                viewport.get("height"),
            )
        )
    if observed_modes != EXPECTED_QA_MODES or len(viewports) != len(EXPECTED_QA_MODES):
        raise ValueError("visual-QA JS mode × viewport grid is incomplete")

    evidence = receipt.get("evidence")
    evidence_files = evidence.get("files") if isinstance(evidence, Mapping) else None
    if (
        not isinstance(evidence, Mapping)
        or evidence.get("all_expected_files_present") is not True
        or not isinstance(evidence_files, Mapping)
        or set(evidence_files) != EXPECTED_QA_EVIDENCE
    ):
        raise ValueError("visual-QA evidence-file set is incomplete")
    for name, spec in evidence_files.items():
        if not isinstance(spec, Mapping):
            raise ValueError(f"visual-QA evidence identity is invalid: {name}")
        verify_identity(
            spec,
            receipt_path.parent / name,
            label=f"visual-QA evidence {name}",
        )
    return receipt


def write_ready_last(
    *,
    output_path: Path,
    value: Mapping[str, Any],
) -> Path:
    if output_path.name != "READY.json":
        raise ValueError("--output must be named READY.json")
    sidecar_path = output_path.with_name("READY.json.sha256")
    if output_path.exists() or sidecar_path.exists():
        raise ValueError(
            f"READY outputs already exist; refusing to overwrite: "
            f"{output_path}, {sidecar_path}"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    stage_id = uuid4().hex
    pending_ready = output_path.parent / f".READY.json.pending-{stage_id}"
    pending_sidecar = output_path.parent / f".READY.json.sha256.pending-{stage_id}"
    published_sidecar = False
    published_ready = False
    try:
        with pending_ready.open("x", encoding="utf-8") as handle:
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
        with pending_sidecar.open("x", encoding="ascii") as handle:
            handle.write(f"{sha256_path(pending_ready)}  READY.json\n")
            handle.flush()
            os.fsync(handle.fileno())

        # Publish the sidecar first. READY.json is renamed last and is the only
        # consumer commit marker.
        os.replace(pending_sidecar, sidecar_path)
        published_sidecar = True
        os.replace(pending_ready, output_path)
        published_ready = True
        directory_fd = os.open(output_path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except Exception:
        candidates = [
            path
            for path in (
                pending_ready,
                pending_sidecar,
                output_path if published_ready else None,
                sidecar_path if published_sidecar else None,
            )
            if path is not None and path.exists()
        ]
        if candidates:
            archive = (
                output_path.parent
                / "archive"
                / "failed_finalization"
                / stage_id
            )
            archive.mkdir(parents=True, exist_ok=False)
            for path in candidates:
                os.replace(path, archive / path.name)
        raise
    return sidecar_path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--data-receipt", required=True, type=Path)
    parser.add_argument("--data-sidecar", required=True, type=Path)
    parser.add_argument("--artifact", required=True, type=Path)
    parser.add_argument("--artifact-receipt", required=True, type=Path)
    parser.add_argument("--artifact-sidecar", required=True, type=Path)
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--delivery-receipt", required=True, type=Path)
    parser.add_argument("--visual-qa-receipt", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    paths = {
        key: value.resolve(strict=True)
        for key, value in {
            "data": args.data,
            "data_receipt": args.data_receipt,
            "data_sidecar": args.data_sidecar,
            "artifact": args.artifact,
            "artifact_receipt": args.artifact_receipt,
            "artifact_sidecar": args.artifact_sidecar,
            "html": args.html,
            "delivery_receipt": args.delivery_receipt,
            "visual_qa_receipt": args.visual_qa_receipt,
        }.items()
    }
    output_path = args.output.resolve()
    bundle_root = paths["artifact"].parent
    expected_layout = {
        "data": bundle_root / "data" / "all7_report_data.json",
        "data_receipt": bundle_root / "data" / "receipt.json",
        "data_sidecar": bundle_root / "data" / "receipt.json.sha256",
        "artifact_receipt": bundle_root / "artifact_receipt.json",
        "artifact_sidecar": bundle_root / "artifact_receipt.json.sha256",
        "delivery_receipt": bundle_root / "delivery_receipt.json",
        "visual_qa_receipt": bundle_root / "visual_qa" / "visual_qa_receipt.json",
    }
    if output_path.parent != bundle_root:
        raise ValueError("READY.json must be published in the artifact bundle root")
    if paths["html"].parent != bundle_root:
        raise ValueError("HTML must be published in the artifact bundle root")
    for key, expected in expected_layout.items():
        if paths[key] != expected.resolve():
            raise ValueError(
                f"noncanonical bundle layout for {key}: "
                f"{paths[key]} != {expected}"
            )

    data = validate_data_package(
        data_path=paths["data"],
        receipt_path=paths["data_receipt"],
        sidecar_path=paths["data_sidecar"],
    )
    artifact, chart_ids, table_ids = validate_artifact_package(
        artifact_path=paths["artifact"],
        receipt_path=paths["artifact_receipt"],
        sidecar_path=paths["artifact_sidecar"],
        data_path=paths["data"],
        data_receipt_path=paths["data_receipt"],
        data_sidecar_path=paths["data_sidecar"],
    )
    validate_delivery(
        receipt_path=paths["delivery_receipt"],
        artifact_path=paths["artifact"],
        html_path=paths["html"],
        chart_ids=chart_ids,
    )
    validate_visual_qa(
        receipt_path=paths["visual_qa_receipt"],
        artifact_path=paths["artifact"],
        html_path=paths["html"],
        chart_ids=chart_ids,
        table_ids=table_ids,
    )

    checks = {
        "data_receipt_pass_and_identities_valid": True,
        "artifact_receipt_pass_and_identities_valid": True,
        "delivery_receipt_pass_and_identities_valid": True,
        "visual_qa_receipt_pass_and_identities_valid": True,
        "scope_exact_7_datasets_x_22_autosomes": True,
        "live_and_static_visual_ids_exact": True,
        "javascript_enabled_and_disabled_1440_390_qa_pass": True,
        "ready_is_last_consumer_commit_marker": True,
    }
    ready = {
        "schema_name": "intersubmod.strict_region_all7_bundle_ready",
        "schema_version": "1.0.0",
        "created_at_utc": utc_now(),
        "all_pass": True,
        "scope": {
            "task_type": "B_comprehensive_validation",
            "datasets": list(CANONICAL_DATASETS),
            "chromosomes": list(AUTOSOMES),
            "dataset_chromosome_rows": len(data["chromosomes"]),
        },
        "checks": checks,
        "commit_marker": {
            "file": "READY.json",
            "sidecar": "READY.json.sha256",
            "semantics": (
                "Consumers may use the bundle only when READY.json exists, its "
                "sidecar matches, and all identities in this document still match."
            ),
        },
        "visual_contract": {
            "chart_ids": chart_ids,
            "table_ids": table_ids,
            "javascript_modes": ["enabled", "disabled"],
            "viewport_widths": [1440, 390],
        },
        "bundle_files": {
            key: file_identity(path) for key, path in paths.items()
        },
        "artifact_title": artifact["manifest"]["title"],
        "finalizer": file_identity(Path(__file__)),
    }
    sidecar_path = write_ready_last(output_path=output_path, value=ready)
    print(
        json.dumps(
            {
                "all_pass": True,
                "output": str(output_path),
                "sidecar": str(sidecar_path),
                "datasets": len(CANONICAL_DATASETS),
                "chromosome_rows": len(data["chromosomes"]),
                "charts": len(chart_ids),
                "tables": len(table_ids),
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
