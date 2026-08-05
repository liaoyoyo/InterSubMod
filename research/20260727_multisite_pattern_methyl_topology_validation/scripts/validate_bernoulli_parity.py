#!/usr/bin/env python3
"""Validate recomputed BERNOULLI distances against frozen marker matrices.

Every formal marker is checked.  To bound quadratic work, each marker uses a
deterministic SHA-256-selected subset of read IDs while retaining every CpG.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np


SCHEMA_NAME = "intersubmod.bernoulli_artifact_parity"
SCHEMA_VERSION = "1.0.0"
OUTPUT_FIELDS = (
    "dataset",
    "chrom",
    "position1",
    "status",
    "selected_reads",
    "cpg_count",
    "pair_count",
    "invalid_mask_mismatches",
    "max_absolute_error",
    "mean_absolute_error",
)


class ParityContractError(RuntimeError):
    """Raised when an input or matrix violates the parity contract."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def load_analysis_module() -> Any:
    script = Path(__file__).resolve().with_name("analyze_pattern_methylation.py")
    module_name = "_intersubmod_pattern_methyl_analysis_for_parity"
    loaded = sys.modules.get(module_name)
    if loaded is not None:
        return loaded
    spec = importlib.util.spec_from_file_location(module_name, script)
    if spec is None or spec.loader is None:
        raise ParityContractError(f"cannot import analyzer: {script}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


ANALYSIS = load_analysis_module()


def selected_read_ids(
    methylation_path: Path,
    marker_key: tuple[str, str, str],
    max_reads: int,
) -> list[str]:
    read_ids: list[str] = []
    with methylation_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ParityContractError(
                f"empty methylation file: {methylation_path}"
            ) from exc
        if len(header) < 2 or header[0] != "read_id":
            raise ParityContractError(
                f"invalid methylation header: {methylation_path}"
            )
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ParityContractError(
                    f"methylation width mismatch: {methylation_path}:{line_number}"
                )
            read_ids.append(row[0])
    if len(read_ids) != len(set(read_ids)):
        raise ParityContractError(
            f"duplicate methylation read IDs: {methylation_path}"
        )
    if len(read_ids) < 2:
        raise ParityContractError(
            f"fewer than two methylation reads: {methylation_path}"
        )
    prefix = "\x1f".join(marker_key)
    ranked = sorted(
        read_ids,
        key=lambda read_id: (
            hashlib.sha256(f"{prefix}\x1f{read_id}".encode("utf-8")).digest(),
            read_id,
        ),
    )
    return ranked[: min(max_reads, len(ranked))]


def load_selected_methylation(
    path: Path, selected: Sequence[str]
) -> tuple[np.ndarray, int]:
    selected_set = set(selected)
    values_by_id: dict[str, list[float]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        for line_number, row in enumerate(reader, start=2):
            if row[0] not in selected_set:
                continue
            values: list[float] = []
            for raw in row[1:]:
                if raw in {"", "NA", "NaN", "nan", "."}:
                    values.append(math.nan)
                    continue
                try:
                    probability = float(raw)
                except ValueError as exc:
                    raise ParityContractError(
                        f"invalid methylation value: {path}:{line_number}"
                    ) from exc
                if not 0.0 <= probability <= 1.0:
                    raise ParityContractError(
                        f"methylation outside [0,1]: {path}:{line_number}"
                    )
                values.append(probability)
            values_by_id[row[0]] = values
    missing = set(selected) - set(values_by_id)
    if missing:
        raise ParityContractError(
            f"selected methylation IDs missing: {path} {sorted(missing)[:3]}"
        )
    return np.asarray([values_by_id[read_id] for read_id in selected]), len(header) - 1


def load_selected_bernoulli(path: Path, selected: Sequence[str]) -> np.ndarray:
    selected_set = set(selected)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ParityContractError(f"empty BERNOULLI file: {path}") from exc
        if len(header) < 2 or header[0] != "read_id":
            raise ParityContractError(f"invalid BERNOULLI header: {path}")
        column_index = {read_id: index + 1 for index, read_id in enumerate(header[1:])}
        if len(column_index) != len(header) - 1:
            raise ParityContractError(f"duplicate BERNOULLI columns: {path}")
        missing_columns = selected_set - set(column_index)
        if missing_columns:
            raise ParityContractError(
                f"selected BERNOULLI columns missing: {path} "
                f"{sorted(missing_columns)[:3]}"
            )
        row_values: dict[str, list[float]] = {}
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ParityContractError(
                    f"BERNOULLI width mismatch: {path}:{line_number}"
                )
            if row[0] not in selected_set:
                continue
            try:
                row_values[row[0]] = [
                    (
                        math.nan
                        if row[column_index[read_id]]
                        in {"", "NA", "NaN", "nan", "."}
                        else float(row[column_index[read_id]])
                    )
                    for read_id in selected
                ]
            except ValueError as exc:
                raise ParityContractError(
                    f"invalid BERNOULLI value: {path}:{line_number}"
                ) from exc
    missing_rows = selected_set - set(row_values)
    if missing_rows:
        raise ParityContractError(
            f"selected BERNOULLI rows missing: {path} {sorted(missing_rows)[:3]}"
        )
    return np.asarray([row_values[read_id] for read_id in selected], dtype=float)


def validate_marker(
    row: Mapping[str, str], *, max_reads: int, tolerance: float
) -> dict[str, Any]:
    marker_key = (row["dataset"], row["chrom"], row["position1"])
    methylation_path = Path(row["methylation_path"])
    bernoulli_path = Path(row["bernoulli_path"])
    for path in (methylation_path, bernoulli_path):
        if not path.is_file():
            raise ParityContractError(f"artifact missing: {path}")
    selected = selected_read_ids(methylation_path, marker_key, max_reads)
    methylation, cpg_count = load_selected_methylation(
        methylation_path, selected
    )
    observed = load_selected_bernoulli(bernoulli_path, selected)
    expected = ANALYSIS.bernoulli_distance_matrix(methylation, min_common=3)
    upper = np.triu_indices(len(selected), 1)
    observed_upper = observed[upper]
    expected_upper = expected[upper]
    observed_valid = np.isfinite(observed_upper) & (observed_upper >= 0.0)
    expected_valid = np.isfinite(expected_upper)
    invalid_mismatches = int(np.count_nonzero(observed_valid != expected_valid))
    jointly_valid = observed_valid & expected_valid
    differences = np.abs(observed_upper[jointly_valid] - expected_upper[jointly_valid])
    max_error = float(differences.max()) if len(differences) else math.nan
    mean_error = float(differences.mean()) if len(differences) else math.nan
    passed = (
        invalid_mismatches == 0
        and len(differences) > 0
        and max_error <= tolerance
    )
    return {
        "dataset": marker_key[0],
        "chrom": marker_key[1],
        "position1": marker_key[2],
        "status": "PASS" if passed else "FAIL",
        "selected_reads": len(selected),
        "cpg_count": cpg_count,
        "pair_count": len(upper[0]),
        "invalid_mask_mismatches": invalid_mismatches,
        "max_absolute_error": max_error,
        "mean_absolute_error": mean_error,
    }


def load_catalog(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "dataset",
            "chrom",
            "position1",
            "status",
            "methylation_path",
            "bernoulli_path",
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ParityContractError(
                f"catalog missing columns: {sorted(required - set(reader.fieldnames or []))}"
            )
        rows = [dict(row) for row in reader]
    if not rows:
        raise ParityContractError("artifact catalog is empty")
    if any(row["status"] != "PASS" for row in rows):
        raise ParityContractError("artifact catalog contains non-PASS markers")
    keys = [(row["dataset"], row["chrom"], row["position1"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise ParityContractError("artifact catalog contains duplicate markers")
    return rows


def output_value(value: object) -> str:
    if isinstance(value, float):
        return "" if not math.isfinite(value) else f"{value:.10g}"
    return str(value)


def execute(
    catalog_path: Path,
    output_tsv: Path,
    output_json: Path,
    *,
    max_reads: int,
    tolerance: float,
) -> dict[str, Any]:
    if max_reads < 2:
        raise ParityContractError("max_reads must be at least two")
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ParityContractError("tolerance must be finite and positive")
    for output in (output_tsv, output_json):
        if output.exists():
            raise ParityContractError(f"refusing to overwrite output: {output}")
    rows = load_catalog(catalog_path)
    results = [
        validate_marker(row, max_reads=max_reads, tolerance=tolerance)
        for row in rows
    ]
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for result in results:
            writer.writerow(
                {field: output_value(result[field]) for field in OUTPUT_FIELDS}
            )
    failed = sum(result["status"] != "PASS" for result in results)
    summary = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "all_pass": failed == 0,
        "contract": {
            "scope": "every formal marker",
            "read_selection": "lowest SHA256(dataset,chrom,position1,read_id)",
            "max_reads_per_marker": max_reads,
            "all_cpgs_retained": True,
            "min_common_cpg": 3,
            "absolute_tolerance": tolerance,
        },
        "counts": {
            "markers_total": len(results),
            "markers_pass": len(results) - failed,
            "markers_fail": failed,
            "pair_cells_checked": sum(result["pair_count"] for result in results),
            "invalid_mask_mismatches": sum(
                result["invalid_mask_mismatches"] for result in results
            ),
        },
        "max_absolute_error": max(
            (
                result["max_absolute_error"]
                for result in results
                if math.isfinite(result["max_absolute_error"])
            ),
            default=None,
        ),
        "inputs": {
            "artifact_catalog": {
                "path": str(catalog_path.resolve()),
                "sha256": sha256_file(catalog_path),
            }
        },
        "outputs": {
            "marker_results": {
                "path": str(output_tsv.resolve()),
                "sha256": sha256_file(output_tsv),
                "size_bytes": output_tsv.stat().st_size,
            }
        },
    }
    output_json.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-catalog", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--max-reads", type=int, default=16)
    parser.add_argument("--tolerance", type=float, default=5e-4)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    summary = execute(
        args.artifact_catalog,
        args.output_tsv,
        args.output_json,
        max_reads=args.max_reads,
        tolerance=args.tolerance,
    )
    print(json.dumps(summary["counts"], ensure_ascii=False, sort_keys=True))
    return 0 if summary["all_pass"] else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ParityContractError as exc:
        print(f"FAIL_CLOSED: {exc}", file=sys.stderr)
        raise SystemExit(2)
