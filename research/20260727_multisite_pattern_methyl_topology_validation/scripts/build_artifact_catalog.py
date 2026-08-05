#!/usr/bin/env python3
"""Build a source-bound catalog for the 2026-07-15 all-sSNV artifacts.

The catalog resolves each unique (dataset, chrom, position1) marker to the
frozen reads, methylation, and BERNOULLI matrix artifacts.  The ``hp`` column
in reads.tsv is deliberately ignored: exact raw HP authority comes from the
separate frozen LongPhase-S read-call sidecar used by the downstream join.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCHEMA_NAME = "intersubmod.multisite_pattern_methyl_artifact_catalog"
SCHEMA_VERSION = "1.0.0"
ANALYSIS_SUFFIX = ".longphase_s.recalibrated.pass.autosomal_biallelic_snv"
DEFAULT_PRIMARY_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/"
    "intersubmod_all_ssnv_v2_verification_fix"
)
CANONICAL_CHROMS = frozenset(f"chr{i}" for i in range(1, 23))
SAFE_DATASET_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
SHA_CHUNK_BYTES = 8 * 1024 * 1024

TSV_COLUMNS = [
    "schema_version",
    "marker_index",
    "input_occurrences",
    "dataset",
    "chrom",
    "position1",
    "status",
    "status_detail",
    "region_path",
    "reads_path",
    "reads_size_bytes",
    "reads_sha256",
    "reads_row_count",
    "reads_hp_column_present",
    "reads_hp_authoritative",
    "methylation_path",
    "methylation_size_bytes",
    "methylation_sha256",
    "methylation_row_count",
    "methylation_cpg_count",
    "bernoulli_path",
    "bernoulli_size_bytes",
    "bernoulli_sha256",
    "bernoulli_row_count",
    "bernoulli_column_count",
    "receipt_path",
    "receipt_size_bytes",
    "receipt_sha256",
    "receipt_schema_name",
    "receipt_schema_version",
    "receipt_pass",
    "artifact_bundle_sha256",
]


class CatalogError(RuntimeError):
    """A machine-classifiable fail-closed catalog error."""

    def __init__(self, code: str, detail: str):
        super().__init__(f"{code}: {detail}")
        self.code = code
        self.detail = detail


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(SHA_CHUNK_BYTES)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def _file_binding(path: Path, primary_root: Path) -> dict[str, Any]:
    try:
        resolved = path.resolve(strict=True)
    except FileNotFoundError as exc:
        raise CatalogError("ARTIFACT_MISSING", str(path)) from exc
    try:
        resolved.relative_to(primary_root)
    except ValueError as exc:
        raise CatalogError("PATH_ESCAPES_PRIMARY_ROOT", str(resolved)) from exc
    if not resolved.is_file():
        raise CatalogError("ARTIFACT_NOT_REGULAR_FILE", str(resolved))
    size = resolved.stat().st_size
    if size <= 0:
        raise CatalogError("ARTIFACT_EMPTY", str(resolved))
    return {
        "path": str(resolved),
        "size_bytes": size,
        "sha256": _sha256_file(resolved),
    }


def _require_unique_nonempty(value: str, seen: set[str], label: str) -> None:
    if not value:
        raise CatalogError("ARTIFACT_SCHEMA_INVALID", f"empty {label}")
    if value in seen:
        raise CatalogError("ARTIFACT_SCHEMA_INVALID", f"duplicate {label}: {value}")
    seen.add(value)


def _scan_reads(path: Path) -> dict[str, Any]:
    required = {
        "read_id",
        "read_name",
        "chr",
        "start",
        "end",
        "mapq",
        "alt_support",
        "is_tumor",
        "strand",
    }
    read_ids: list[str] = []
    seen: set[str] = set()
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames
            if not fields or len(fields) != len(set(fields)):
                raise CatalogError("ARTIFACT_SCHEMA_INVALID", "reads header missing or duplicated")
            missing = sorted(required.difference(fields))
            if missing:
                raise CatalogError(
                    "ARTIFACT_SCHEMA_INVALID", f"reads missing columns: {','.join(missing)}"
                )
            for row_number, row in enumerate(reader, start=2):
                if None in row:
                    raise CatalogError(
                        "ARTIFACT_SCHEMA_INVALID", f"reads row {row_number} has extra cells"
                    )
                read_id = (row.get("read_id") or "").strip()
                _require_unique_nonempty(read_id, seen, "reads read_id")
                read_ids.append(read_id)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise CatalogError("ARTIFACT_PARSE_ERROR", f"reads: {exc}") from exc
    if not read_ids:
        raise CatalogError("ARTIFACT_SCHEMA_INVALID", "reads contains no data rows")
    return {
        "row_count": len(read_ids),
        "read_ids": read_ids,
        "hp_column_present": "hp" in fields,
    }


def _scan_methylation(path: Path) -> dict[str, Any]:
    read_ids: list[str] = []
    seen: set[str] = set()
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            try:
                header = next(reader)
            except StopIteration as exc:
                raise CatalogError("ARTIFACT_SCHEMA_INVALID", "methylation is empty") from exc
            if len(header) < 2 or header[0] != "read_id":
                raise CatalogError(
                    "ARTIFACT_SCHEMA_INVALID", "methylation header must start with read_id"
                )
            cpg_labels = header[1:]
            if len(cpg_labels) != len(set(cpg_labels)):
                raise CatalogError("ARTIFACT_SCHEMA_INVALID", "duplicate methylation CpG columns")
            try:
                cpg_positions = [int(value) for value in cpg_labels]
            except ValueError as exc:
                raise CatalogError(
                    "ARTIFACT_SCHEMA_INVALID", "non-integer methylation CpG column"
                ) from exc
            if any(position <= 0 for position in cpg_positions):
                raise CatalogError(
                    "ARTIFACT_SCHEMA_INVALID", "methylation CpG positions must be positive"
                )
            for row_number, row in enumerate(reader, start=2):
                if len(row) != len(header):
                    raise CatalogError(
                        "ARTIFACT_SCHEMA_INVALID",
                        f"methylation row {row_number} width {len(row)} != {len(header)}",
                    )
                read_id = row[0].strip()
                _require_unique_nonempty(read_id, seen, "methylation read_id")
                read_ids.append(read_id)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise CatalogError("ARTIFACT_PARSE_ERROR", f"methylation: {exc}") from exc
    if not read_ids:
        raise CatalogError("ARTIFACT_SCHEMA_INVALID", "methylation contains no data rows")
    return {
        "row_count": len(read_ids),
        "cpg_count": len(cpg_labels),
        "read_ids": read_ids,
    }


def _scan_bernoulli(path: Path) -> dict[str, Any]:
    row_ids: list[str] = []
    seen: set[str] = set()
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            try:
                header = next(reader)
            except StopIteration as exc:
                raise CatalogError("ARTIFACT_SCHEMA_INVALID", "BERNOULLI matrix is empty") from exc
            if len(header) < 2 or header[0] != "read_id":
                raise CatalogError(
                    "ARTIFACT_SCHEMA_INVALID", "BERNOULLI header must start with read_id"
                )
            column_ids = [value.strip() for value in header[1:]]
            if any(not value for value in column_ids) or len(column_ids) != len(set(column_ids)):
                raise CatalogError(
                    "ARTIFACT_SCHEMA_INVALID", "BERNOULLI column IDs are empty or duplicated"
                )
            for row_number, row in enumerate(reader, start=2):
                if len(row) != len(header):
                    raise CatalogError(
                        "ARTIFACT_SCHEMA_INVALID",
                        f"BERNOULLI row {row_number} width {len(row)} != {len(header)}",
                    )
                read_id = row[0].strip()
                _require_unique_nonempty(read_id, seen, "BERNOULLI row read_id")
                row_ids.append(read_id)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise CatalogError("ARTIFACT_PARSE_ERROR", f"BERNOULLI: {exc}") from exc
    if not row_ids:
        raise CatalogError("ARTIFACT_SCHEMA_INVALID", "BERNOULLI contains no data rows")
    if len(row_ids) != len(column_ids) or set(row_ids) != set(column_ids):
        raise CatalogError(
            "BERNOULLI_SHAPE_INVALID",
            f"rows={len(row_ids)} columns={len(column_ids)} or ID sets differ",
        )
    return {
        "row_count": len(row_ids),
        "column_count": len(column_ids),
        "read_ids": row_ids,
    }


def _command_option(command: Sequence[str], option: str) -> str | None:
    try:
        index = command.index(option)
    except ValueError:
        return None
    if index + 1 >= len(command):
        return None
    return command[index + 1]


def _load_receipt(
    dataset: str, dataset_dir: Path, primary_root: Path, window_bp: int
) -> dict[str, Any]:
    binding = _file_binding(dataset_dir / "run_receipt.json", primary_root)
    receipt_path = Path(binding["path"])
    try:
        payload = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise CatalogError("RECEIPT_INVALID", str(exc)) from exc
    if not isinstance(payload, dict):
        raise CatalogError("RECEIPT_INVALID", "receipt root is not an object")
    if payload.get("schema_name") != "intersubmod.all_ssnv_site_run":
        raise CatalogError("RECEIPT_INVALID", f"unexpected schema_name={payload.get('schema_name')!r}")
    if payload.get("sample") != dataset:
        raise CatalogError(
            "RECEIPT_INVALID", f"sample={payload.get('sample')!r} expected={dataset!r}"
        )
    if payload.get("pass") is not True:
        raise CatalogError("RECEIPT_INVALID", "receipt pass is not true")
    validation = payload.get("validation")
    if not isinstance(validation, dict) or validation.get("pass") is not True:
        raise CatalogError("RECEIPT_INVALID", "validation.pass is not true")
    try:
        output_dir = Path(str(payload["output_dir"])).resolve(strict=True)
    except (KeyError, FileNotFoundError) as exc:
        raise CatalogError("RECEIPT_INVALID", "output_dir missing or unavailable") from exc
    if output_dir != dataset_dir:
        raise CatalogError(
            "RECEIPT_INVALID", f"output_dir={output_dir} expected={dataset_dir}"
        )
    command = payload.get("command")
    if not isinstance(command, list) or not all(isinstance(item, str) for item in command):
        raise CatalogError("RECEIPT_INVALID", "command is not a string array")
    if _command_option(command, "--distance-metric") != "BERNOULLI":
        raise CatalogError("RECEIPT_INVALID", "command is not bound to BERNOULLI")
    if _command_option(command, "-w") != str(window_bp):
        raise CatalogError("RECEIPT_INVALID", f"command is not bound to window {window_bp}")
    expected = validation.get("expected_vcf_sites")
    counts = [
        validation.get("reads_files"),
        validation.get("methylation_files"),
        validation.get("bernoulli_matrix_files"),
    ]
    if (
        isinstance(expected, bool)
        or not isinstance(expected, int)
        or expected <= 0
        or any(isinstance(value, bool) or not isinstance(value, int) for value in counts)
        or any(value != expected for value in counts)
    ):
        raise CatalogError(
            "RECEIPT_INVALID",
            f"artifact counts are not positive and equal: expected={expected!r}, counts={counts!r}",
        )
    return {
        "receipt_path": binding["path"],
        "receipt_size_bytes": binding["size_bytes"],
        "receipt_sha256": binding["sha256"],
        "receipt_schema_name": payload["schema_name"],
        "receipt_schema_version": payload.get("schema_version", ""),
        "receipt_pass": 1,
    }


def _empty_record(index: int, marker: Mapping[str, Any]) -> dict[str, Any]:
    record = {column: "" for column in TSV_COLUMNS}
    record.update(
        {
            "schema_version": SCHEMA_VERSION,
            "marker_index": index,
            "input_occurrences": marker["input_occurrences"],
            "dataset": marker["dataset"],
            "chrom": marker["chrom"],
            "position1": marker["position1"],
            "status": "NOT_RUN",
            "reads_hp_authoritative": 0,
            "receipt_pass": 0,
        }
    )
    return record


def _bundle_digest(bindings: Iterable[Mapping[str, Any]]) -> str:
    canonical = "\n".join(
        f"{item['path']}\t{item['size_bytes']}\t{item['sha256']}" for item in bindings
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _resolve_marker(
    index: int,
    marker: Mapping[str, Any],
    primary_root: Path,
    window_bp: int,
    receipt_cache: dict[str, dict[str, Any] | CatalogError],
) -> dict[str, Any]:
    record = _empty_record(index, marker)
    dataset = str(marker["dataset"])
    chrom = str(marker["chrom"])
    position1 = int(marker["position1"])
    try:
        dataset_dir = (primary_root / dataset).resolve(strict=True)
        try:
            dataset_dir.relative_to(primary_root)
        except ValueError as exc:
            raise CatalogError("PATH_ESCAPES_PRIMARY_ROOT", str(dataset_dir)) from exc
        if not dataset_dir.is_dir():
            raise CatalogError("DATASET_DIR_INVALID", str(dataset_dir))

        if dataset not in receipt_cache:
            try:
                receipt_cache[dataset] = _load_receipt(
                    dataset, dataset_dir, primary_root, window_bp
                )
            except CatalogError as exc:
                receipt_cache[dataset] = exc
        receipt = receipt_cache[dataset]
        if isinstance(receipt, CatalogError):
            raise receipt
        record.update(receipt)

        analysis_dir = dataset_dir / f"{dataset}{ANALYSIS_SUFFIX}"
        marker_dir = analysis_dir / chrom / f"{chrom}_{position1}"
        region_dir = marker_dir / f"{chrom}_{position1 - window_bp}_{position1 + window_bp}"
        try:
            region_resolved = region_dir.resolve(strict=True)
        except FileNotFoundError as exc:
            raise CatalogError("REGION_DIR_MISSING", str(region_dir)) from exc
        try:
            region_resolved.relative_to(primary_root)
        except ValueError as exc:
            raise CatalogError("PATH_ESCAPES_PRIMARY_ROOT", str(region_resolved)) from exc
        if not region_resolved.is_dir():
            raise CatalogError("REGION_DIR_INVALID", str(region_resolved))
        record["region_path"] = str(region_resolved)

        reads_path = region_resolved / "reads" / "reads.tsv"
        methylation_path = region_resolved / "methylation" / "methylation.csv"
        bernoulli_path = region_resolved / "distance" / "BERNOULLI" / "matrix.csv"
        reads_binding = _file_binding(reads_path, primary_root)
        methylation_binding = _file_binding(methylation_path, primary_root)
        bernoulli_binding = _file_binding(bernoulli_path, primary_root)

        reads = _scan_reads(Path(reads_binding["path"]))
        methylation = _scan_methylation(Path(methylation_binding["path"]))
        bernoulli = _scan_bernoulli(Path(bernoulli_binding["path"]))
        read_set = set(reads["read_ids"])
        if set(methylation["read_ids"]) != read_set:
            raise CatalogError(
                "READ_ID_MISMATCH", "reads.tsv and methylation.csv read IDs differ"
            )
        if set(bernoulli["read_ids"]) != read_set:
            raise CatalogError(
                "READ_ID_MISMATCH", "reads.tsv and BERNOULLI matrix read IDs differ"
            )
        if not (
            reads["row_count"]
            == methylation["row_count"]
            == bernoulli["row_count"]
            == bernoulli["column_count"]
        ):
            raise CatalogError(
                "ROW_COUNT_MISMATCH",
                "reads, methylation, and BERNOULLI dimensions are not identical",
            )

        record.update(
            {
                "status": "PASS",
                "status_detail": "",
                "reads_path": reads_binding["path"],
                "reads_size_bytes": reads_binding["size_bytes"],
                "reads_sha256": reads_binding["sha256"],
                "reads_row_count": reads["row_count"],
                "reads_hp_column_present": int(reads["hp_column_present"]),
                "reads_hp_authoritative": 0,
                "methylation_path": methylation_binding["path"],
                "methylation_size_bytes": methylation_binding["size_bytes"],
                "methylation_sha256": methylation_binding["sha256"],
                "methylation_row_count": methylation["row_count"],
                "methylation_cpg_count": methylation["cpg_count"],
                "bernoulli_path": bernoulli_binding["path"],
                "bernoulli_size_bytes": bernoulli_binding["size_bytes"],
                "bernoulli_sha256": bernoulli_binding["sha256"],
                "bernoulli_row_count": bernoulli["row_count"],
                "bernoulli_column_count": bernoulli["column_count"],
                "artifact_bundle_sha256": _bundle_digest(
                    [
                        reads_binding,
                        methylation_binding,
                        bernoulli_binding,
                        {
                            "path": receipt["receipt_path"],
                            "size_bytes": receipt["receipt_size_bytes"],
                            "sha256": receipt["receipt_sha256"],
                        },
                    ]
                ),
            }
        )
    except FileNotFoundError as exc:
        record["status"] = "FAIL_DATASET_DIR_MISSING"
        record["status_detail"] = str(exc)
    except CatalogError as exc:
        record["status"] = f"FAIL_{exc.code}"
        record["status_detail"] = exc.detail
    except (OSError, ValueError) as exc:
        record["status"] = "FAIL_UNEXPECTED_IO_OR_VALUE_ERROR"
        record["status_detail"] = str(exc)
    return record


def _load_markers(path: Path) -> tuple[list[dict[str, Any]], int]:
    required = {"dataset", "chrom", "position1"}
    unique: "OrderedDict[tuple[str, str, int], dict[str, Any]]" = OrderedDict()
    source_rows = 0
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames
            if not fields or len(fields) != len(set(fields)):
                raise CatalogError("MARKER_TSV_INVALID", "header missing or duplicated")
            missing = sorted(required.difference(fields))
            if missing:
                raise CatalogError(
                    "MARKER_TSV_INVALID", f"missing columns: {','.join(missing)}"
                )
            for row_number, row in enumerate(reader, start=2):
                source_rows += 1
                if None in row:
                    raise CatalogError(
                        "MARKER_TSV_INVALID", f"row {row_number} has extra cells"
                    )
                dataset = (row.get("dataset") or "").strip()
                chrom = (row.get("chrom") or "").strip()
                raw_position = (row.get("position1") or "").strip()
                if not SAFE_DATASET_RE.fullmatch(dataset):
                    raise CatalogError(
                        "MARKER_TSV_INVALID", f"row {row_number} unsafe dataset={dataset!r}"
                    )
                if chrom not in CANONICAL_CHROMS:
                    raise CatalogError(
                        "MARKER_TSV_INVALID", f"row {row_number} invalid chrom={chrom!r}"
                    )
                try:
                    position1 = int(raw_position)
                except ValueError as exc:
                    raise CatalogError(
                        "MARKER_TSV_INVALID",
                        f"row {row_number} invalid position1={raw_position!r}",
                    ) from exc
                if position1 <= 0:
                    raise CatalogError(
                        "MARKER_TSV_INVALID",
                        f"row {row_number} position1 must be positive",
                    )
                key = (dataset, chrom, position1)
                if key not in unique:
                    unique[key] = {
                        "dataset": dataset,
                        "chrom": chrom,
                        "position1": position1,
                        "input_occurrences": 1,
                    }
                else:
                    unique[key]["input_occurrences"] += 1
    except CatalogError:
        raise
    except (OSError, UnicodeError, csv.Error) as exc:
        raise CatalogError("MARKER_TSV_INVALID", str(exc)) from exc
    if source_rows == 0:
        raise CatalogError("MARKER_TSV_INVALID", "no marker rows")
    return list(unique.values()), source_rows


def build_catalog(
    marker_tsv: Path | str,
    primary_root: Path | str = DEFAULT_PRIMARY_ROOT,
    *,
    window_bp: int = 5000,
    created_at_utc: str | None = None,
) -> dict[str, Any]:
    """Resolve markers and return the complete JSON catalog document."""

    marker_path = Path(marker_tsv)
    root_path = Path(primary_root)
    if isinstance(window_bp, bool) or not isinstance(window_bp, int) or window_bp <= 0:
        raise CatalogError("ARGUMENT_INVALID", "window_bp must be a positive integer")
    try:
        marker_resolved = marker_path.resolve(strict=True)
    except FileNotFoundError as exc:
        raise CatalogError("MARKER_TSV_INVALID", str(marker_path)) from exc
    try:
        root_resolved = root_path.resolve(strict=True)
    except FileNotFoundError as exc:
        raise CatalogError("PRIMARY_ROOT_INVALID", str(root_path)) from exc
    if not root_resolved.is_dir():
        raise CatalogError("PRIMARY_ROOT_INVALID", str(root_resolved))

    markers, source_rows = _load_markers(marker_resolved)
    receipt_cache: dict[str, dict[str, Any] | CatalogError] = {}
    records = [
        _resolve_marker(index, marker, root_resolved, window_bp, receipt_cache)
        for index, marker in enumerate(markers, start=1)
    ]
    passed = sum(record["status"] == "PASS" for record in records)
    failed = len(records) - passed
    return {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": created_at_utc
        or datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "contract": {
            "primary_root": str(root_resolved),
            "analysis_suffix": ANALYSIS_SUFFIX,
            "window_bp": window_bp,
            "artifact_paths": {
                "reads": "reads/reads.tsv",
                "methylation": "methylation/methylation.csv",
                "bernoulli": "distance/BERNOULLI/matrix.csv",
                "receipt": "<dataset>/run_receipt.json",
            },
            "reads_hp_authority": "NOT_AUTHORITATIVE_AND_NOT_CONSUMED",
            "exact_raw_hp_authority": "SEPARATE_FROZEN_LONGPHASE_S_READ_CALL_SIDECAR",
            "fail_closed": True,
        },
        "source": {
            "marker_tsv": str(marker_resolved),
            "marker_tsv_size_bytes": marker_resolved.stat().st_size,
            "marker_tsv_sha256": _sha256_file(marker_resolved),
            "input_rows": source_rows,
            "unique_markers": len(markers),
            "duplicate_rows_collapsed": source_rows - len(markers),
        },
        "summary": {
            "markers_total": len(records),
            "markers_pass": passed,
            "markers_fail": failed,
        },
        "pass": failed == 0,
        "status": "PASS" if failed == 0 else "FAIL_CLOSED",
        "records": records,
    }


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    )
    temporary = Path(handle.name)
    try:
        with handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        if temporary.exists():
            archive = path.parent / "_failed_staging_archive"
            archive.mkdir(parents=True, exist_ok=True)
            os.replace(temporary, archive / temporary.name)
        raise


def write_catalog(document: Mapping[str, Any], output_tsv: Path | str, output_json: Path | str) -> None:
    """Write JSON and TSV projections atomically per file."""

    tsv_path = Path(output_tsv)
    json_path = Path(output_json)
    if tsv_path.resolve() == json_path.resolve():
        raise CatalogError("ARGUMENT_INVALID", "TSV and JSON output paths must differ")

    from io import StringIO

    buffer = StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=TSV_COLUMNS,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    for record in document["records"]:
        writer.writerow(record)
    json_text = json.dumps(document, ensure_ascii=False, indent=2, sort_keys=False) + "\n"
    _atomic_write_text(tsv_path, buffer.getvalue())
    _atomic_write_text(json_path, json_text)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve marker TSV rows to source-bound reads/methylation/BERNOULLI "
            "artifacts and emit v1 TSV/JSON catalogs."
        )
    )
    parser.add_argument("--markers", required=True, type=Path, help="Input TSV")
    parser.add_argument(
        "--primary-root",
        type=Path,
        default=DEFAULT_PRIMARY_ROOT,
        help=f"All-sSNV primary root (default: {DEFAULT_PRIMARY_ROOT})",
    )
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--window-bp", type=int, default=5000)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        document = build_catalog(
            args.markers,
            args.primary_root,
            window_bp=args.window_bp,
        )
        write_catalog(document, args.output_tsv, args.output_json)
    except CatalogError as exc:
        print(f"FAIL_CLOSED {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"FAIL_CLOSED OUTPUT_IO_ERROR: {exc}", file=sys.stderr)
        return 2
    summary = document["summary"]
    print(
        json.dumps(
            {
                "status": document["status"],
                **summary,
                "output_tsv": str(args.output_tsv.resolve()),
                "output_json": str(args.output_json.resolve()),
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if document["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
