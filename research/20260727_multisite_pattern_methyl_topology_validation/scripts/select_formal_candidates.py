#!/usr/bin/env python3
"""Freeze formal_n5 units and their deduplicated marker universe."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Sequence


SCHEMA_NAME = "intersubmod.pattern_methyl_formal_candidate_selection"
SCHEMA_VERSION = "1.0.0"


class SelectionError(RuntimeError):
    pass


def truth(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fields,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)
        Path(temporary_name).replace(path)
    except Exception:
        temporary = Path(temporary_name)
        if temporary.exists():
            archive = path.parent / "_failed_staging_archive"
            archive.mkdir(parents=True, exist_ok=True)
            temporary.replace(archive / temporary.name)
        raise


def select_candidates(
    pattern_counts: Path, marker_universe: Path
) -> tuple[list[str], list[dict[str, str]], list[dict[str, str]], dict[str, int]]:
    with pattern_counts.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "formal_n5" not in reader.fieldnames:
            raise SelectionError("pattern counts missing formal_n5")
        count_fields = list(reader.fieldnames)
        formal_rows = [dict(row) for row in reader if truth(row["formal_n5"])]
    if not formal_rows:
        raise SelectionError("no formal_n5 rows")
    formal_regions = {
        (row["dataset"], row["chrom"], row["region_id"]) for row in formal_rows
    }

    marker_records: dict[tuple[str, str, int], dict[str, str]] = {}
    with marker_universe.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"dataset", "chrom", "region_id", "position1"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise SelectionError(
                f"marker universe missing fields: {sorted(required - set(reader.fieldnames or []))}"
            )
        for row in reader:
            region_key = (row["dataset"], row["chrom"], row["region_id"])
            if region_key not in formal_regions:
                continue
            marker_key = (row["dataset"], row["chrom"], int(row["position1"]))
            marker_records.setdefault(
                marker_key,
                {
                    "dataset": row["dataset"],
                    "chrom": row["chrom"],
                    "position1": row["position1"],
                },
            )
    if not marker_records:
        raise SelectionError("formal units resolved to no markers")
    marker_rows = [
        marker_records[key]
        for key in sorted(
            marker_records,
            key=lambda item: (
                item[0],
                int(item[1][3:]) if item[1].startswith("chr") else 10**9,
                item[2],
            ),
        )
    ]
    strata = Counter(
        "pair_full4"
        if truth(row.get("pair_full4"))
        else "long_k_ge_3"
        if truth(row.get("k_ge_3"))
        else "secondary_pair"
        for row in formal_rows
    )
    summary = {
        "formal_units": len(formal_rows),
        "formal_regions": len(formal_regions),
        "unique_markers": len(marker_rows),
        **dict(sorted(strata.items())),
    }
    return count_fields, formal_rows, marker_rows, summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pattern-counts", type=Path, required=True)
    parser.add_argument("--marker-universe", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    fields, formal_rows, marker_rows, counts = select_candidates(
        args.pattern_counts, args.marker_universe
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    formal_path = args.output_dir / "formal_pattern_counts.tsv"
    markers_path = args.output_dir / "formal_marker_universe.tsv"
    receipt_path = args.output_dir / "formal_selection.receipt.json"
    for path in (formal_path, markers_path, receipt_path):
        if path.exists():
            raise SelectionError(f"refusing to overwrite: {path}")
    atomic_write_tsv(formal_path, fields, formal_rows)
    atomic_write_tsv(markers_path, ["dataset", "chrom", "position1"], marker_rows)
    receipt = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "all_pass": True,
        "counts": counts,
        "inputs": {
            "pattern_counts": {
                "path": str(args.pattern_counts.resolve()),
                "sha256": sha256_file(args.pattern_counts),
            },
            "marker_universe": {
                "path": str(args.marker_universe.resolve()),
                "sha256": sha256_file(args.marker_universe),
            },
        },
        "outputs": {
            "formal_pattern_counts": {
                "path": str(formal_path.resolve()),
                "sha256": sha256_file(formal_path),
            },
            "formal_marker_universe": {
                "path": str(markers_path.resolve()),
                "sha256": sha256_file(markers_path),
            },
        },
    }
    receipt_path.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(counts, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SelectionError as exc:
        print(f"FAIL_CLOSED: {exc}", file=__import__("sys").stderr)
        raise SystemExit(2)
