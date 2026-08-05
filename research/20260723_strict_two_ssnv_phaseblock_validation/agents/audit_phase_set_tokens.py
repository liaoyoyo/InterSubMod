#!/usr/bin/env python3
"""Count missing/sentinel phase_set tokens in raw and strict artifacts."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from concurrent.futures import ProcessPoolExecutor
import gzip
import json
from pathlib import Path
import re


SENTINELS = {
    ".",
    "NA",
    "N/A",
    "NAN",
    "NONE",
    "NULL",
    "UNKNOWN",
    "UNK",
    "MISSING",
    "NOT_AVAILABLE",
    "NOT APPLICABLE",
}
POSITIVE_INTEGER = re.compile(r"^[0-9]+$")


def scan(task: tuple[str, str, str]) -> tuple[str, str, Counter[str]]:
    dataset, layer, path_string = task
    path = Path(path_string)
    counts: Counter[str] = Counter()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        index = header.index("phase_set")
        for row in reader:
            counts["rows"] += 1
            token = row[index]
            if token == "":
                counts["empty"] += 1
            else:
                normalized = token.strip().upper()
                if normalized in SENTINELS:
                    counts[f"sentinel:{normalized}"] += 1
                elif not POSITIVE_INTEGER.fullmatch(token) or int(token) <= 0:
                    counts["non_positive_integer"] += 1
    return dataset, layer, counts


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample-root",
        action="append",
        required=True,
        help="DATASET=/absolute/path/to/strict_regions_root",
    )
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    tasks: list[tuple[str, str, str]] = []
    for spec in args.sample_root:
        dataset, path_string = spec.split("=", 1)
        root = Path(path_string)
        for receipt_path in sorted((root / "chromosomes").glob("chr*/receipt.json")):
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            tasks.extend(
                (
                    dataset,
                    layer,
                    source["path"],
                )
                for layer, source in (
                    ("raw_molecule_rows", receipt["inputs"]["molecule_calls"]),
                    ("strict_membership_rows", receipt["outputs"]["membership"]),
                    ("strict_edge_rows", receipt["outputs"]["edges"]),
                )
            )

    aggregate: dict[str, dict[str, Counter[str]]] = {}
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for dataset, layer, counts in executor.map(scan, tasks):
            aggregate.setdefault(dataset, {}).setdefault(layer, Counter()).update(counts)

    output = {
        dataset: {
            layer: dict(sorted(counts.items()))
            for layer, counts in sorted(layers.items())
        }
        for dataset, layers in sorted(aggregate.items())
    }
    print(json.dumps(output, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
