#!/usr/bin/env python3
"""Select pre-gated confirmatory permutation-floor units for refinement."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any, Sequence


SCHEMA_NAME = "intersubmod.pattern_methyl_adaptive_refinement_selection"
SCHEMA_VERSION = "1.0.0"
FROZEN_SCREEN_PERMUTATIONS = 999
FROZEN_REFINEMENT_PERMUTATIONS = 49999
KEY_FIELDS = ("dataset", "chrom", "region_id", "hp_raw")
OUTPUT_FIELDS = KEY_FIELDS + (
    "screen_p",
    "screen_r2",
    "best_pair_distance_contrast",
    "best_pair_standardized_effect",
    "permdisp_p",
    "max_geometry_smd",
    "equal_n_retention",
    "rarefaction_retention",
    "refinement_permutations",
)


class RefinementSelectionError(RuntimeError):
    """Raised when screen evidence violates the frozen selection contract."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def truth(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def finite_float(value: object, field: str) -> float:
    try:
        number = float(str(value))
    except (TypeError, ValueError) as exc:
        raise RefinementSelectionError(f"{field} is not numeric: {value!r}") from exc
    if not math.isfinite(number):
        raise RefinementSelectionError(f"{field} is not finite: {value!r}")
    return number


def positive_integer(value: object, field: str) -> int:
    try:
        number = int(str(value))
    except (TypeError, ValueError) as exc:
        raise RefinementSelectionError(
            f"{field} is not an integer: {value!r}"
        ) from exc
    if number <= 0:
        raise RefinementSelectionError(f"{field} must be positive: {number}")
    return number


def open_evidence(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(
        "r", encoding="utf-8", newline=""
    )


def select_rows(
    evidence_path: Path,
    *,
    screen_permutations: int,
    refinement_permutations: int,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    if screen_permutations != FROZEN_SCREEN_PERMUTATIONS:
        raise RefinementSelectionError(
            "screen permutations must equal the frozen budget "
            f"{FROZEN_SCREEN_PERMUTATIONS}: {screen_permutations}"
        )
    if refinement_permutations != FROZEN_REFINEMENT_PERMUTATIONS:
        raise RefinementSelectionError(
            "refinement permutations must equal the frozen budget "
            f"{FROZEN_REFINEMENT_PERMUTATIONS}: {refinement_permutations}"
        )
    required = set(KEY_FIELDS) | {
        "pair_full4",
        "k_ge_3",
        "evaluation_status",
        "permanova_p",
        "permanova_r2",
        "permanova_permutations_requested",
        "permanova_permutations_realized",
        "permdisp_p",
        "best_pair_distance_contrast",
        "best_pair_standardized_effect",
        "max_geometry_smd",
        "all_states_n8",
        "equal_n_retention",
        "rarefaction_retention",
        "multiplicity_family",
    }
    selected: list[dict[str, Any]] = []
    total = 0
    evaluable = 0
    floor = 1.0 / (screen_permutations + 1)
    seen: set[tuple[str, str, str, str]] = set()
    with open_evidence(evidence_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise RefinementSelectionError(
                f"screen evidence missing fields: "
                f"{sorted(required - set(reader.fieldnames or []))}"
            )
        for line_number, row in enumerate(reader, start=2):
            total += 1
            key = tuple(row[field].strip() for field in KEY_FIELDS)
            if any(not value for value in key) or key in seen:
                raise RefinementSelectionError(
                    f"invalid or duplicate exact key at line {line_number}: {key}"
                )
            seen.add(key)
            confirmatory = truth(row["pair_full4"]) or truth(row["k_ge_3"])
            if (
                not confirmatory
                or row["multiplicity_family"] != "CONFIRMATORY_FULL4_OR_LONG"
            ):
                raise RefinementSelectionError(
                    f"non-confirmatory row in screen evidence: {key}"
                )
            if row["evaluation_status"] != "EVALUABLE":
                continue
            evaluable += 1
            requested = positive_integer(
                row["permanova_permutations_requested"],
                "permanova_permutations_requested",
            )
            realized = positive_integer(
                row["permanova_permutations_realized"],
                "permanova_permutations_realized",
            )
            if requested != screen_permutations or realized != screen_permutations:
                raise RefinementSelectionError(
                    f"screen permutation mismatch at {key}: "
                    f"{requested}/{realized}"
                )
            p_value = finite_float(row["permanova_p"], "permanova_p")
            metrics = {
                "screen_r2": finite_float(row["permanova_r2"], "permanova_r2"),
                "best_pair_distance_contrast": finite_float(
                    row["best_pair_distance_contrast"],
                    "best_pair_distance_contrast",
                ),
                "best_pair_standardized_effect": finite_float(
                    row["best_pair_standardized_effect"],
                    "best_pair_standardized_effect",
                ),
                "permdisp_p": finite_float(row["permdisp_p"], "permdisp_p"),
                "max_geometry_smd": finite_float(
                    row["max_geometry_smd"], "max_geometry_smd"
                ),
                "equal_n_retention": finite_float(
                    row["equal_n_retention"], "equal_n_retention"
                ),
                "rarefaction_retention": finite_float(
                    row["rarefaction_retention"], "rarefaction_retention"
                ),
            }
            eligible = (
                math.isclose(p_value, floor, rel_tol=0.0, abs_tol=1e-12)
                and metrics["screen_r2"] >= 0.10
                and metrics["best_pair_distance_contrast"] >= 0.10
                and metrics["best_pair_standardized_effect"] >= 0.50
                and metrics["permdisp_p"] >= 0.05
                and metrics["max_geometry_smd"] < 0.50
                and truth(row["all_states_n8"])
                and metrics["equal_n_retention"] >= 0.50
                and metrics["rarefaction_retention"] >= 0.50
            )
            if eligible:
                selected.append(
                    {
                        **dict(zip(KEY_FIELDS, key)),
                        "screen_p": p_value,
                        **metrics,
                        "refinement_permutations": refinement_permutations,
                    }
                )
    selected.sort(
        key=lambda row: (
            row["dataset"],
            int(row["chrom"][3:]),
            row["region_id"],
            row["hp_raw"],
        )
    )
    return selected, {
        "screen_units": total,
        "screen_evaluable": evaluable,
        "selected_for_refinement": len(selected),
    }


def output_value(value: object) -> str:
    return f"{value:.10g}" if isinstance(value, float) else str(value)


def execute(
    evidence_path: Path,
    output_tsv: Path,
    output_json: Path,
    *,
    screen_permutations: int,
    refinement_permutations: int,
) -> dict[str, Any]:
    if not evidence_path.is_file():
        raise RefinementSelectionError(f"screen evidence missing: {evidence_path}")
    for path in (output_tsv, output_json):
        if path.exists():
            raise RefinementSelectionError(f"refusing to overwrite output: {path}")
    rows, counts = select_rows(
        evidence_path,
        screen_permutations=screen_permutations,
        refinement_permutations=refinement_permutations,
    )
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: output_value(row[field]) for field in OUTPUT_FIELDS})
    summary = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "all_pass": True,
        "contract": {
            "family": "CONFIRMATORY_FULL4_OR_LONG",
            "screen_permutations": screen_permutations,
            "screen_floor": 1.0 / (screen_permutations + 1),
            "refinement_permutations": refinement_permutations,
            "pre_gates": {
                "permanova_r2_min": 0.10,
                "best_pair_distance_contrast_min": 0.10,
                "best_pair_standardized_effect_min": 0.50,
                "permdisp_p_min": 0.05,
                "max_geometry_smd_strict_max": 0.50,
                "all_states_n8": True,
                "equal_n_retention_min": 0.50,
                "rarefaction_retention_min": 0.50,
            },
        },
        "counts": counts,
        "inputs": {
            "screen_evidence": {
                "path": str(evidence_path.resolve()),
                "sha256": sha256_file(evidence_path),
                "size_bytes": evidence_path.stat().st_size,
            }
        },
        "outputs": {
            "unit_keys": {
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
    parser.add_argument("--screen-evidence", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument(
        "--screen-permutations",
        type=int,
        default=FROZEN_SCREEN_PERMUTATIONS,
    )
    parser.add_argument(
        "--refinement-permutations",
        type=int,
        default=FROZEN_REFINEMENT_PERMUTATIONS,
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    summary = execute(
        args.screen_evidence,
        args.output_tsv,
        args.output_json,
        screen_permutations=args.screen_permutations,
        refinement_permutations=args.refinement_permutations,
    )
    print(json.dumps(summary["counts"], ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RefinementSelectionError as exc:
        print(f"FAIL_CLOSED: {exc}", file=sys.stderr)
        raise SystemExit(2)
