#!/usr/bin/env python3
"""Quantify annotation-stratified HCC1395 topology reproducibility.

The input grain is one exact-coordinate, complete-both HCC1395 technical-pair
region.  Cancer-gene and drug annotations are context labels, never topology
truth.  The stratified permutation null preserves chromosome and region-length
decile so the result asks only whether annotation-labelled regions have a
different reproducibility rate from comparable regions.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import math
from pathlib import Path

import numpy as np


SEED = 20260712
N_PERMUTATIONS = 5000
EXPECTED_ROWS = 5720

FEATURES = [
    ("body_gene_any", "GENCODE gene body"),
    ("cgc_body_any", "COSMIC CGC v104 gene body"),
    ("dgidb_interaction_body_any", "DGIdb interaction gene body"),
    ("dgidb_approved_body_any", "DGIdb approved-drug gene body"),
    (
        "dgidb_approved_antineoplastic_body_any",
        "DGIdb approved + antineoplastic gene body",
    ),
    (
        "clp_all_variant_any",
        "COSMIC CLP v104 HCC1395 allele-containing region (all statuses)",
    ),
    (
        "clp_confirmed_somatic_variant_any",
        "COSMIC CLP v104 HCC1395 confirmed-somatic allele-containing region",
    ),
]


def parse_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes"}


def parse_region(region: str) -> tuple[str, int, int]:
    chrom, span = region.split(":", 1)
    start_text, end_text = span.split("-", 1)
    return chrom, int(start_text), int(end_text)


def wilson(success: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total == 0:
        return math.nan, math.nan
    p = success / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return center - half, center + half


def kappa(rows: list[dict]) -> float:
    if not rows:
        return math.nan
    categories = sorted({row["category_a"] for row in rows} | {row["category_b"] for row in rows})
    count_a = collections.Counter(row["category_a"] for row in rows)
    count_b = collections.Counter(row["category_b"] for row in rows)
    n = len(rows)
    observed = sum(row["agree"] for row in rows) / n
    expected = sum((count_a[value] / n) * (count_b[value] / n) for value in categories)
    return (observed - expected) / (1 - expected) if expected < 1 else math.nan


def length_deciles(rows: list[dict]) -> None:
    lengths = np.asarray([row["length_bp"] for row in rows], dtype=np.float64)
    cutpoints = np.unique(np.quantile(lengths, np.linspace(0.1, 0.9, 9)))
    for row in rows:
        row["length_decile"] = int(np.searchsorted(cutpoints, row["length_bp"], side="right")) + 1


def permutation_difference(rows: list[dict], field: str) -> dict:
    feature = np.asarray([bool(row[field]) for row in rows], dtype=bool)
    outcome = np.asarray([bool(row["agree"]) for row in rows], dtype=np.float64)
    n_present = int(feature.sum())
    n_absent = len(rows) - n_present
    if not n_present or not n_absent:
        return {
            "seed": SEED,
            "permutations": 0,
            "stratification": "chromosome + global region-length decile",
            "null_mean": math.nan,
            "null_q025": math.nan,
            "null_q975": math.nan,
            "p_two_sided": math.nan,
        }

    observed = float(outcome[feature].mean() - outcome[~feature].mean())
    # Conditional on each stratum's feature count, a label permutation is
    # exactly a hypergeometric draw of agreement successes into the present
    # group.  This avoids millions of small array permutations while preserving
    # the same chromosome + length-decile null distribution.
    groups = []
    for key in sorted({(row["chrom"], row["length_decile"]) for row in rows}):
        indices = np.asarray(
            [index for index, row in enumerate(rows) if (row["chrom"], row["length_decile"]) == key],
            dtype=np.int64,
        )
        n_group = int(len(indices))
        successes = int(outcome[indices].sum())
        present_in_group = int(feature[indices].sum())
        groups.append((successes, n_group - successes, present_in_group))
    rng = np.random.default_rng(SEED)
    present_success = np.zeros(N_PERMUTATIONS, dtype=np.int64)
    for good, bad, sample_size in groups:
        if sample_size:
            present_success += rng.hypergeometric(good, bad, sample_size, size=N_PERMUTATIONS)
    total_success = int(outcome.sum())
    null = present_success / n_present - (total_success - present_success) / n_absent
    null_mean = float(null.mean())
    return {
        "seed": SEED,
        "permutations": N_PERMUTATIONS,
        "stratification": "chromosome + global region-length decile",
        "null_mean": null_mean,
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
        "p_two_sided": float(
            (1 + np.count_nonzero(np.abs(null - null_mean) >= abs(observed - null_mean)))
            / (N_PERMUTATIONS + 1)
        ),
    }


def load_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        raw_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(raw_rows) != EXPECTED_ROWS:
        raise RuntimeError(f"expected {EXPECTED_ROWS} exact complete-both rows, observed {len(raw_rows)}")
    rows = []
    for raw in raw_rows:
        if raw.get("scenario") != "exact_coordinate":
            raise RuntimeError("annotation input contains non-exact scenario")
        if raw.get("region_a") != raw.get("region_b"):
            raise RuntimeError("exact-coordinate annotation row has non-identical region keys")
        chrom, start, end = parse_region(raw["region_a"])
        row = dict(raw)
        row.update(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "length_bp": end - start + 1,
                "agree": parse_bool(raw["category_agree"]),
            }
        )
        for field, _ in FEATURES:
            if field in raw:
                row[field] = parse_bool(raw[field])
            elif field == "clp_all_variant_any" and "clp_all_variant_count" in raw:
                row[field] = int(raw["clp_all_variant_count"] or 0) > 0
            elif field == "clp_confirmed_somatic_variant_any" and "clp_confirmed_somatic_variant_count" in raw:
                row[field] = int(raw["clp_confirmed_somatic_variant_count"] or 0) > 0
            else:
                raise RuntimeError(f"required annotation flag is missing: {field}")
        rows.append(row)
    if len({row["region_a"] for row in rows}) != len(rows):
        raise RuntimeError("exact annotation region keys are not unique")
    length_deciles(rows)
    return rows


def summarise(rows: list[dict]) -> tuple[list[dict], dict]:
    baseline_success = sum(row["agree"] for row in rows)
    baseline_low, baseline_high = wilson(baseline_success, len(rows))
    output = [
        {
            "feature": "ALL",
            "feature_label": "All exact complete-both regions",
            "n_present": len(rows),
            "n_absent": 0,
            "present_agree_n": baseline_success,
            "present_agreement": baseline_success / len(rows),
            "present_ci95_low": baseline_low,
            "present_ci95_high": baseline_high,
            "present_kappa": kappa(rows),
            "absent_agree_n": 0,
            "absent_agreement": math.nan,
            "absent_ci95_low": math.nan,
            "absent_ci95_high": math.nan,
            "absent_kappa": math.nan,
            "difference_present_minus_absent_pp": math.nan,
            "permutation_null_mean_pp": math.nan,
            "permutation_null_q025_pp": math.nan,
            "permutation_null_q975_pp": math.nan,
            "permutation_p_two_sided": math.nan,
        }
    ]
    permutation_meta = {}
    for field, label in FEATURES:
        present = [row for row in rows if row[field]]
        absent = [row for row in rows if not row[field]]
        present_success = sum(row["agree"] for row in present)
        absent_success = sum(row["agree"] for row in absent)
        present_low, present_high = wilson(present_success, len(present))
        absent_low, absent_high = wilson(absent_success, len(absent))
        present_rate = present_success / len(present) if present else math.nan
        absent_rate = absent_success / len(absent) if absent else math.nan
        null = permutation_difference(rows, field)
        permutation_meta[field] = null
        output.append(
            {
                "feature": field,
                "feature_label": label,
                "n_present": len(present),
                "n_absent": len(absent),
                "present_agree_n": present_success,
                "present_agreement": present_rate,
                "present_ci95_low": present_low,
                "present_ci95_high": present_high,
                "present_kappa": kappa(present),
                "absent_agree_n": absent_success,
                "absent_agreement": absent_rate,
                "absent_ci95_low": absent_low,
                "absent_ci95_high": absent_high,
                "absent_kappa": kappa(absent),
                "difference_present_minus_absent_pp": (present_rate - absent_rate) * 100,
                "permutation_null_mean_pp": null["null_mean"] * 100,
                "permutation_null_q025_pp": null["null_q025"] * 100,
                "permutation_null_q975_pp": null["null_q975"] * 100,
                "permutation_p_two_sided": null["p_two_sided"],
            }
        )
    metadata = {
        "schema_version": "1.0",
        "scope": "HCC1395 vs HCC1395_DORADO exact-coordinate complete-both regions",
        "grain": "one exact matched region",
        "n_regions": len(rows),
        "baseline_agreement": baseline_success / len(rows),
        "baseline_kappa": kappa(rows),
        "permutation": {
            "seed": SEED,
            "n_permutations": N_PERMUTATIONS,
            "stratification": "chromosome + global region-length decile",
            "test": "two-sided difference in category agreement; annotation is context, not truth",
        },
        "features": permutation_meta,
    }
    return output, metadata


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def json_safe(value: object) -> object:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    args = parser.parse_args()

    rows = load_rows(args.input)
    summary, metadata = summarise(rows)
    write_tsv(args.output_tsv, summary)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(json_safe({**metadata, "rows": summary}), ensure_ascii=False, indent=2, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )
    print(f"INPUT -> {args.input}")
    print(f"OUTPUT TSV -> {args.output_tsv}")
    print(f"OUTPUT JSON -> {args.output_json}")
    print(
        "RESULT -> "
        f"n={metadata['n_regions']}; agreement={metadata['baseline_agreement']:.6f}; "
        f"kappa={metadata['baseline_kappa']:.6f}; checks=PASS"
    )


if __name__ == "__main__":
    main()
