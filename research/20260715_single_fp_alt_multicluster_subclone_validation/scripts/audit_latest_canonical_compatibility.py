#!/usr/bin/env python3
"""Audit compatibility of latest and legacy canonical outputs at common FP sites."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from scipy.stats import binomtest


ENDPOINTS = [
    "stable_null_multigroup",
    "residual_unexplained_multigroup",
    "phase_anchored_robust_epigenetic_candidate",
]


def load_rows(path: Path) -> dict[tuple[str, str, int, str, str], dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {
        (row["sample"], row["chrom"], int(row["pos"]), row["ref"], row["alt"]): row for row in rows
    }


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def compare_file_pair(paths: tuple[Path, Path]) -> bool:
    first, second = paths
    if first.stat().st_size != second.stat().st_size:
        return False
    return digest(first) == digest(second)


def outcome(row: dict[str, str], field: str) -> bool:
    return row.get(field, "").lower() == "true"


def endpoint_comparison(
    keys: set[tuple[str, str, int, str, str]],
    latest: dict[tuple[str, str, int, str, str], dict[str, str]],
    canonical: dict[tuple[str, str, int, str, str], dict[str, str]],
    field: str,
) -> dict[str, Any]:
    latest_yes = {key for key in keys if outcome(latest[key], field)}
    canonical_yes = {key for key in keys if outcome(canonical[key], field)}
    both = len(latest_yes & canonical_yes)
    latest_only = len(latest_yes - canonical_yes)
    canonical_only = len(canonical_yes - latest_yes)
    union = len(latest_yes | canonical_yes)
    discordant = latest_only + canonical_only
    return {
        "n_common_sites": len(keys),
        "latest_n": len(latest_yes),
        "latest_fraction_all_common": len(latest_yes) / len(keys) if keys else None,
        "legacy_canonical_n": len(canonical_yes),
        "legacy_canonical_fraction_all_common": len(canonical_yes) / len(keys) if keys else None,
        "both_positive": both,
        "latest_only": latest_only,
        "legacy_canonical_only": canonical_only,
        "jaccard": both / union if union else None,
        "exact_mcnemar_p": (
            binomtest(latest_only, discordant, 0.5, alternative="two-sided").pvalue
            if discordant
            else 1.0
        ),
    }


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--latest-sites",
        type=Path,
        default=root / "results/focal_alt_multicluster/latest_full_v1/latest_site_results.tsv",
    )
    parser.add_argument(
        "--canonical-sites",
        type=Path,
        default=root / "results/focal_alt_multicluster/canonical_reference_v1/canonical_site_results.tsv",
    )
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument(
        "--summary", type=Path, default=root / "results/latest_canonical_compatibility_audit.json"
    )
    args = parser.parse_args()

    latest = load_rows(args.latest_sites)
    canonical = load_rows(args.canonical_sites)
    latest_keys = set(latest)
    canonical_keys = set(canonical)
    common = latest_keys & canonical_keys

    file_types = {
        "reads_tsv": Path("reads/reads.tsv"),
        "methylation_matrix": Path("methylation/methylation.csv"),
        "bernoulli_distance_matrix": Path("distance/BERNOULLI/matrix.csv"),
    }
    pair_jobs: list[tuple[str, str, tuple[Path, Path]]] = []
    for key in sorted(common):
        latest_root = Path(latest[key]["region_dir"])
        canonical_root = Path(canonical[key]["region_dir"])
        site = f"{key[0]}:{key[1]}:{key[2]}:{key[3]}>{key[4]}"
        for name, relative in file_types.items():
            pair_jobs.append((site, name, (latest_root / relative, canonical_root / relative)))

    identical_by_type: Counter[str] = Counter()
    different_examples: dict[str, list[str]] = defaultdict(list)
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        matches = executor.map(lambda job: compare_file_pair(job[2]), pair_jobs)
        for (site, name, _), identical in zip(pair_jobs, matches):
            if identical:
                identical_by_type[name] += 1
            elif len(different_examples[name]) < 10:
                different_examples[name].append(site)

    same_readset = sum(
        latest[key]["alt_readset_sha256"] == canonical[key]["alt_readset_sha256"] for key in common
    )
    per_sample: dict[str, Any] = {}
    for sample in sorted({key[0] for key in common}):
        sample_keys = {key for key in common if key[0] == sample}
        per_sample[sample] = {
            "n_common_sites": len(sample_keys),
            "n_same_alt_readset": sum(
                latest[key]["alt_readset_sha256"] == canonical[key]["alt_readset_sha256"]
                for key in sample_keys
            ),
            "endpoints": {
                field: endpoint_comparison(sample_keys, latest, canonical, field) for field in ENDPOINTS
            },
        }

    summary = {
        "schema_name": "intersubmod.latest_legacy_canonical_compatibility_audit",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {"latest_sites": str(args.latest_sites), "canonical_sites": str(args.canonical_sites)},
        "site_sets": {
            "latest": len(latest_keys),
            "legacy_canonical": len(canonical_keys),
            "common": len(common),
            "latest_only": len(latest_keys - canonical_keys),
            "legacy_canonical_only": len(canonical_keys - latest_keys),
        },
        "n_same_alt_readset": same_readset,
        "same_alt_readset_fraction_common": same_readset / len(common) if common else None,
        "file_identity": {
            name: {
                "n_identical": identical_by_type[name],
                "n_different": len(common) - identical_by_type[name],
                "fraction_identical": identical_by_type[name] / len(common) if common else None,
                "different_examples": different_examples[name],
            }
            for name in file_types
        },
        "endpoint_agreement": {
            field: endpoint_comparison(common, latest, canonical, field) for field in ENDPOINTS
        },
        "per_sample": per_sample,
        "interpretation_guardrail": (
            "Legacy canonical distance matrices were generated by an older InterSubMod command/binary without "
            "the current explicit min-common-coverage=3 and nan-distance-strategy=SKIP contract. Legacy endpoint "
            "rates are sensitivity context, not a directly exchangeable estimate for the latest analysis."
        ),
        "pass": len(latest) == 7745 and len(canonical) == 3458 and len(common) == 3302,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(args.summary), "site_sets": summary["site_sets"], "pass": summary["pass"]}, indent=2))


if __name__ == "__main__":
    main()
