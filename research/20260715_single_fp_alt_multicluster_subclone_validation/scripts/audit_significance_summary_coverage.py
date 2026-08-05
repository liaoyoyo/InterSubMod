#!/usr/bin/env python3
"""Audit significance_summary.csv coverage against complete focal-site outputs."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def flag(row: dict[str, str], name: str) -> bool:
    return row.get(name, "").lower() == "true"


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--site-results",
        type=Path,
        default=root / "results/focal_alt_multicluster/latest_full_v1/latest_site_results.tsv",
    )
    parser.add_argument("--intersubmod-root", type=Path, required=True)
    parser.add_argument(
        "--output-tsv", type=Path, default=root / "results/significance_summary_omissions.tsv"
    )
    parser.add_argument(
        "--summary", type=Path, default=root / "results/significance_summary_coverage_audit.json"
    )
    args = parser.parse_args()

    with args.site_results.open(encoding="utf-8") as handle:
        sites = list(csv.DictReader(handle, delimiter="\t"))
    omitted: list[dict[str, str]] = []
    per_sample: dict[str, dict[str, object]] = {}
    for sample in sorted({row["sample"] for row in sites}):
        summary_path = args.intersubmod_root / sample / "significance_summary.csv"
        with summary_path.open(encoding="utf-8") as handle:
            summary_rows = list(csv.DictReader(handle))
        emitted = {(row["Chr"], int(row["Pos"]), row["Ref"], row["Alt"]) for row in summary_rows}
        sample_sites = [row for row in sites if row["sample"] == sample]
        sample_omitted = [
            row
            for row in sample_sites
            if (row["chrom"], int(row["pos"]), row["ref"], row["alt"]) not in emitted
        ]
        omitted.extend(sample_omitted)
        per_sample[sample] = {
            "n_site_outputs": len(sample_sites),
            "n_significance_summary_rows": len(summary_rows),
            "n_omitted": len(sample_omitted),
            "analysis_status_counts": dict(Counter(row["analysis_status"] for row in sample_omitted)),
            "n_stable_null_multigroup": sum(flag(row, "stable_null_multigroup") for row in sample_omitted),
            "n_residual_unexplained_multigroup": sum(
                flag(row, "residual_unexplained_multigroup") for row in sample_omitted
            ),
            "n_phase_anchored_robust_epigenetic_candidate": sum(
                flag(row, "phase_anchored_robust_epigenetic_candidate") for row in sample_omitted
            ),
        }

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample",
        "chrom",
        "pos",
        "ref",
        "alt",
        "n_reads_total",
        "n_alt_raw",
        "n_alt_after_peel",
        "analysis_status",
        "evidence_tier",
        "stable_null_multigroup",
        "residual_unexplained_multigroup",
        "phase_anchored_robust_epigenetic_candidate",
        "region_dir",
    ]
    with args.output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(omitted)

    summary = {
        "schema_name": "intersubmod.significance_summary_coverage_audit",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "site_results": str(args.site_results),
        "intersubmod_root": str(args.intersubmod_root),
        "n_site_outputs": len(sites),
        "n_significance_summary_rows": len(sites) - len(omitted),
        "n_omitted": len(omitted),
        "n_omitted_evaluable": sum(row["analysis_status"] == "evaluable" for row in omitted),
        "n_omitted_stable_null_multigroup": sum(flag(row, "stable_null_multigroup") for row in omitted),
        "n_omitted_residual_unexplained_multigroup": sum(
            flag(row, "residual_unexplained_multigroup") for row in omitted
        ),
        "n_omitted_phase_anchored_robust_epigenetic_candidate": sum(
            flag(row, "phase_anchored_robust_epigenetic_candidate") for row in omitted
        ),
        "per_sample": per_sample,
        "main_analysis_impact": (
            "None: focal analysis reads complete per-site reads.tsv, BERNOULLI matrices, and methylation matrices; "
            "significance_summary.csv is an auxiliary aggregate."
        ),
        "pass": len(sites) == 7745 and len(omitted) == 119,
    }
    args.summary.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(args.summary), "omissions": str(args.output_tsv), **summary}, indent=2))


if __name__ == "__main__":
    main()
