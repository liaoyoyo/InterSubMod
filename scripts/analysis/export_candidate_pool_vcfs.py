#!/usr/bin/env python3
"""Export candidate-eligible lost-TP / removed-FP VCFs from a rescue candidate pool."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Dict, Iterable, Set, TextIO

import pandas as pd

from research_common import ensure_dir, write_tsv_rows


SUMMARY_FIELDS = [
    "sample",
    "caller_vcf",
    "candidate_pool_tsv",
    "downstream_status",
    "target_records",
    "exported_records",
    "missing_records",
    "output_vcf",
    "missing_tsv",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-pool-tsv", required=True, help="TSV from extract_borderline_rescue_candidates.py")
    parser.add_argument("--caller-vcf", required=True, help="Caller final VCF(.gz)")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--sample", default="", help="Optional sample label")
    parser.add_argument(
        "--force-pass-filter",
        action="store_true",
        help="Rewrite exported candidate records to FILTER=PASS for downstream analysis-only VCF inputs",
    )
    return parser.parse_args()


def open_vcf(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def build_targets(df: pd.DataFrame, downstream_status: str) -> Set[str]:
    subset = df[(df["downstream_status"] == downstream_status) & (df["candidate_eligible"])].copy()
    return set(subset["region_key"].astype(str))


def write_missing(path: Path, keys: Iterable[str]) -> None:
    rows = [{"region_key": key} for key in sorted(keys)]
    write_tsv_rows(path, ["region_key"], rows)


def export_vcf_subset(
    caller_vcf: Path,
    target_keys: Set[str],
    output_vcf: Path,
    force_pass_filter: bool,
) -> Set[str]:
    exported: Set[str] = set()
    ensure_dir(output_vcf.parent)
    with open_vcf(caller_vcf) as src, output_vcf.open("w", encoding="utf-8") as dst:
        for raw in src:
            if raw.startswith("#"):
                dst.write(raw)
                continue
            fields = raw.rstrip("\n").split("\t")
            key = f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}"
            if key in target_keys:
                if force_pass_filter and len(fields) >= 7:
                    fields[6] = "PASS"
                    dst.write("\t".join(fields) + "\n")
                else:
                    dst.write(raw)
                exported.add(key)
    return exported


def main() -> None:
    args = parse_args()
    candidate_pool_tsv = Path(args.candidate_pool_tsv).resolve()
    caller_vcf = Path(args.caller_vcf).resolve()
    output_dir = ensure_dir(Path(args.output_dir).resolve())

    df = pd.read_csv(candidate_pool_tsv, sep="\t")
    required_columns = {"region_key", "downstream_status", "candidate_eligible"}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(f"candidate pool missing columns: {sorted(missing_columns)}")
    df["candidate_eligible"] = df["candidate_eligible"].astype(bool)

    spec: Dict[str, str] = {
        "caller_lost_tp": "candidate_lost_tp.vcf",
        "caller_removed_fp": "candidate_removed_fp.vcf",
    }

    summary_rows = []
    for downstream_status, filename in spec.items():
        target_keys = build_targets(df, downstream_status)
        output_vcf = output_dir / filename
        exported_keys = export_vcf_subset(caller_vcf, target_keys, output_vcf, args.force_pass_filter)
        missing_keys = target_keys.difference(exported_keys)
        missing_tsv = output_dir / f"{output_vcf.stem}_missing.tsv"
        write_missing(missing_tsv, missing_keys)

        summary_rows.append(
            {
                "sample": args.sample,
                "caller_vcf": str(caller_vcf),
                "candidate_pool_tsv": str(candidate_pool_tsv),
                "downstream_status": downstream_status,
                "target_records": len(target_keys),
                "exported_records": len(exported_keys),
                "missing_records": len(missing_keys),
                "output_vcf": str(output_vcf),
                "missing_tsv": str(missing_tsv),
            }
        )

    summary_tsv = output_dir / "candidate_vcf_export_summary.tsv"
    write_tsv_rows(summary_tsv, SUMMARY_FIELDS, summary_rows)

    lines = [
        "# Candidate VCF Export Summary",
        "",
        f"- sample: `{args.sample}`" if args.sample else "- sample: `<not set>`",
        f"- caller_vcf: `{caller_vcf}`",
        f"- candidate_pool_tsv: `{candidate_pool_tsv}`",
        f"- force_pass_filter: `{args.force_pass_filter}`",
        "",
    ]
    for row in summary_rows:
        lines.extend(
            [
                f"## {row['downstream_status']}",
                "",
                f"- target_records: `{row['target_records']}`",
                f"- exported_records: `{row['exported_records']}`",
                f"- missing_records: `{row['missing_records']}`",
                f"- output_vcf: `{row['output_vcf']}`",
                f"- missing_tsv: `{row['missing_tsv']}`",
                "",
            ]
        )

    (output_dir / "candidate_vcf_export_summary.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
