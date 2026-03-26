#!/usr/bin/env python3
"""Summarize Phase 1A prediction errors from benchmark prediction tables."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd

from research_common import ensure_dir, write_json, write_tsv_rows


DEFAULT_INPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1a_read_classifier_benchmark_sample200_v1/benchmark_predictions.tsv"
)

DEFAULT_OUTPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1a_to_pure_failure_diagnosis_v1"
)


SUMMARY_FIELDS = [
    "model_name",
    "evaluation_split",
    "harmonization_group",
    "truth_status",
    "VerificationClass",
    "PassedGating",
    "phase1a_read_label",
    "rows_total",
    "error_count",
    "error_rate",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions-tsv", default=str(DEFAULT_INPUT), help="Input benchmark_predictions.tsv path.")
    parser.add_argument("--model-name", default="logistic_methyl_context", help="Filter by model_name.")
    parser.add_argument("--evaluation-split", default="external_validation", help="Filter by evaluation split.")
    parser.add_argument("--harmonization-group", default="ONT_Dorado|to-pure", help="Filter by harmonization group.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT), help="Output directory.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = ensure_dir(Path(args.output_dir).resolve())
    df = pd.read_csv(args.predictions_tsv, sep="\t", low_memory=False)

    filtered = df[
        (df["model_name"] == args.model_name)
        & (df["evaluation_split"] == args.evaluation_split)
        & (df["harmonization_group"] == args.harmonization_group)
    ].copy()

    group_columns = [
        "model_name",
        "evaluation_split",
        "harmonization_group",
        "truth_status",
        "VerificationClass",
        "PassedGating",
        "phase1a_read_label",
    ]

    rows: List[Dict[str, object]] = []
    if not filtered.empty:
        grouped = filtered.groupby(group_columns, dropna=False, sort=True)
        for keys, sub in grouped:
            row = {column: value for column, value in zip(group_columns, keys)}
            error_count = int(sub["is_error"].astype(bool).sum())
            row.update(
                {
                    "rows_total": int(len(sub.index)),
                    "error_count": error_count,
                    "error_rate": float(error_count / len(sub.index)),
                }
            )
            rows.append(row)

    total_row = {
        "model_name": args.model_name,
        "evaluation_split": args.evaluation_split,
        "harmonization_group": args.harmonization_group,
        "rows_total": int(len(filtered.index)),
        "error_count": int(filtered["is_error"].astype(bool).sum()) if not filtered.empty else 0,
        "error_rate": float(filtered["is_error"].astype(bool).mean()) if not filtered.empty else 0.0,
    }

    write_tsv_rows(output_dir / "failure_breakdown.tsv", SUMMARY_FIELDS, rows)
    write_tsv_rows(
        output_dir / "failure_summary.tsv",
        ["model_name", "evaluation_split", "harmonization_group", "rows_total", "error_count", "error_rate"],
        [total_row],
    )
    write_json(
        output_dir / "run_context.json",
        {
            "task": "Phase 1A prediction failure summary",
            "predictions_tsv": args.predictions_tsv,
            "model_name": args.model_name,
            "evaluation_split": args.evaluation_split,
            "harmonization_group": args.harmonization_group,
            "output_dir": str(output_dir),
        },
    )

    notes = [
        "# Phase 1A Failure Diagnosis",
        "",
        f"- predictions_tsv: `{args.predictions_tsv}`",
        f"- model_name: `{args.model_name}`",
        f"- evaluation_split: `{args.evaluation_split}`",
        f"- harmonization_group: `{args.harmonization_group}`",
        "",
        "## Outputs",
        "",
        "- `failure_breakdown.tsv`",
        "- `failure_summary.tsv`",
        "- `run_context.json`",
    ]
    (output_dir / "round_summary.md").write_text("\n".join(notes) + "\n", encoding="utf-8")

    print(f"[phase1a-failure] wrote {output_dir / 'failure_breakdown.tsv'}")
    print(f"[phase1a-failure] wrote {output_dir / 'failure_summary.tsv'}")


if __name__ == "__main__":
    main()
