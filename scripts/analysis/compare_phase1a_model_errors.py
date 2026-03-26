#!/usr/bin/env python3
"""Compare Phase 1A model error rates across datasets and verification buckets."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd

from research_common import ensure_dir, write_json, write_tsv_rows


DEFAULT_INPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1a_read_classifier_benchmark_paired_multibio_sample637_v1/benchmark_predictions.tsv"
)

DEFAULT_OUTPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1a_model_error_compare_v1"
)


PAIR_KEYS = [
    "evaluation_split",
    "dataset_id",
    "dataset_label",
    "dataset_role",
    "harmonization_group",
    "region_key",
    "read_id",
    "phase1a_read_label",
    "truth_status",
    "VerificationClass",
    "PassedGating",
]


BUCKET_FIELDS = [
    "evaluation_split",
    "dataset_id",
    "dataset_label",
    "harmonization_group",
    "truth_status",
    "VerificationClass",
    "PassedGating",
    "rows_total",
    "model_a",
    "model_a_error_count",
    "model_a_error_rate",
    "model_b",
    "model_b_error_count",
    "model_b_error_rate",
    "delta_error_count",
    "delta_error_rate",
    "improved_rows",
    "worsened_rows",
]


DATASET_FIELDS = [
    "evaluation_split",
    "dataset_id",
    "dataset_label",
    "harmonization_group",
    "rows_total",
    "model_a",
    "model_a_error_count",
    "model_a_error_rate",
    "model_b",
    "model_b_error_count",
    "model_b_error_rate",
    "delta_error_count",
    "delta_error_rate",
    "improved_rows",
    "worsened_rows",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions-tsv", default=str(DEFAULT_INPUT), help="Input benchmark_predictions.tsv path.")
    parser.add_argument("--evaluation-split", default="external_validation", help="Filter by evaluation split.")
    parser.add_argument("--model-a", default="logistic_context_only", help="Baseline model name.")
    parser.add_argument("--model-b", default="logistic_methyl_context", help="Comparison model name.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT), help="Output directory.")
    return parser.parse_args()


def build_pairwise_table(df: pd.DataFrame, model_a: str, model_b: str) -> pd.DataFrame:
    filtered = df[
        (df["model_name"].isin([model_a, model_b]))
    ].copy()

    pivot = filtered.pivot_table(
        index=PAIR_KEYS,
        columns="model_name",
        values="is_error",
        aggfunc="first",
    ).reset_index()
    if model_a not in pivot.columns or model_b not in pivot.columns:
        return pivot.iloc[0:0].copy()

    pivot["model_a_error"] = pivot[model_a].astype(bool)
    pivot["model_b_error"] = pivot[model_b].astype(bool)
    pivot["improved"] = pivot["model_a_error"] & (~pivot["model_b_error"])
    pivot["worsened"] = (~pivot["model_a_error"]) & pivot["model_b_error"]
    return pivot


def summarize_group(
    sub: pd.DataFrame,
    model_a: str,
    model_b: str,
) -> Dict[str, object]:
    rows_total = int(len(sub.index))
    model_a_error_count = int(sub["model_a_error"].sum())
    model_b_error_count = int(sub["model_b_error"].sum())
    return {
        "rows_total": rows_total,
        "model_a": model_a,
        "model_a_error_count": model_a_error_count,
        "model_a_error_rate": float(model_a_error_count / rows_total) if rows_total else 0.0,
        "model_b": model_b,
        "model_b_error_count": model_b_error_count,
        "model_b_error_rate": float(model_b_error_count / rows_total) if rows_total else 0.0,
        "delta_error_count": model_b_error_count - model_a_error_count,
        "delta_error_rate": float((model_b_error_count - model_a_error_count) / rows_total) if rows_total else 0.0,
        "improved_rows": int(sub["improved"].sum()),
        "worsened_rows": int(sub["worsened"].sum()),
    }


def main() -> None:
    args = parse_args()
    output_dir = ensure_dir(Path(args.output_dir).resolve())
    df = pd.read_csv(args.predictions_tsv, sep="\t", low_memory=False)
    df = df[df["evaluation_split"] == args.evaluation_split].copy()

    pairwise = build_pairwise_table(df, args.model_a, args.model_b)

    bucket_rows: List[Dict[str, object]] = []
    if not pairwise.empty:
        grouped = pairwise.groupby(
            [
                "evaluation_split",
                "dataset_id",
                "dataset_label",
                "harmonization_group",
                "truth_status",
                "VerificationClass",
                "PassedGating",
            ],
            sort=True,
            dropna=False,
        )
        for keys, sub in grouped:
            row = {field: value for field, value in zip(
                [
                    "evaluation_split",
                    "dataset_id",
                    "dataset_label",
                    "harmonization_group",
                    "truth_status",
                    "VerificationClass",
                    "PassedGating",
                ],
                keys,
            )}
            row.update(summarize_group(sub, args.model_a, args.model_b))
            bucket_rows.append(row)

    dataset_rows: List[Dict[str, object]] = []
    if not pairwise.empty:
        grouped = pairwise.groupby(
            ["evaluation_split", "dataset_id", "dataset_label", "harmonization_group"],
            sort=True,
            dropna=False,
        )
        for keys, sub in grouped:
            row = {field: value for field, value in zip(
                ["evaluation_split", "dataset_id", "dataset_label", "harmonization_group"],
                keys,
            )}
            row.update(summarize_group(sub, args.model_a, args.model_b))
            dataset_rows.append(row)

    write_tsv_rows(output_dir / "model_error_bucket_summary.tsv", BUCKET_FIELDS, bucket_rows)
    write_tsv_rows(output_dir / "model_error_dataset_summary.tsv", DATASET_FIELDS, dataset_rows)
    write_json(
        output_dir / "run_context.json",
        {
            "task": "Phase 1A model error comparison",
            "predictions_tsv": args.predictions_tsv,
            "evaluation_split": args.evaluation_split,
            "model_a": args.model_a,
            "model_b": args.model_b,
            "output_dir": str(output_dir),
        },
    )

    notes = [
        "# Phase 1A Model Error Comparison",
        "",
        f"- predictions_tsv: `{args.predictions_tsv}`",
        f"- evaluation_split: `{args.evaluation_split}`",
        f"- model_a: `{args.model_a}`",
        f"- model_b: `{args.model_b}`",
        "",
        "## Outputs",
        "",
        "- `model_error_bucket_summary.tsv`",
        "- `model_error_dataset_summary.tsv`",
        "- `run_context.json`",
    ]
    (output_dir / "round_summary.md").write_text("\n".join(notes) + "\n", encoding="utf-8")

    print(f"[phase1a-compare] wrote {output_dir / 'model_error_bucket_summary.tsv'}")
    print(f"[phase1a-compare] wrote {output_dir / 'model_error_dataset_summary.tsv'}")


if __name__ == "__main__":
    main()
