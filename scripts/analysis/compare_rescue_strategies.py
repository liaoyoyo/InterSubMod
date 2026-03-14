#!/usr/bin/env python3
"""Compare rescue strategies across one or more rescue_rule_comparison.tsv files."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import pandas as pd

from research_common import write_tsv_rows


COMPARISON_FIELDS = [
    "sample",
    "platform",
    "caller",
    "mode",
    "strategy",
    "rule_id",
    "tp_rescued",
    "fp_reintroduced",
    "fp_per_tp",
    "f1",
    "delta_f1_vs_baseline",
    "meets_safety",
    "rank_within_strategy",
    "notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rule-tsv", action="append", required=True, help="Input rescue_rule_comparison.tsv (repeatable)")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    return parser.parse_args()


def ranking_key(row: pd.Series):
    meets_safety = bool(row["meets_safety"])
    tp_rescued = int(row["tp_rescued"])
    fp_per_tp = float(row["fp_per_tp"]) if str(row["fp_per_tp"]) != "inf" else float("inf")
    delta_f1 = float(row["delta_f1_vs_baseline"])
    return (
        not meets_safety,
        -tp_rescued,
        fp_per_tp,
        -delta_f1,
    )


def to_output_rows(df: pd.DataFrame) -> List[dict]:
    rows: List[dict] = []
    for record in df.to_dict("records"):
        rows.append({field: record.get(field, "") for field in COMPARISON_FIELDS})
    return rows


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    frames = [pd.read_csv(path, sep="\t") for path in args.rule_tsv]
    combined = pd.concat(frames, ignore_index=True)

    for column in ("tp_rescued", "fp_reintroduced"):
        combined[column] = pd.to_numeric(combined[column], errors="coerce").fillna(0).astype(int)
    for column in ("fp_per_tp", "f1", "delta_f1_vs_baseline"):
        combined[column] = pd.to_numeric(combined[column], errors="coerce")
    combined["meets_safety"] = combined["meets_safety"].fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})

    comparison_rows: List[dict] = []
    best_rows: List[dict] = []

    grouped = combined.groupby(["sample", "platform", "caller", "mode", "strategy"], dropna=False)
    for _group_key, group_df in grouped:
        sorted_group = group_df.sort_values(
            by=["meets_safety", "tp_rescued", "fp_per_tp", "delta_f1_vs_baseline"],
            ascending=[False, False, True, False],
            na_position="last",
        ).copy()
        sorted_group["rank_within_strategy"] = range(1, len(sorted_group) + 1)
        comparison_rows.extend(to_output_rows(sorted_group))
        best_rows.append({field: sorted_group.iloc[0].get(field, "") for field in COMPARISON_FIELDS})

    comparison_df = pd.DataFrame(comparison_rows)
    best_df = pd.DataFrame(best_rows).sort_values(
        by=["sample", "mode", "strategy", "rank_within_strategy"], ascending=[True, True, True, True]
    )

    write_tsv_rows(output_dir / "strategy_comparison.tsv", COMPARISON_FIELDS, comparison_rows)
    write_tsv_rows(output_dir / "strategy_best_by_mode.tsv", COMPARISON_FIELDS, best_rows)

    md_lines = [
        "# Rescue Strategy Comparison",
        "",
        "## Best rule per sample / mode / strategy",
        "",
        "| sample | mode | strategy | rule_id | tp_rescued | fp_reintroduced | fp_per_tp | delta_f1_vs_baseline | meets_safety |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in best_df.to_dict("records"):
        md_lines.append(
            f"| {row['sample']} | {row['mode']} | {row['strategy']} | {row['rule_id']} | {row['tp_rescued']} | {row['fp_reintroduced']} | {row['fp_per_tp']} | {row['delta_f1_vs_baseline']} | {row['meets_safety']} |"
        )
    (output_dir / "strategy_summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"[compare_rescue_strategies] Wrote {output_dir / 'strategy_comparison.tsv'}")
    print(f"[compare_rescue_strategies] Wrote {output_dir / 'strategy_best_by_mode.tsv'}")


if __name__ == "__main__":
    main()
