#!/usr/bin/env python3
"""Summarize label/cluster class shifts stratified by TP/FP across round bundles."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List

import pandas as pd

from research_common import infer_platform, read_json, write_tsv_rows


SUMMARY_FIELDS = [
    "sample",
    "platform",
    "sample_dir",
    "class_shift",
    "source_scope",
    "count",
    "scope_total",
    "rate_within_scope",
    "tp_fraction_within_shift",
    "fp_fraction_within_shift",
    "median_effect_size",
    "median_fisher_p",
    "median_cluster_permanova_p",
    "median_label_hp_permanova_p",
    "median_label_allele_permanova_p",
]

TOP_FIELDS = [
    "sample",
    "platform",
    "sample_dir",
    "class_shift",
    "source_scope",
    "region_key",
    "cluster_class",
    "label_class",
    "agreement_type",
    "effect_size",
    "fisher_p",
    "cluster_permanova_p",
    "label_hp_permanova_p",
    "label_allele_permanova_p",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", action="append", required=True, help="Round sample bundle directory")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument(
        "--shift",
        action="append",
        default=["Weak->Strong", "Noise->Strong", "Weak->Subclone", "Noise->Subclone"],
        help="Class shift to summarize",
    )
    parser.add_argument("--top-n", type=int, default=20, help="Top regions per shift/scope by effect size")
    return parser.parse_args()


def median_or_blank(series: pd.Series) -> str:
    clean = pd.to_numeric(series, errors="coerce").dropna()
    if clean.empty:
        return ""
    return f"{clean.median():.6f}"


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    shifts = list(dict.fromkeys(args.shift))
    summary_rows: List[Dict[str, object]] = []
    top_rows: List[Dict[str, object]] = []

    for sample_dir_str in args.sample_dir:
        sample_dir = Path(sample_dir_str).resolve()
        agreement_path = sample_dir / "label_cluster_agreement.tsv"
        if not agreement_path.exists():
            continue

        context = read_json(sample_dir / "round_context.json") if (sample_dir / "round_context.json").exists() else {}
        sample = str(context.get("sample") or sample_dir.name)
        platform = str(context.get("platform") or infer_platform(sample))
        df = pd.read_csv(agreement_path, sep="\t")
        if df.empty:
            continue

        scope_totals = df["source_scope"].value_counts().to_dict()

        for shift in shifts:
            shift_df = df[df["class_shift"] == shift].copy()
            if shift_df.empty:
                for scope in ("tp", "fp"):
                    summary_rows.append(
                        {
                            "sample": sample,
                            "platform": platform,
                            "sample_dir": str(sample_dir),
                            "class_shift": shift,
                            "source_scope": scope,
                            "count": 0,
                            "scope_total": int(scope_totals.get(scope, 0)),
                            "rate_within_scope": "0.000000",
                            "tp_fraction_within_shift": "",
                            "fp_fraction_within_shift": "",
                            "median_effect_size": "",
                            "median_fisher_p": "",
                            "median_cluster_permanova_p": "",
                            "median_label_hp_permanova_p": "",
                            "median_label_allele_permanova_p": "",
                        }
                    )
                continue

            shift_scope_counts = shift_df["source_scope"].value_counts().to_dict()
            tp_shift_total = int(shift_scope_counts.get("tp", 0))
            fp_shift_total = int(shift_scope_counts.get("fp", 0))
            shift_total = max(tp_shift_total + fp_shift_total, 1)

            for scope in ("tp", "fp"):
                scope_df = shift_df[shift_df["source_scope"] == scope].copy()
                scope_count = len(scope_df.index)
                scope_total = int(scope_totals.get(scope, 0))
                summary_rows.append(
                    {
                        "sample": sample,
                        "platform": platform,
                        "sample_dir": str(sample_dir),
                        "class_shift": shift,
                        "source_scope": scope,
                        "count": scope_count,
                        "scope_total": scope_total,
                        "rate_within_scope": f"{(scope_count / scope_total):.6f}" if scope_total else "0.000000",
                        "tp_fraction_within_shift": f"{(tp_shift_total / shift_total):.6f}",
                        "fp_fraction_within_shift": f"{(fp_shift_total / shift_total):.6f}",
                        "median_effect_size": median_or_blank(scope_df["effect_size"]) if scope_count else "",
                        "median_fisher_p": median_or_blank(scope_df["fisher_p"]) if scope_count else "",
                        "median_cluster_permanova_p": median_or_blank(scope_df["cluster_permanova_p"]) if scope_count else "",
                        "median_label_hp_permanova_p": median_or_blank(scope_df["label_hp_permanova_p"]) if scope_count else "",
                        "median_label_allele_permanova_p": median_or_blank(scope_df["label_allele_permanova_p"]) if scope_count else "",
                    }
                )

                if scope_count == 0:
                    continue
                ranked = scope_df.copy()
                ranked["effect_size_numeric"] = pd.to_numeric(ranked["effect_size"], errors="coerce")
                ranked["fisher_p_numeric"] = pd.to_numeric(ranked["fisher_p"], errors="coerce")
                ranked = ranked.sort_values(
                    by=["effect_size_numeric", "fisher_p_numeric"],
                    ascending=[False, True],
                    na_position="last",
                ).head(args.top_n)
                for _, row in ranked.iterrows():
                    top_rows.append(
                        {
                            "sample": sample,
                            "platform": platform,
                            "sample_dir": str(sample_dir),
                            "class_shift": shift,
                            "source_scope": scope,
                            "region_key": row.get("region_key", ""),
                            "cluster_class": row.get("cluster_class", ""),
                            "label_class": row.get("label_class", ""),
                            "agreement_type": row.get("agreement_type", ""),
                            "effect_size": row.get("effect_size", ""),
                            "fisher_p": row.get("fisher_p", ""),
                            "cluster_permanova_p": row.get("cluster_permanova_p", ""),
                            "label_hp_permanova_p": row.get("label_hp_permanova_p", ""),
                            "label_allele_permanova_p": row.get("label_allele_permanova_p", ""),
                        }
                    )

    write_tsv_rows(output_dir / "shift_tp_fp_summary.tsv", SUMMARY_FIELDS, summary_rows)
    write_tsv_rows(output_dir / "shift_top_regions.tsv", TOP_FIELDS, top_rows)

    md_lines = [
        "# Shift Pattern Summary",
        "",
        "- 目的：觀察 `Weak->Strong / Noise->Strong` 等升級型 shift 是否偏向 TP 或 FP。",
        "- 主表：`shift_tp_fp_summary.tsv`",
        "- Top regions：`shift_top_regions.tsv`",
        "",
    ]
    grouped: Dict[str, List[Dict[str, object]]] = {}
    for row in summary_rows:
        grouped.setdefault(str(row["sample"]), []).append(row)
    for sample, rows in grouped.items():
        md_lines.append(f"## {sample}")
        md_lines.append("")
        md_lines.append("| class_shift | source_scope | count | scope_total | rate_within_scope | tp_fraction_within_shift | fp_fraction_within_shift |")
        md_lines.append("| --- | --- | --- | --- | --- | --- | --- |")
        for row in rows:
            md_lines.append(
                f"| {row['class_shift']} | {row['source_scope']} | {row['count']} | {row['scope_total']} | {row['rate_within_scope']} | {row['tp_fraction_within_shift']} | {row['fp_fraction_within_shift']} |"
            )
        md_lines.append("")
    (output_dir / "shift_pattern_summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"[analyze_shift_patterns] Wrote {output_dir / 'shift_tp_fp_summary.tsv'}")
    print(f"[analyze_shift_patterns] Wrote {output_dir / 'shift_top_regions.tsv'}")


if __name__ == "__main__":
    main()
