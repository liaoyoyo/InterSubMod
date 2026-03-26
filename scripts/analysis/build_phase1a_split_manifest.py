#!/usr/bin/env python3
"""Build a Phase 1A region-level split manifest from the baseline Phase 1 training manifest."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd

from build_phase1_training_manifest import dataset_role, harmonization_group
from research_common import ensure_dir, write_json, write_tsv_rows


DEFAULT_INPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1_training_manifest_v1/phase1_training_manifest.tsv"
)

DEFAULT_OUTPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1a_split_manifest_v1"
)


SPLIT_FIELDS = [
    "dataset_id",
    "dataset_label",
    "sample",
    "platform",
    "mode",
    "dataset_role",
    "harmonization_group",
    "phase1a_task",
    "split_role",
    "region_key",
    "truth_status",
    "source_scope",
    "VerificationClass",
    "Quality_Score",
    "PairwiseMedianDist",
    "PassedGating",
    "region_resolve_status",
]


SUMMARY_FIELDS = [
    "dataset_id",
    "dataset_role",
    "split_role",
    "regions_total",
    "tp_regions",
    "fp_regions",
    "resolved_regions",
    "passed_gating_true",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-tsv", default=str(DEFAULT_INPUT), help="Input Phase 1 training manifest TSV.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT), help="Output directory.")
    return parser.parse_args()


def split_role_for_platform(platform: str) -> str:
    role = dataset_role(platform)
    if role == "discovery":
        return "discovery"
    if role == "validation":
        return "external_validation"
    return "unknown"


def build_split_rows(manifest_df: pd.DataFrame) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for _, row in manifest_df.iterrows():
        platform = str(row.get("platform", ""))
        mode = str(row.get("mode", ""))
        rows.append(
            {
                "dataset_id": row.get("dataset_id", ""),
                "dataset_label": row.get("dataset_label", ""),
                "sample": row.get("sample", ""),
                "platform": platform,
                "mode": mode,
                "dataset_role": dataset_role(platform),
                "harmonization_group": harmonization_group(platform, mode),
                "phase1a_task": "within_tumor_alt_support",
                "split_role": split_role_for_platform(platform),
                "region_key": row.get("region_key", ""),
                "truth_status": row.get("truth_status", ""),
                "source_scope": row.get("source_scope", ""),
                "VerificationClass": row.get("VerificationClass", ""),
                "Quality_Score": row.get("Quality_Score", ""),
                "PairwiseMedianDist": row.get("PairwiseMedianDist", ""),
                "PassedGating": row.get("PassedGating", ""),
                "region_resolve_status": row.get("region_resolve_status", ""),
            }
        )
    return rows


def build_summary_rows(split_df: pd.DataFrame) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    grouped = split_df.groupby(["dataset_id", "dataset_role", "split_role"], dropna=False, sort=True)
    for (dataset_id, dataset_role_value, split_role), sub in grouped:
        passed_gating = sub["PassedGating"].astype(str).str.lower().isin({"true", "1"})
        rows.append(
            {
                "dataset_id": dataset_id,
                "dataset_role": dataset_role_value,
                "split_role": split_role,
                "regions_total": int(len(sub.index)),
                "tp_regions": int((sub["truth_status"] == "TP").sum()),
                "fp_regions": int((sub["truth_status"] == "FP").sum()),
                "resolved_regions": int((sub["region_resolve_status"] == "resolved").sum()),
                "passed_gating_true": int(passed_gating.sum()),
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    output_dir = ensure_dir(Path(args.output_dir).resolve())
    manifest_df = pd.read_csv(args.manifest_tsv, sep="\t", low_memory=False)

    split_rows = build_split_rows(manifest_df)
    split_df = pd.DataFrame(split_rows)
    summary_rows = build_summary_rows(split_df)

    write_tsv_rows(output_dir / "phase1a_split_manifest.tsv", SPLIT_FIELDS, split_rows)
    write_tsv_rows(output_dir / "phase1a_split_summary.tsv", SUMMARY_FIELDS, summary_rows)
    write_json(
        output_dir / "run_context.json",
        {
            "task": "Phase 1A split manifest build",
            "manifest_tsv": args.manifest_tsv,
            "output_dir": str(output_dir),
        },
    )

    notes = [
        "# Phase 1A Split Manifest",
        "",
        f"- manifest_tsv: `{args.manifest_tsv}`",
        "",
        "## Outputs",
        "",
        "- `phase1a_split_manifest.tsv`",
        "- `phase1a_split_summary.tsv`",
        "- `run_context.json`",
    ]
    (output_dir / "round_summary.md").write_text("\n".join(notes) + "\n", encoding="utf-8")

    print(f"[phase1a-split] wrote {output_dir / 'phase1a_split_manifest.tsv'}")
    print(f"[phase1a-split] wrote {output_dir / 'phase1a_split_summary.tsv'}")


if __name__ == "__main__":
    main()
