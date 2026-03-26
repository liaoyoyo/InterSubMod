#!/usr/bin/env python3
"""Build a Phase 1A head-to-head baseline summary table."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd

from build_phase1_training_manifest import dataset_role, harmonization_group
from research_common import ensure_dir, write_json, write_tsv_rows


DEFAULT_MANIFEST = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1_training_manifest_v1/phase1_training_manifest.tsv"
)

DEFAULT_READS = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1_manifest_shard_export_sample80/phase1_shard_read_training_table.tsv"
)

DEFAULT_OUTPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1a_head_to_head_baseline_v1"
)


BASELINE_FIELDS = [
    "dataset_id",
    "dataset_role",
    "harmonization_group",
    "truth_status",
    "regions_total",
    "regions_passed_gating",
    "regions_passed_gating_rate",
    "quality_score_mean",
    "quality_score_median",
    "pairwise_median_dist_mean",
    "pairwise_median_dist_median",
    "verification_noise_regions",
    "verification_weak_regions",
    "verification_strong_regions",
    "verification_subclone_regions",
    "sampled_reads_total",
    "sampled_alt_reads",
    "sampled_ref_reads",
    "sampled_alt_fraction",
    "sampled_num_cpg_observed_mean",
    "sampled_num_cpg_observed_median",
    "sampled_methyl_mean_mean",
    "sampled_methyl_mean_median",
    "sampled_methyl_std_mean",
    "sampled_methyl_std_median",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-tsv", default=str(DEFAULT_MANIFEST), help="Phase 1 training manifest TSV.")
    parser.add_argument("--reads-tsv", default=str(DEFAULT_READS), help="Phase 1A sampled shard read TSV.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT), help="Output directory.")
    return parser.parse_args()


def normalize_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1"})


def build_manifest_stats(manifest_df: pd.DataFrame) -> pd.DataFrame:
    working = manifest_df.copy()
    working["dataset_role"] = working["platform"].map(dataset_role)
    working["harmonization_group"] = [
        harmonization_group(str(platform), str(mode))
        for platform, mode in zip(working["platform"], working["mode"])
    ]
    working["PassedGatingNorm"] = normalize_bool(working["PassedGating"])

    grouped = working.groupby(["dataset_id", "dataset_role", "harmonization_group", "truth_status"], sort=True)
    rows: List[Dict[str, object]] = []
    for (dataset_id, dataset_role_value, group_name, truth_status), sub in grouped:
        verification_counts = sub["VerificationClass"].fillna("NA").value_counts()
        rows.append(
            {
                "dataset_id": dataset_id,
                "dataset_role": dataset_role_value,
                "harmonization_group": group_name,
                "truth_status": truth_status,
                "regions_total": int(len(sub.index)),
                "regions_passed_gating": int(sub["PassedGatingNorm"].sum()),
                "regions_passed_gating_rate": float(sub["PassedGatingNorm"].mean()),
                "quality_score_mean": float(sub["Quality_Score"].mean()),
                "quality_score_median": float(sub["Quality_Score"].median()),
                "pairwise_median_dist_mean": float(sub["PairwiseMedianDist"].mean()),
                "pairwise_median_dist_median": float(sub["PairwiseMedianDist"].median()),
                "verification_noise_regions": int(verification_counts.get("Noise", 0)),
                "verification_weak_regions": int(verification_counts.get("Weak", 0)),
                "verification_strong_regions": int(verification_counts.get("Strong", 0)),
                "verification_subclone_regions": int(verification_counts.get("Subclone", 0)),
            }
        )
    return pd.DataFrame(rows)


def build_read_stats(reads_df: pd.DataFrame) -> pd.DataFrame:
    working = reads_df.copy()
    for column in ["num_cpg_observed", "methyl_mean", "methyl_std"]:
        working[column] = pd.to_numeric(working[column], errors="coerce")

    grouped = working.groupby(["dataset_id", "dataset_role", "harmonization_group", "truth_status"], sort=True)
    rows: List[Dict[str, object]] = []
    for (dataset_id, dataset_role_value, group_name, truth_status), sub in grouped:
        label_counts = sub["phase1a_read_label"].fillna("UNKNOWN").value_counts()
        sampled_total = int(len(sub.index))
        sampled_alt = int(label_counts.get("ALT", 0))
        rows.append(
            {
                "dataset_id": dataset_id,
                "dataset_role": dataset_role_value,
                "harmonization_group": group_name,
                "truth_status": truth_status,
                "sampled_reads_total": sampled_total,
                "sampled_alt_reads": sampled_alt,
                "sampled_ref_reads": int(label_counts.get("REF", 0)),
                "sampled_alt_fraction": float(sampled_alt / sampled_total) if sampled_total else 0.0,
                "sampled_num_cpg_observed_mean": float(sub["num_cpg_observed"].mean()),
                "sampled_num_cpg_observed_median": float(sub["num_cpg_observed"].median()),
                "sampled_methyl_mean_mean": float(sub["methyl_mean"].mean()),
                "sampled_methyl_mean_median": float(sub["methyl_mean"].median()),
                "sampled_methyl_std_mean": float(sub["methyl_std"].mean()),
                "sampled_methyl_std_median": float(sub["methyl_std"].median()),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    output_dir = ensure_dir(Path(args.output_dir).resolve())
    manifest_df = pd.read_csv(args.manifest_tsv, sep="\t", low_memory=False)
    reads_df = pd.read_csv(args.reads_tsv, sep="\t", low_memory=False)

    manifest_stats = build_manifest_stats(manifest_df)
    read_stats = build_read_stats(reads_df)
    baseline_df = manifest_stats.merge(
        read_stats,
        on=["dataset_id", "dataset_role", "harmonization_group", "truth_status"],
        how="left",
    )
    baseline_rows = baseline_df.to_dict("records")

    summary_rows = [
        {
            "manifest_groups": int(len(manifest_stats.index)),
            "read_groups": int(len(read_stats.index)),
            "baseline_rows": int(len(baseline_df.index)),
            "manifest_tsv": args.manifest_tsv,
            "reads_tsv": args.reads_tsv,
        }
    ]

    write_tsv_rows(output_dir / "phase1a_head_to_head_baseline.tsv", BASELINE_FIELDS, baseline_rows)
    write_tsv_rows(
        output_dir / "phase1a_head_to_head_summary.tsv",
        ["manifest_groups", "read_groups", "baseline_rows", "manifest_tsv", "reads_tsv"],
        summary_rows,
    )
    write_json(
        output_dir / "run_context.json",
        {
            "task": "Phase 1A head-to-head baseline build",
            "manifest_tsv": args.manifest_tsv,
            "reads_tsv": args.reads_tsv,
            "output_dir": str(output_dir),
        },
    )

    notes = [
        "# Phase 1A Head-to-Head Baseline",
        "",
        f"- manifest_tsv: `{args.manifest_tsv}`",
        f"- reads_tsv: `{args.reads_tsv}`",
        "",
        "## Outputs",
        "",
        "- `phase1a_head_to_head_baseline.tsv`",
        "- `phase1a_head_to_head_summary.tsv`",
        "- `run_context.json`",
    ]
    (output_dir / "round_summary.md").write_text("\n".join(notes) + "\n", encoding="utf-8")

    print(f"[phase1a-baseline] wrote {output_dir / 'phase1a_head_to_head_baseline.tsv'}")
    print(f"[phase1a-baseline] wrote {output_dir / 'phase1a_head_to_head_summary.tsv'}")


if __name__ == "__main__":
    main()
