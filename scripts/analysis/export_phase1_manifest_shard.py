#!/usr/bin/env python3
"""Export read-level shards from a Phase 1 training manifest."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List

import pandas as pd

from build_phase1_training_manifest import (
    DATASET_CONFIG,
    READ_LEVEL_FIELDS,
    load_methylation,
    load_reads,
    phase1_read_context,
    resolve_region_dir,
    summarize_methylation_row,
)
from research_common import ensure_dir, write_json, write_tsv_rows

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (
    VERIFICATION_PROVENANCE_COLUMNS,
    extract_provenance_frame,
)


OUTPUT_ROOT_DEFAULT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1_manifest_shard_export_v1"
)


SHARD_MANIFEST_FIELDS = list(dict.fromkeys([
    "dataset_id",
    "dataset_label",
    "region_key",
    "truth_status",
    "source_scope",
    "Quality_Score",
    *VERIFICATION_PROVENANCE_COLUMNS,
    "region_resolve_status",
    "region_dir",
    "summary_num_reads",
    "summary_num_cpgs",
    "read_rows_exported",
]))


READ_OUTPUT_FIELDS = list(dict.fromkeys([
    *READ_LEVEL_FIELDS,
    *VERIFICATION_PROVENANCE_COLUMNS,
]))


SHARD_SUMMARY_FIELDS = [
    "dataset_id",
    "source_scope",
    "selected_regions",
    "resolved_regions",
    "missing_regions",
    "read_rows_exported",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-tsv", required=True, help="Input Phase 1 manifest TSV.")
    parser.add_argument(
        "--dataset",
        action="append",
        choices=sorted(DATASET_CONFIG.keys()),
        help="Restrict to dataset id. Default: all datasets in manifest.",
    )
    parser.add_argument(
        "--source-scope",
        action="append",
        choices=["tp", "fp"],
        help="Restrict to source scope. Default: tp + fp.",
    )
    parser.add_argument(
        "--max-regions-per-dataset-scope",
        type=int,
        default=1,
        help="How many regions to export per dataset per scope.",
    )
    parser.add_argument("--output-dir", default=str(OUTPUT_ROOT_DEFAULT), help="Output directory.")
    return parser.parse_args()


def attach_verification_provenance(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate complete v2 provenance and retain it for every downstream row."""
    provenance = extract_provenance_frame(frame)
    prepared = frame.copy()
    for column in VERIFICATION_PROVENANCE_COLUMNS:
        prepared[column] = provenance[column]
    prepared.attrs["verification_provenance"] = {
        "schema_status": provenance.attrs["schema_status"],
        "columns": list(VERIFICATION_PROVENANCE_COLUMNS),
        "row_count": len(prepared),
    }
    return prepared


def provenance_payload(row: pd.Series) -> Dict[str, object]:
    return {column: row.get(column, "") for column in VERIFICATION_PROVENANCE_COLUMNS}


def select_manifest_rows(
    manifest_df: pd.DataFrame,
    dataset_ids: List[str],
    scopes: List[str],
    limit: int,
) -> pd.DataFrame:
    df = manifest_df.copy()
    if dataset_ids:
        df = df[df["dataset_id"].isin(dataset_ids)].copy()
    if scopes:
        df = df[df["source_scope"].isin(scopes)].copy()

    chosen: List[pd.DataFrame] = []
    for dataset_id in df["dataset_id"].drop_duplicates().tolist():
        for scope in scopes:
            sub = df[(df["dataset_id"] == dataset_id) & (df["source_scope"] == scope)].copy()
            if sub.empty:
                continue
            sub = sub.sort_values(
                ["Quality_Score", "summary_num_reads", "PairwiseMedianDist"],
                ascending=[False, False, False],
            )
            chosen.append(sub.head(limit))
    if not chosen:
        return df.iloc[0:0].copy()
    return pd.concat(chosen, ignore_index=True)


def manifest_region_dir(manifest_row: pd.Series) -> Path | None:
    resolve_status = str(manifest_row.get("region_resolve_status", "") or "")
    region_dir_value = manifest_row.get("region_dir", "")
    if resolve_status != "resolved" or pd.isna(region_dir_value):
        return None

    region_dir_text = str(region_dir_value).strip()
    if not region_dir_text or region_dir_text.lower() == "nan":
        return None

    region_dir = Path(region_dir_text)
    if not region_dir.exists():
        return None
    return region_dir


def export_shard_rows(selected_df: pd.DataFrame) -> tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    read_rows: List[Dict[str, object]] = []
    shard_manifest_rows: List[Dict[str, object]] = []

    for _, manifest_row in selected_df.iterrows():
        dataset_id = str(manifest_row["dataset_id"])
        cfg = DATASET_CONFIG[dataset_id]
        scope = str(manifest_row["source_scope"])
        region_key = str(manifest_row["region_key"])
        search_root = cfg.tp_region_root if scope == "tp" else cfg.fp_region_root

        region_dir = manifest_region_dir(manifest_row)
        if region_dir is None:
            region_dir = resolve_region_dir(search_root, region_key)
        if region_dir is None or not region_dir.exists():
            shard_manifest_rows.append(
                {
                    "dataset_id": dataset_id,
                    "dataset_label": manifest_row["dataset_label"],
                    "region_key": region_key,
                    "truth_status": manifest_row["truth_status"],
                    "source_scope": scope,
                    "Quality_Score": manifest_row["Quality_Score"],
                    **provenance_payload(manifest_row),
                    "region_resolve_status": "missing",
                    "region_dir": "",
                    "summary_num_reads": manifest_row["summary_num_reads"],
                    "summary_num_cpgs": manifest_row["summary_num_cpgs"],
                    "read_rows_exported": 0,
                }
            )
            continue

        reads = load_reads(region_dir / "reads" / "reads.tsv")
        methyl, cpg_columns = load_methylation(region_dir / "methylation" / "methylation.csv")
        merged = reads.merge(methyl, on="read_id", how="left", copy=False)

        for _, read_row in merged.iterrows():
            payload = {field: manifest_row.get(field, "") for field in READ_OUTPUT_FIELDS if field not in {
                "region_dir",
                "read_id",
                "read_name",
                "chr",
                "start",
                "end",
                "mapq",
                "hp",
                "alt_support",
                "is_tumor",
                "strand",
                "dataset_role",
                "harmonization_group",
                "phase1a_task",
                "phase1a_region_label",
                "phase1a_read_label",
                "phase1b_ready",
                "phase1b_blocker",
                "num_cpg_total",
                "num_cpg_observed",
                "methyl_na_fraction",
                "methyl_mean",
                "methyl_std",
                "methyl_median",
                "methyl_high_fraction",
                "methyl_low_fraction",
            }}
            payload["region_dir"] = str(region_dir)
            payload.update(provenance_payload(manifest_row))
            payload.update(
                {
                    "read_id": str(read_row.get("read_id", "")),
                    "read_name": read_row.get("read_name", ""),
                    "chr": read_row.get("chr", ""),
                    "start": read_row.get("start", ""),
                    "end": read_row.get("end", ""),
                    "mapq": read_row.get("mapq", ""),
                    "hp": read_row.get("hp", ""),
                    "alt_support": read_row.get("alt_support", ""),
                    "is_tumor": read_row.get("is_tumor", ""),
                    "strand": read_row.get("strand", ""),
                }
            )
            payload.update(
                phase1_read_context(
                    platform=str(manifest_row.get("platform", "")),
                    mode=str(manifest_row.get("mode", "")),
                    truth_status=manifest_row.get("truth_status", ""),
                    alt_support=read_row.get("alt_support", ""),
                )
            )
            payload.update(summarize_methylation_row(read_row[cpg_columns].tolist(), len(cpg_columns)))
            read_rows.append(payload)

        shard_manifest_rows.append(
            {
                "dataset_id": dataset_id,
                "dataset_label": manifest_row["dataset_label"],
                "region_key": region_key,
                "truth_status": manifest_row["truth_status"],
                "source_scope": scope,
                "Quality_Score": manifest_row["Quality_Score"],
                **provenance_payload(manifest_row),
                "region_resolve_status": "resolved",
                "region_dir": str(region_dir),
                "summary_num_reads": manifest_row["summary_num_reads"],
                "summary_num_cpgs": manifest_row["summary_num_cpgs"],
                "read_rows_exported": int(len(merged.index)),
            }
        )

    return shard_manifest_rows, read_rows


def build_summary_rows(shard_manifest_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    if not shard_manifest_rows:
        return []

    df = pd.DataFrame(shard_manifest_rows)
    summary_rows: List[Dict[str, object]] = []
    grouped = df.groupby(["dataset_id", "source_scope"], dropna=False, sort=True)
    for (dataset_id, source_scope), sub in grouped:
        summary_rows.append(
            {
                "dataset_id": dataset_id,
                "source_scope": source_scope,
                "selected_regions": int(len(sub.index)),
                "resolved_regions": int((sub["region_resolve_status"] == "resolved").sum()),
                "missing_regions": int((sub["region_resolve_status"] == "missing").sum()),
                "read_rows_exported": int(sub["read_rows_exported"].fillna(0).sum()),
            }
        )
    return summary_rows


def main() -> None:
    args = parse_args()
    output_dir = ensure_dir(Path(args.output_dir).resolve())
    manifest_df = attach_verification_provenance(
        pd.read_csv(args.manifest_tsv, sep="\t", low_memory=False)
    )
    verification_provenance = dict(manifest_df.attrs["verification_provenance"])

    dataset_ids = args.dataset or manifest_df["dataset_id"].drop_duplicates().tolist()
    scopes = args.source_scope or ["tp", "fp"]
    selected_df = select_manifest_rows(
        manifest_df=manifest_df,
        dataset_ids=dataset_ids,
        scopes=scopes,
        limit=args.max_regions_per_dataset_scope,
    )
    shard_manifest_rows, read_rows = export_shard_rows(selected_df)
    shard_summary_rows = build_summary_rows(shard_manifest_rows)

    write_tsv_rows(output_dir / "phase1_shard_manifest.tsv", SHARD_MANIFEST_FIELDS, shard_manifest_rows)
    write_tsv_rows(output_dir / "phase1_shard_summary.tsv", SHARD_SUMMARY_FIELDS, shard_summary_rows)
    write_tsv_rows(output_dir / "phase1_shard_read_training_table.tsv", READ_OUTPUT_FIELDS, read_rows)
    write_json(
        output_dir / "run_context.json",
        {
            "task": "Phase 1 manifest-driven shard export",
            "manifest_tsv": args.manifest_tsv,
            "dataset_ids": dataset_ids,
            "scopes": scopes,
            "max_regions_per_dataset_scope": args.max_regions_per_dataset_scope,
            "output_dir": str(output_dir),
            "verification_provenance": verification_provenance,
        },
    )

    notes_lines = [
        "# Phase 1 Manifest Shard Export",
        "",
        f"- manifest_tsv: `{args.manifest_tsv}`",
        f"- dataset_ids: `{', '.join(dataset_ids)}`",
        f"- scopes: `{', '.join(scopes)}`",
        f"- max_regions_per_dataset_scope: `{args.max_regions_per_dataset_scope}`",
        "",
        "## Outputs",
        "",
        "- `phase1_shard_manifest.tsv`",
        "- `phase1_shard_summary.tsv`",
        "- `phase1_shard_read_training_table.tsv`",
        "- `run_context.json`",
    ]
    (output_dir / "round_summary.md").write_text("\n".join(notes_lines) + "\n", encoding="utf-8")

    print(f"[phase1-shard] wrote {output_dir / 'phase1_shard_manifest.tsv'}")
    print(f"[phase1-shard] wrote {output_dir / 'phase1_shard_read_training_table.tsv'}")


if __name__ == "__main__":
    main()
