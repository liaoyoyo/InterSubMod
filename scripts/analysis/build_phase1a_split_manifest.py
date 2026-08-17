#!/usr/bin/env python3
"""Build a Phase 1A region-level split manifest from the baseline Phase 1 training manifest."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.analysis.build_phase1_training_manifest import (  # noqa: E402
    PROVENANCE_EXPORT_FIELDS,
    dataset_role,
    harmonization_group,
)
from scripts.analysis.research_common import (  # noqa: E402
    ensure_dir,
    write_json,
    write_tsv_rows,
)
from scripts.lib.verification_schema_contract import (  # noqa: E402
    UNKNOWN_CURRENT_CLASS,
    VERIFICATION_PROVENANCE_COLUMNS,
    SchemaContractError,
    extract_provenance_frame,
    read_evidence,
    select_current_view,
    select_legacy_view,
    select_loh_legacy,
)


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
    *PROVENANCE_EXPORT_FIELDS,
    "Quality_Score",
    "PairwiseMedianDist",
    "PassedGating",
    "region_resolve_status",
]


SUMMARY_FIELDS = [
    "dataset_id",
    "dataset_role",
    "split_role",
    "VerificationProvenanceStatus",
    "VerificationSchemaVersion",
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
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help="Explicitly authorize an H1 manifest marked UNVERSIONED_V1.",
    )
    return parser.parse_args()


def split_role_for_platform(platform: str) -> str:
    role = dataset_role(platform)
    if role == "discovery":
        return "discovery"
    if role == "validation":
        return "external_validation"
    return "unknown"


def validate_manifest_provenance(
    manifest_df: pd.DataFrame,
    allow_unversioned_v1: bool = False,
) -> str:
    """Validate the derived-manifest receipt without reconstructing missing provenance."""

    missing = [field for field in PROVENANCE_EXPORT_FIELDS if field not in manifest_df.columns]
    if missing:
        raise SchemaContractError(
            "phase1a split provenance: input manifest is missing pass-through fields: "
            + ", ".join(missing)
        )
    statuses = sorted(
        manifest_df["VerificationProvenanceStatus"].dropna().astype(str).unique().tolist()
    )
    if len(statuses) != 1:
        raise SchemaContractError(
            f"phase1a split provenance: expected one schema status, observed {statuses}"
        )
    status = statuses[0]
    if status == "V2":
        current = select_current_view(manifest_df)
        select_legacy_view(manifest_df)
        read_evidence(manifest_df)
        select_loh_legacy(manifest_df)
        known_mask = current.values != UNKNOWN_CURRENT_CLASS
        if known_mask.any():
            extract_provenance_frame(manifest_df.loc[known_mask, list(VERIFICATION_PROVENANCE_COLUMNS)])
        return status
    if status != "UNVERSIONED_V1":
        raise SchemaContractError(f"phase1a split provenance: unsupported schema status {status!r}")
    if not allow_unversioned_v1:
        raise SchemaContractError(
            "phase1a split provenance: H1 input requires --allow-unversioned-v1 authorization"
        )
    if manifest_df["VerificationSchemaVersion"].notna().any() and (
        manifest_df["VerificationSchemaVersion"].astype(str).str.strip() != ""
    ).any():
        raise SchemaContractError(
            "phase1a split provenance: UNVERSIONED_V1 rows cannot carry a schema version"
        )
    historical = manifest_df[["VerificationClass", "LOH_Subtype"]].copy()
    extract_provenance_frame(historical, allow_unversioned_v1=True)
    return status


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
                **{field: row.get(field, "") for field in PROVENANCE_EXPORT_FIELDS},
                "Quality_Score": row.get("Quality_Score", ""),
                "PairwiseMedianDist": row.get("PairwiseMedianDist", ""),
                "PassedGating": row.get("PassedGating", ""),
                "region_resolve_status": row.get("region_resolve_status", ""),
            }
        )
    return rows


def build_summary_rows(split_df: pd.DataFrame) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    grouped = split_df.groupby(
        [
            "dataset_id",
            "dataset_role",
            "split_role",
            "VerificationProvenanceStatus",
            "VerificationSchemaVersion",
        ],
        dropna=False,
        sort=True,
    )
    for (dataset_id, dataset_role_value, split_role, provenance_status, schema_version), sub in grouped:
        passed_gating = sub["PassedGating"].astype(str).str.lower().isin({"true", "1"})
        rows.append(
            {
                "dataset_id": dataset_id,
                "dataset_role": dataset_role_value,
                "split_role": split_role,
                "VerificationProvenanceStatus": provenance_status,
                "VerificationSchemaVersion": schema_version,
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
    manifest_df = pd.read_csv(
        args.manifest_tsv,
        sep="\t",
        keep_default_na=False,
        low_memory=False,
    )
    provenance_status = validate_manifest_provenance(
        manifest_df,
        allow_unversioned_v1=args.allow_unversioned_v1,
    )

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
            "verification_provenance_status": provenance_status,
            "allow_unversioned_v1": args.allow_unversioned_v1,
            "verification_provenance_fields": PROVENANCE_EXPORT_FIELDS,
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
