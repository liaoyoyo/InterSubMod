#!/usr/bin/env python3
"""Build a Phase 1 baseline training manifest from four baseline InterSubMod outputs."""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import pandas as pd

from research_common import ensure_dir, to_float, write_json, write_tsv_rows

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

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


OUTPUT_ROOT_DEFAULT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1_training_manifest_v1"
)

PROVENANCE_METADATA_FIELDS = [
    "VerificationProvenanceStatus",
    "VerificationProvenanceSourceField",
    "VerificationProvenanceWarnings",
]
PROVENANCE_EXPORT_FIELDS = list(VERIFICATION_PROVENANCE_COLUMNS) + PROVENANCE_METADATA_FIELDS


@dataclass(frozen=True)
class BaselineDatasetConfig:
    dataset_id: str
    label: str
    sample: str
    platform: str
    mode: str
    scope_note: str
    tp_summary_tsv: Path
    fp_summary_tsv: Path
    tp_region_root: Path
    fp_region_root: Path


DATASET_CONFIG: Dict[str, BaselineDatasetConfig] = {
    "hcc1395_5khz_paired": BaselineDatasetConfig(
        dataset_id="hcc1395_5khz_paired",
        label="HCC1395 5kHz paired",
        sample="HCC1395",
        platform="ONT_5kHz",
        mode="paired-pure",
        scope_note="canonical paired full complete matrix baseline",
        tp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
            "20260314_HCC1395_paired_full_full_complete_matrix/intersubmod_tp/significance_summary.csv"
        ),
        fp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
            "20260314_HCC1395_paired_full_full_complete_matrix/intersubmod_fp/significance_summary.csv"
        ),
        tp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
            "20260314_HCC1395_paired_full_full_complete_matrix/intersubmod_tp/filtered_snv_tp"
        ),
        fp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
            "20260314_HCC1395_paired_full_full_complete_matrix/intersubmod_fp/filtered_snv_fp"
        ),
    ),
    "hcc1395_5khz_to": BaselineDatasetConfig(
        dataset_id="hcc1395_5khz_to",
        label="HCC1395 5kHz TO",
        sample="HCC1395",
        platform="ONT_5kHz",
        mode="to-pure",
        scope_note="synthesis TO pilot baseline",
        tp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
            "20260315_hcc1395_to_pilot/step05_intersubmod/intersubmod_tp/significance_summary.csv"
        ),
        fp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
            "20260315_hcc1395_to_pilot/step05_intersubmod/intersubmod_fp/significance_summary.csv"
        ),
        tp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
            "20260315_hcc1395_to_pilot/step05_intersubmod/intersubmod_tp/filtered_snv_tp"
        ),
        fp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
            "20260315_hcc1395_to_pilot/step05_intersubmod/intersubmod_fp/filtered_snv_fp"
        ),
    ),
    "hcc1395_dorado_paired": BaselineDatasetConfig(
        dataset_id="hcc1395_dorado_paired",
        label="HCC1395 DORADO paired",
        sample="HCC1395_DORADO",
        platform="ONT_Dorado",
        mode="paired-pure",
        scope_note="canonical paired full complete matrix baseline",
        tp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395_DORADO/paired_full/"
            "20260315_HCC1395_DORADO_paired_full_full_complete_matrix/intersubmod_tp/significance_summary.csv"
        ),
        fp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395_DORADO/paired_full/"
            "20260315_HCC1395_DORADO_paired_full_full_complete_matrix/intersubmod_fp/significance_summary.csv"
        ),
        tp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395_DORADO/paired_full/"
            "20260315_HCC1395_DORADO_paired_full_full_complete_matrix/intersubmod_tp/filtered_snv_tp"
        ),
        fp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395_DORADO/paired_full/"
            "20260315_HCC1395_DORADO_paired_full_full_complete_matrix/intersubmod_fp/filtered_snv_fp"
        ),
    ),
    "hcc1395_dorado_to": BaselineDatasetConfig(
        dataset_id="hcc1395_dorado_to",
        label="HCC1395 DORADO TO",
        sample="HCC1395_DORADO",
        platform="ONT_Dorado",
        mode="to-pure",
        scope_note="synthesis DORADO TO pilot baseline",
        tp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
            "20260315_hcc1395_dorado_to_pilot/step05_intersubmod/intersubmod_tp/significance_summary.csv"
        ),
        fp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
            "20260315_hcc1395_dorado_to_pilot/step05_intersubmod/intersubmod_fp/significance_summary.csv"
        ),
        tp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
            "20260315_hcc1395_dorado_to_pilot/step05_intersubmod/intersubmod_tp/filtered_snv_tp"
        ),
        fp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
            "20260315_hcc1395_dorado_to_pilot/step05_intersubmod/intersubmod_fp/filtered_snv_fp"
        ),
    ),
    "colo829_pao_paired": BaselineDatasetConfig(
        dataset_id="colo829_pao_paired",
        label="COLO829 PAO paired",
        sample="COLO829",
        platform="ONT_PAO",
        mode="paired-pure",
        scope_note="canonical paired full complete matrix baseline",
        tp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/COLO829/paired_full/"
            "20260315_COLO829_paired_full_full_complete_matrix/intersubmod_tp/significance_summary.csv"
        ),
        fp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/COLO829/paired_full/"
            "20260315_COLO829_paired_full_full_complete_matrix/intersubmod_fp/significance_summary.csv"
        ),
        tp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/COLO829/paired_full/"
            "20260315_COLO829_paired_full_full_complete_matrix/intersubmod_tp/filtered_snv_tp"
        ),
        fp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/COLO829/paired_full/"
            "20260315_COLO829_paired_full_full_complete_matrix/intersubmod_fp/filtered_snv_fp"
        ),
    ),
    "h1437_google_paired": BaselineDatasetConfig(
        dataset_id="h1437_google_paired",
        label="H1437 Google paired",
        sample="H1437",
        platform="ONT_Google",
        mode="paired-pure",
        scope_note="canonical paired full complete matrix baseline",
        tp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/H1437/paired_full/"
            "20260315_H1437_paired_full_full_complete_matrix/intersubmod_tp/significance_summary.csv"
        ),
        fp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/H1437/paired_full/"
            "20260315_H1437_paired_full_full_complete_matrix/intersubmod_fp/significance_summary.csv"
        ),
        tp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/H1437/paired_full/"
            "20260315_H1437_paired_full_full_complete_matrix/intersubmod_tp/filtered_snv_tp"
        ),
        fp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/H1437/paired_full/"
            "20260315_H1437_paired_full_full_complete_matrix/intersubmod_fp/filtered_snv_fp"
        ),
    ),
    "h2009_google_paired": BaselineDatasetConfig(
        dataset_id="h2009_google_paired",
        label="H2009 Google paired",
        sample="H2009",
        platform="ONT_Google",
        mode="paired-pure",
        scope_note="canonical paired full complete matrix baseline",
        tp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/H2009/paired_full/"
            "20260315_H2009_paired_full_full_complete_matrix/intersubmod_tp/significance_summary.csv"
        ),
        fp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/H2009/paired_full/"
            "20260315_H2009_paired_full_full_complete_matrix/intersubmod_fp/significance_summary.csv"
        ),
        tp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/H2009/paired_full/"
            "20260315_H2009_paired_full_full_complete_matrix/intersubmod_tp/filtered_snv_tp"
        ),
        fp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/H2009/paired_full/"
            "20260315_H2009_paired_full_full_complete_matrix/intersubmod_fp/filtered_snv_fp"
        ),
    ),
    "hcc1937_google_paired": BaselineDatasetConfig(
        dataset_id="hcc1937_google_paired",
        label="HCC1937 Google paired",
        sample="HCC1937",
        platform="ONT_Google",
        mode="paired-pure",
        scope_note="canonical paired full complete matrix baseline",
        tp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1937/paired_full/"
            "20260315_HCC1937_paired_full_full_complete_matrix/intersubmod_tp/significance_summary.csv"
        ),
        fp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1937/paired_full/"
            "20260315_HCC1937_paired_full_full_complete_matrix/intersubmod_fp/significance_summary.csv"
        ),
        tp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1937/paired_full/"
            "20260315_HCC1937_paired_full_full_complete_matrix/intersubmod_tp/filtered_snv_tp"
        ),
        fp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1937/paired_full/"
            "20260315_HCC1937_paired_full_full_complete_matrix/intersubmod_fp/filtered_snv_fp"
        ),
    ),
    "hcc1954_google_paired": BaselineDatasetConfig(
        dataset_id="hcc1954_google_paired",
        label="HCC1954 Google paired",
        sample="HCC1954",
        platform="ONT_Google",
        mode="paired-pure",
        scope_note="canonical paired full complete matrix baseline",
        tp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1954/paired_full/"
            "20260315_HCC1954_paired_full_full_complete_matrix/intersubmod_tp/significance_summary.csv"
        ),
        fp_summary_tsv=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1954/paired_full/"
            "20260315_HCC1954_paired_full_full_complete_matrix/intersubmod_fp/significance_summary.csv"
        ),
        tp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1954/paired_full/"
            "20260315_HCC1954_paired_full_full_complete_matrix/intersubmod_tp/filtered_snv_tp"
        ),
        fp_region_root=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1954/paired_full/"
            "20260315_HCC1954_paired_full_full_complete_matrix/intersubmod_fp/filtered_snv_fp"
        ),
    ),
}


MANIFEST_FIELDS = [
    "dataset_id",
    "dataset_label",
    "sample",
    "platform",
    "mode",
    "source_profile",
    "scope_note",
    "region_key",
    "region_id",
    "truth_status",
    "downstream_status",
    "source_scope",
    "region_resolve_status",
    "region_dir",
    "summary_num_reads",
    "summary_num_cpgs",
    "GlobalP",
    "CramersV",
    "HeuristicScore",
    "PassedGating",
    "PairwiseMeanDist",
    "PairwiseMedianDist",
    "AlleleDelta",
    "DominantLabel",
    "Potential_LOH",
    "Coverage_Category",
    *PROVENANCE_EXPORT_FIELDS,
    "Quality_Score",
    "Quality_Tier",
    "Significant",
    "SuggestFilter",
    "read_export_selected",
    "read_rows_exported",
]


READ_LEVEL_FIELDS = [
    "dataset_id",
    "dataset_label",
    "sample",
    "platform",
    "mode",
    "source_profile",
    "scope_note",
    "dataset_role",
    "harmonization_group",
    "phase1a_task",
    "phase1a_region_label",
    "phase1a_read_label",
    "phase1b_ready",
    "phase1b_blocker",
    "region_key",
    "region_id",
    "truth_status",
    "downstream_status",
    "source_scope",
    "region_dir",
    "summary_num_reads",
    "summary_num_cpgs",
    "GlobalP",
    "CramersV",
    "HeuristicScore",
    "PassedGating",
    "PairwiseMeanDist",
    "PairwiseMedianDist",
    "AlleleDelta",
    "DominantLabel",
    "Potential_LOH",
    "Coverage_Category",
    *PROVENANCE_EXPORT_FIELDS,
    "Quality_Score",
    "Quality_Tier",
    "Significant",
    "SuggestFilter",
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
    "num_cpg_total",
    "num_cpg_observed",
    "methyl_na_fraction",
    "methyl_mean",
    "methyl_std",
    "methyl_median",
    "methyl_high_fraction",
    "methyl_low_fraction",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        action="append",
        choices=sorted(DATASET_CONFIG.keys()),
        help="Dataset id to include. Default: all four baseline datasets.",
    )
    parser.add_argument(
        "--read-export-limit-per-scope",
        type=int,
        default=0,
        help="Number of regions to expand into read-level rows per dataset per scope. Default: 0.",
    )
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help="Explicitly authorize legacy H1 input without VerificationSchemaVersion.",
    )
    parser.add_argument("--output-dir", default=str(OUTPUT_ROOT_DEFAULT), help="Output directory.")
    return parser.parse_args()


def parse_metadata(path: Path) -> Dict[str, str]:
    info: Dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("Region ID:"):
            info["region_id"] = line.split(":", 1)[1].strip()
        elif line.startswith("SNV Position:"):
            info["snv_position"] = line.split(":", 1)[1].strip()
        elif line.startswith("SNV:"):
            info["snv_change"] = line.split(":", 1)[1].strip().replace(" ", "")
    return info


def region_key_from_metadata(info: Dict[str, str]) -> Optional[str]:
    snv_position = info.get("snv_position")
    snv_change = info.get("snv_change")
    if not snv_position or not snv_change:
        return None
    chrom, pos = snv_position.split(":")
    ref, alt = [part.strip() for part in snv_change.split("->")]
    return f"{chrom}:{pos}:{ref}:{alt}"


def resolve_region_dir(root: Path, region_key: str) -> Optional[Path]:
    chrom, pos, _ref, _alt = region_key.split(":")
    pos_dir = root / chrom / f"{chrom}_{pos}"
    if pos_dir.exists():
        candidate_dirs = [path for path in pos_dir.iterdir() if path.is_dir()]
        if len(candidate_dirs) == 1:
            return candidate_dirs[0]
        for candidate_dir in candidate_dirs:
            metadata_path = candidate_dir / "metadata.txt"
            if metadata_path.exists() and region_key_from_metadata(parse_metadata(metadata_path)) == region_key:
                return candidate_dir

    for metadata_path in root.rglob("metadata.txt"):
        if region_key_from_metadata(parse_metadata(metadata_path)) == region_key:
            return metadata_path.parent
    return None


def attach_input_provenance(df: pd.DataFrame, allow_unversioned_v1: bool = False) -> pd.DataFrame:
    """Validate one source table and attach explicit pass-through receipt columns."""

    selected = df.copy()
    if "VerificationSchemaVersion" in selected.columns:
        current = select_current_view(selected)
        select_legacy_view(selected)
        read_evidence(selected)
        select_loh_legacy(selected)
        missing = [field for field in VERIFICATION_PROVENANCE_COLUMNS if field not in selected.columns]
        if missing:
            raise SchemaContractError(
                "phase1 manifest provenance: missing required columns: " + ", ".join(missing)
            )
        known_mask = current.values != UNKNOWN_CURRENT_CLASS
        if known_mask.any():
            extract_provenance_frame(selected.loc[known_mask])
        schema_status = "V2"
        source_field = current.field
        warning_payload = {
            "unknown_current_counts": current.unknown_counts,
            "messages": list(current.warning_messages),
        }
    else:
        if not allow_unversioned_v1:
            raise SchemaContractError(
                "phase1 manifest provenance: VerificationSchemaVersion missing; "
                "rerun with --allow-unversioned-v1 only for authorized historical H1 input"
            )
        select_loh_legacy(selected, allow_unversioned_v1=True)
        provenance = extract_provenance_frame(selected, allow_unversioned_v1=True)
        schema_status = "UNVERSIONED_V1"
        source_field = str(provenance.attrs.get("selection_field", "VerificationClass"))
        warning_payload = {
            "authorization": "--allow-unversioned-v1",
            "messages": list(provenance.attrs.get("warnings", [])),
        }
        for field in VERIFICATION_PROVENANCE_COLUMNS:
            if field not in selected.columns:
                selected[field] = ""

    selected["VerificationProvenanceStatus"] = schema_status
    selected["VerificationProvenanceSourceField"] = source_field
    selected["VerificationProvenanceWarnings"] = json.dumps(
        warning_payload, ensure_ascii=False, sort_keys=True
    )
    return selected


def provenance_export_value(row: pd.Series, field: str) -> object:
    """Serialize validated missing evidence without emitting the invalid literal ``nan``."""
    value = row.get(field, "")
    if not pd.isna(value):
        return value
    if (
        row.get("VerificationProvenanceStatus") == "V2"
        and field in {"WithinHPSupport", "DispersionWarning"}
    ):
        return "NA"
    return ""


def load_summary(path: Path, scope: str, allow_unversioned_v1: bool = False) -> pd.DataFrame:
    df = pd.read_csv(path, keep_default_na=False, low_memory=False)
    df = attach_input_provenance(df, allow_unversioned_v1=allow_unversioned_v1)
    df["region_key"] = df.apply(lambda row: f"{row['Chr']}:{row['Pos']}:{row['Ref']}:{row['Alt']}", axis=1)
    df["source_scope"] = scope
    df["truth_status"] = "TP" if scope == "tp" else "FP"
    df["downstream_status"] = "baseline_kept_tp" if scope == "tp" else "baseline_kept_fp"
    return df


def build_manifest_for_dataset(
    cfg: BaselineDatasetConfig,
    allow_unversioned_v1: bool = False,
) -> List[Dict[str, object]]:
    summary = pd.concat(
        [
            load_summary(cfg.tp_summary_tsv, "tp", allow_unversioned_v1),
            load_summary(cfg.fp_summary_tsv, "fp", allow_unversioned_v1),
        ],
        ignore_index=True,
    )
    provenance_statuses = sorted(summary["VerificationProvenanceStatus"].dropna().unique().tolist())
    if len(provenance_statuses) != 1:
        raise SchemaContractError(
            f"phase1 manifest provenance: mixed schema statuses within dataset: {provenance_statuses}"
        )
    rows: List[Dict[str, object]] = []
    for _, row in summary.iterrows():
        rows.append(
            {
                "dataset_id": cfg.dataset_id,
                "dataset_label": cfg.label,
                "sample": cfg.sample,
                "platform": cfg.platform,
                "mode": cfg.mode,
                "source_profile": "baseline_full",
                "scope_note": cfg.scope_note,
                "region_key": row.get("region_key", ""),
                "region_id": row.get("RegionID", ""),
                "truth_status": row.get("truth_status", ""),
                "downstream_status": row.get("downstream_status", ""),
                "source_scope": row.get("source_scope", ""),
                "region_resolve_status": "not_attempted",
                "region_dir": "",
                "summary_num_reads": row.get("NumReads", ""),
                "summary_num_cpgs": row.get("NumCpGs", ""),
                "GlobalP": to_float(row.get("GlobalP")),
                "CramersV": to_float(row.get("CramersV")),
                "HeuristicScore": to_float(row.get("HeuristicScore")),
                "PassedGating": row.get("PassedGating", ""),
                "PairwiseMeanDist": to_float(row.get("PairwiseMeanDist")),
                "PairwiseMedianDist": to_float(row.get("PairwiseMedianDist")),
                "AlleleDelta": to_float(row.get("AlleleDelta")),
                "DominantLabel": row.get("DominantLabel", ""),
                "Potential_LOH": row.get("Potential_LOH", ""),
                "Coverage_Category": row.get("Coverage_Category", ""),
                **{field: provenance_export_value(row, field) for field in PROVENANCE_EXPORT_FIELDS},
                "Quality_Score": to_float(row.get("Quality_Score")),
                "Quality_Tier": row.get("Quality_Tier", ""),
                "Significant": row.get("Significant", ""),
                "SuggestFilter": row.get("SuggestFilter", ""),
                "read_export_selected": False,
                "read_rows_exported": 0,
            }
        )
    return rows


def choose_regions_for_read_export(manifest_df: pd.DataFrame, limit_per_scope: int) -> pd.DataFrame:
    if limit_per_scope <= 0:
        return manifest_df.iloc[0:0].copy()
    chosen_frames: List[pd.DataFrame] = []
    for dataset_id in manifest_df["dataset_id"].drop_duplicates().tolist():
        for scope in ["tp", "fp"]:
            sub = manifest_df[
                (manifest_df["dataset_id"] == dataset_id) & (manifest_df["source_scope"] == scope)
            ].copy()
            if sub.empty:
                continue
            sub = sub.sort_values(
                ["Quality_Score", "summary_num_reads", "PairwiseMedianDist"],
                ascending=[False, False, False],
            )
            chosen_frames.append(sub.head(limit_per_scope))
    if not chosen_frames:
        return manifest_df.iloc[0:0].copy()
    return pd.concat(chosen_frames, ignore_index=True)


def load_reads(reads_path: Path) -> pd.DataFrame:
    reads = pd.read_csv(reads_path, sep="\t", dtype={"read_id": str})
    reads["read_id"] = reads["read_id"].astype(str)
    return reads


def load_methylation(methylation_path: Path) -> tuple[pd.DataFrame, List[str]]:
    methyl = pd.read_csv(methylation_path, dtype={"read_id": str})
    methyl["read_id"] = methyl["read_id"].astype(str)
    cpg_columns = [column for column in methyl.columns if column != "read_id"]
    return methyl, cpg_columns


def summarize_methylation_row(values: Sequence[object], num_cpg_total: int) -> Dict[str, object]:
    series = pd.Series(values, dtype="float64")
    observed = int(series.notna().sum())
    observed_values = series.dropna()
    if observed_values.empty:
        return {
            "num_cpg_total": num_cpg_total,
            "num_cpg_observed": 0,
            "methyl_na_fraction": 1.0 if num_cpg_total else math.nan,
            "methyl_mean": math.nan,
            "methyl_std": math.nan,
            "methyl_median": math.nan,
            "methyl_high_fraction": math.nan,
            "methyl_low_fraction": math.nan,
        }
    return {
        "num_cpg_total": num_cpg_total,
        "num_cpg_observed": observed,
        "methyl_na_fraction": (num_cpg_total - observed) / num_cpg_total if num_cpg_total else math.nan,
        "methyl_mean": float(observed_values.mean()),
        "methyl_std": float(observed_values.std(ddof=0)),
        "methyl_median": float(observed_values.median()),
        "methyl_high_fraction": float((observed_values >= 0.8).mean()),
        "methyl_low_fraction": float((observed_values <= 0.2).mean()),
    }


def dataset_role(platform: str) -> str:
    if platform == "ONT_5kHz":
        return "discovery"
    if platform in {"ONT_Dorado", "ONT_PAO", "ONT_Google"}:
        return "validation"
    return "unknown"


def harmonization_group(platform: str, mode: str) -> str:
    return f"{platform}|{mode}"


def normalize_alt_support(value: object) -> str:
    text = str(value).strip().upper()
    if text in {"ALT", "REF"}:
        return text
    return "UNKNOWN"


def phase1_read_context(platform: str, mode: str, truth_status: object, alt_support: object) -> Dict[str, object]:
    return {
        "dataset_role": dataset_role(platform),
        "harmonization_group": harmonization_group(platform, mode),
        "phase1a_task": "within_tumor_alt_support",
        "phase1a_region_label": str(truth_status or ""),
        "phase1a_read_label": normalize_alt_support(alt_support),
        "phase1b_ready": False,
        "phase1b_blocker": "normal_reads_not_exported",
    }


def resolve_and_expand_reads(
    cfg: BaselineDatasetConfig,
    selected_manifest_df: pd.DataFrame,
) -> tuple[List[Dict[str, object]], Dict[str, Dict[str, object]]]:
    read_rows: List[Dict[str, object]] = []
    region_updates: Dict[str, Dict[str, object]] = {}

    for _, manifest_row in selected_manifest_df.iterrows():
        region_key = str(manifest_row["region_key"])
        scope = str(manifest_row["source_scope"])
        search_root = cfg.tp_region_root if scope == "tp" else cfg.fp_region_root
        region_dir = resolve_region_dir(search_root, region_key)
        if region_dir is None:
            region_updates[region_key] = {
                "region_resolve_status": "missing",
                "region_dir": "",
                "read_rows_exported": 0,
            }
            continue

        reads = load_reads(region_dir / "reads" / "reads.tsv")
        methyl, cpg_columns = load_methylation(region_dir / "methylation" / "methylation.csv")
        merged = reads.merge(methyl, on="read_id", how="left", copy=False)

        row_count = 0
        for _, read_row in merged.iterrows():
            payload = {field: manifest_row.get(field, "") for field in READ_LEVEL_FIELDS if field not in {
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
            row_count += 1

        region_updates[region_key] = {
            "region_resolve_status": "resolved",
            "region_dir": str(region_dir),
            "read_rows_exported": row_count,
        }

    return read_rows, region_updates


def build_summary_rows(manifest_df: pd.DataFrame) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for dataset_id in manifest_df["dataset_id"].drop_duplicates().tolist():
        sub = manifest_df[manifest_df["dataset_id"] == dataset_id].copy()
        rows.append(
            {
                "dataset_id": dataset_id,
                "dataset_label": sub["dataset_label"].iloc[0],
                "regions_total": int(len(sub.index)),
                "tp_regions": int((sub["source_scope"] == "tp").sum()),
                "fp_regions": int((sub["source_scope"] == "fp").sum()),
                "resolved_regions": int((sub["region_resolve_status"] == "resolved").sum()),
                "missing_regions": int((sub["region_resolve_status"] == "missing").sum()),
                "not_attempted_regions": int((sub["region_resolve_status"] == "not_attempted").sum()),
                "read_export_regions": int((sub["read_export_selected"] == True).sum()),
                "read_rows_exported": int(pd.to_numeric(sub["read_rows_exported"], errors="coerce").fillna(0).sum()),
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    output_dir = ensure_dir(Path(args.output_dir).resolve())
    dataset_ids = args.dataset or list(DATASET_CONFIG.keys())

    manifest_rows: List[Dict[str, object]] = []
    for dataset_id in dataset_ids:
        manifest_rows.extend(
            build_manifest_for_dataset(
                DATASET_CONFIG[dataset_id],
                allow_unversioned_v1=args.allow_unversioned_v1,
            )
        )

    manifest_df = pd.DataFrame(manifest_rows)
    provenance_statuses = sorted(
        manifest_df["VerificationProvenanceStatus"].dropna().unique().tolist()
    )
    if len(provenance_statuses) != 1:
        raise SchemaContractError(
            f"phase1 manifest provenance: mixed schema statuses across datasets: {provenance_statuses}"
        )
    selected_for_read_export = choose_regions_for_read_export(manifest_df, args.read_export_limit_per_scope)

    read_rows: List[Dict[str, object]] = []
    for dataset_id in dataset_ids:
        cfg = DATASET_CONFIG[dataset_id]
        selected_sub = selected_for_read_export[selected_for_read_export["dataset_id"] == dataset_id].copy()
        if selected_sub.empty:
            continue
        dataset_read_rows, region_updates = resolve_and_expand_reads(cfg, selected_sub)
        read_rows.extend(dataset_read_rows)
        for region_key, update in region_updates.items():
            mask = (manifest_df["dataset_id"] == dataset_id) & (manifest_df["region_key"] == region_key)
            manifest_df.loc[mask, "read_export_selected"] = True
            manifest_df.loc[mask, "region_resolve_status"] = update["region_resolve_status"]
            manifest_df.loc[mask, "region_dir"] = update["region_dir"]
            manifest_df.loc[mask, "read_rows_exported"] = update["read_rows_exported"]

    summary_rows = build_summary_rows(manifest_df)

    write_tsv_rows(
        output_dir / "phase1_training_manifest.tsv",
        MANIFEST_FIELDS,
        manifest_df[MANIFEST_FIELDS].to_dict(orient="records"),
    )
    write_tsv_rows(
        output_dir / "phase1_manifest_summary.tsv",
        [
            "dataset_id",
            "dataset_label",
            "regions_total",
            "tp_regions",
            "fp_regions",
            "resolved_regions",
            "missing_regions",
            "not_attempted_regions",
            "read_export_regions",
            "read_rows_exported",
        ],
        summary_rows,
    )
    if read_rows:
        write_tsv_rows(output_dir / "phase1_read_training_table.tsv", READ_LEVEL_FIELDS, read_rows)

    write_json(
        output_dir / "run_context.json",
        {
            "task": "Phase 1 baseline training manifest build",
            "datasets": dataset_ids,
            "read_export_limit_per_scope": args.read_export_limit_per_scope,
            "output_dir": str(output_dir),
            "strategy": "summary-first manifest + selected-region resolve",
            "verification_provenance_status": provenance_statuses[0],
            "allow_unversioned_v1": args.allow_unversioned_v1,
            "verification_provenance_fields": PROVENANCE_EXPORT_FIELDS,
        },
    )

    notes_lines = [
        "# Phase 1 Baseline Training Manifest",
        "",
        f"- datasets: `{', '.join(dataset_ids)}`",
        f"- read_export_limit_per_scope: `{args.read_export_limit_per_scope}`",
        "",
        "## Outputs",
        "",
        "- `phase1_training_manifest.tsv`",
        "- `phase1_manifest_summary.tsv`",
        "- `run_context.json`",
    ]
    if read_rows:
        notes_lines.append("- `phase1_read_training_table.tsv`")
    notes_lines.extend(
        [
            "",
            "## Notes",
            "",
            "- 這份 manifest 以四個 baseline dataset 的 `significance_summary.csv` 為主軸，先完成全量 region 清單。",
            "- `region_dir` 不做全量預先掃描，只對 read export 選中的少量 region 進行 resolve。",
            "- `truth_status` 由 `tp/fp` scope 直接繼承。",
            "- `downstream_status` 在 baseline manifest 中固定為 `baseline_kept_tp` 或 `baseline_kept_fp`。",
        ]
    )
    (output_dir / "round_summary.md").write_text("\n".join(notes_lines) + "\n", encoding="utf-8")

    print(f"[phase1-manifest] wrote {output_dir / 'phase1_training_manifest.tsv'}")
    print(f"[phase1-manifest] wrote {output_dir / 'phase1_manifest_summary.tsv'}")
    if read_rows:
        print(f"[phase1-manifest] wrote {output_dir / 'phase1_read_training_table.tsv'}")


if __name__ == "__main__":
    main()
