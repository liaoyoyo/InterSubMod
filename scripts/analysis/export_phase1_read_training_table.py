#!/usr/bin/env python3
"""Export a Phase 1 pilot read-level training table from candidate-specific InterSubMod outputs."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

from research_common import as_bool, ensure_dir, to_float, write_json, write_tsv_rows


OUTPUT_ROOT_DEFAULT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1_read_training_exporter_pilot"
)


@dataclass(frozen=True)
class DatasetConfig:
    dataset_id: str
    label: str
    sample: str
    platform: str
    mode: str
    joined_tsv: Path
    round_root: Path
    scope_note: str


DATASET_CONFIG: Dict[str, DatasetConfig] = {
    "hcc1395_5khz_paired": DatasetConfig(
        dataset_id="hcc1395_5khz_paired",
        label="HCC1395 5kHz paired",
        sample="HCC1395",
        platform="ONT_5kHz",
        mode="paired-pure",
        joined_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
            "20260308_borderline_rescue_analysis/HCC1395_paired/eval/rescue_joined_features.tsv"
        ),
        round_root=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
            "20260308_borderline_rescue_analysis/HCC1395_paired"
        ),
        scope_note="candidate-specific rescue join with lower methylation coverage ceiling",
    ),
    "hcc1395_5khz_to": DatasetConfig(
        dataset_id="hcc1395_5khz_to",
        label="HCC1395 5kHz TO",
        sample="HCC1395",
        platform="ONT_5kHz",
        mode="to-pure",
        joined_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
            "20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_joined_features.tsv"
        ),
        round_root=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
            "20260308_hcc1395_to_candidate_rescue_methylation"
        ),
        scope_note="full tagged BAM + candidate-specific rescue join",
    ),
    "hcc1395_dorado_paired": DatasetConfig(
        dataset_id="hcc1395_dorado_paired",
        label="HCC1395 DORADO paired",
        sample="HCC1395_DORADO",
        platform="ONT_Dorado",
        mode="paired-pure",
        joined_tsv=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260309_hcc1395_dorado_paired_candidate_rescue/eval/rescue_joined_features.tsv"
        ),
        round_root=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260309_hcc1395_dorado_paired_candidate_rescue"
        ),
        scope_note="candidate-specific rescue join with lower methylation coverage ceiling",
    ),
    "hcc1395_dorado_to": DatasetConfig(
        dataset_id="hcc1395_dorado_to",
        label="HCC1395 DORADO TO",
        sample="HCC1395_DORADO",
        platform="ONT_Dorado",
        mode="to",
        joined_tsv=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260311_hcc1395_dorado_to_candidate_rescue/eval/rescue_joined_features.tsv"
        ),
        round_root=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260311_hcc1395_dorado_to_candidate_rescue"
        ),
        scope_note="candidate-window subset tagged BAM + candidate-specific rescue join",
    ),
}


READ_LEVEL_FIELDS = [
    "dataset_id",
    "dataset_label",
    "sample",
    "platform",
    "mode",
    "scope_note",
    "dataset_role",
    "harmonization_group",
    "phase1a_task",
    "phase1a_region_label",
    "phase1a_read_label",
    "phase1b_ready",
    "phase1b_blocker",
    "region_key",
    "region_dir",
    "truth_status",
    "downstream_status",
    "source_scope",
    "candidate_eligible",
    "qual",
    "gq",
    "dp",
    "af",
    "ad_ref",
    "ad_alt",
    "VerificationClass",
    "DominantLabel",
    "agreement_type",
    "class_shift",
    "cluster_class",
    "label_class",
    "PairwiseMeanDist",
    "PairwiseMedianDist",
    "AlleleDelta",
    "CramersV",
    "GlobalP",
    "Quality_Score",
    "PassedGating",
    "hp_assign_rate",
    "allele_assign_rate",
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
        help="Dataset id to export. Default: all supported datasets.",
    )
    parser.add_argument(
        "--downstream-status",
        action="append",
        choices=["caller_lost_tp", "caller_removed_fp"],
        help="Restrict to specific downstream statuses. Default: caller_lost_tp + caller_removed_fp.",
    )
    parser.add_argument("--max-regions-per-dataset", type=int, default=10, help="Max exported regions per dataset.")
    parser.add_argument("--output-dir", default=str(OUTPUT_ROOT_DEFAULT), help="Output directory.")
    parser.add_argument(
        "--include-methyl-vector",
        action="store_true",
        help="Include compact methyl vector and CpG position list columns.",
    )
    return parser.parse_args()


def parse_metadata(path: Path) -> Dict[str, str]:
    info: Dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("Region ID:"):
            info["region_id"] = line.split(":", 1)[1].strip()
        elif line.startswith("Region:"):
            info["region"] = line.split(":", 1)[1].strip()
        elif line.startswith("SNV Position:"):
            info["snv_position"] = line.split(":", 1)[1].strip()
        elif line.startswith("SNV:"):
            info["snv_change"] = line.split(":", 1)[1].strip().replace(" ", "")
        elif line.startswith("Num Reads:"):
            info["num_reads"] = line.split(":", 1)[1].strip()
        elif line.startswith("Num CpG Sites:"):
            info["num_cpgs"] = line.split(":", 1)[1].strip()
    return info


def region_key_from_metadata(info: Dict[str, str]) -> Optional[str]:
    snv_position = info.get("snv_position")
    snv_change = info.get("snv_change")
    if not snv_position or not snv_change:
        return None
    chrom, pos = snv_position.split(":")
    ref, alt = [part.strip() for part in snv_change.split("->")]
    return f"{chrom}:{pos}:{ref}:{alt}"


def infer_source_scope(row: pd.Series) -> Optional[str]:
    source_scope = str(row.get("source_scope", "")).strip().lower()
    if source_scope in {"tp", "fp"}:
        return source_scope

    status = str(row.get("downstream_status", "")).strip()
    if status == "caller_lost_tp":
        return "tp"
    if status == "caller_removed_fp":
        return "fp"
    return None


def dataset_role(platform: str) -> str:
    if platform == "ONT_5kHz":
        return "discovery"
    if platform == "ONT_Dorado":
        return "validation"
    return "unknown"


def harmonization_group(platform: str, mode: str) -> str:
    return f"{platform}|{mode}"


def normalize_alt_support(value: object) -> str:
    text = str(value).strip().upper()
    if text in {"ALT", "REF"}:
        return text
    return "UNKNOWN"


def load_joined(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "candidate_eligible" in df.columns:
        df["candidate_eligible"] = df["candidate_eligible"].map(as_bool)
    else:
        df["candidate_eligible"] = False

    numeric_columns = [
        "qual",
        "gq",
        "dp",
        "af",
        "ad_ref",
        "ad_alt",
        "PairwiseMeanDist",
        "PairwiseMedianDist",
        "AlleleDelta",
        "CramersV",
        "GlobalP",
        "Quality_Score",
        "hp_assign_rate",
        "allele_assign_rate",
    ]
    for column in numeric_columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

    return df


def build_region_index(round_root: Path) -> Dict[Tuple[str, str], Path]:
    index: Dict[Tuple[str, str], Path] = {}
    scope_roots = {
        "tp": round_root / "intersubmod" / "intersubmod_tp",
        "fp": round_root / "intersubmod" / "intersubmod_fp",
    }
    for scope, root in scope_roots.items():
        if not root.exists():
            continue
        for metadata_path in root.rglob("metadata.txt"):
            key = region_key_from_metadata(parse_metadata(metadata_path))
            if key is None:
                continue
            index[(scope, key)] = metadata_path.parent
    return index


def load_reads(reads_path: Path) -> pd.DataFrame:
    reads = pd.read_csv(reads_path, sep="\t", dtype={"read_id": str})
    reads["read_id"] = reads["read_id"].astype(str)
    return reads


def load_methylation(methylation_path: Path) -> Tuple[pd.DataFrame, List[str]]:
    methyl = pd.read_csv(methylation_path, dtype={"read_id": str})
    methyl["read_id"] = methyl["read_id"].astype(str)
    cpg_columns = [column for column in methyl.columns if column != "read_id"]
    return methyl, cpg_columns


def summarize_methylation_row(
    values: Sequence[object], num_cpg_total: int, include_vector: bool
) -> Dict[str, object]:
    series = pd.Series(values, dtype="float64")
    observed = int(series.notna().sum())
    observed_values = series.dropna()

    if observed_values.empty:
        payload = {
            "num_cpg_total": num_cpg_total,
            "num_cpg_observed": 0,
            "methyl_na_fraction": 1.0 if num_cpg_total else math.nan,
            "methyl_mean": math.nan,
            "methyl_std": math.nan,
            "methyl_median": math.nan,
            "methyl_high_fraction": math.nan,
            "methyl_low_fraction": math.nan,
        }
    else:
        payload = {
            "num_cpg_total": num_cpg_total,
            "num_cpg_observed": observed,
            "methyl_na_fraction": (num_cpg_total - observed) / num_cpg_total if num_cpg_total else math.nan,
            "methyl_mean": float(observed_values.mean()),
            "methyl_std": float(observed_values.std(ddof=0)),
            "methyl_median": float(observed_values.median()),
            # 0.8/0.2 thresholds give a simple first-pass high/low methylation summary.
            "methyl_high_fraction": float((observed_values >= 0.8).mean()),
            "methyl_low_fraction": float((observed_values <= 0.2).mean()),
        }

    if include_vector:
        rounded = ["NA" if pd.isna(value) else f"{float(value):.4f}" for value in series.tolist()]
        payload["methyl_vector_compact"] = ",".join(rounded)

    return payload


def select_regions(
    df: pd.DataFrame,
    downstream_statuses: Sequence[str],
    max_regions_per_dataset: int,
) -> pd.DataFrame:
    selected = df[df["downstream_status"].isin(downstream_statuses)].copy()
    selected = selected[selected["candidate_eligible"] == True].copy()
    if "gq" in selected.columns:
        selected = selected.sort_values(["downstream_status", "gq", "qual", "af"], ascending=[True, False, False, False])
    return selected.head(max_regions_per_dataset)


def row_to_base_payload(cfg: DatasetConfig, row: pd.Series, region_dir: Path) -> Dict[str, object]:
    return {
        "dataset_id": cfg.dataset_id,
        "dataset_label": cfg.label,
        "sample": cfg.sample,
        "platform": cfg.platform,
        "mode": cfg.mode,
        "scope_note": cfg.scope_note,
        "dataset_role": dataset_role(cfg.platform),
        "harmonization_group": harmonization_group(cfg.platform, cfg.mode),
        "phase1a_task": "within_tumor_alt_support",
        "phase1a_region_label": str(row.get("truth_status", "") or ""),
        "phase1a_read_label": "",
        "phase1b_ready": False,
        "phase1b_blocker": "normal_reads_not_exported",
        "region_key": row.get("region_key", ""),
        "region_dir": str(region_dir),
        "truth_status": row.get("truth_status", ""),
        "downstream_status": row.get("downstream_status", ""),
        "source_scope": infer_source_scope(row) or "",
        "candidate_eligible": bool(row.get("candidate_eligible", False)),
        "qual": to_float(row.get("qual")),
        "gq": to_float(row.get("gq")),
        "dp": to_float(row.get("dp")),
        "af": to_float(row.get("af")),
        "ad_ref": to_float(row.get("ad_ref")),
        "ad_alt": to_float(row.get("ad_alt")),
        "VerificationClass": row.get("VerificationClass", ""),
        "DominantLabel": row.get("DominantLabel", ""),
        "agreement_type": row.get("agreement_type", ""),
        "class_shift": row.get("class_shift", ""),
        "cluster_class": row.get("cluster_class", ""),
        "label_class": row.get("label_class", ""),
        "PairwiseMeanDist": to_float(row.get("PairwiseMeanDist")),
        "PairwiseMedianDist": to_float(row.get("PairwiseMedianDist")),
        "AlleleDelta": to_float(row.get("AlleleDelta")),
        "CramersV": to_float(row.get("CramersV")),
        "GlobalP": to_float(row.get("GlobalP")),
        "Quality_Score": to_float(row.get("Quality_Score")),
        "PassedGating": row.get("PassedGating", ""),
        "hp_assign_rate": to_float(row.get("hp_assign_rate")),
        "allele_assign_rate": to_float(row.get("allele_assign_rate")),
    }


def export_dataset(
    cfg: DatasetConfig,
    downstream_statuses: Sequence[str],
    max_regions_per_dataset: int,
    include_vector: bool,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], Dict[str, object]]:
    joined = load_joined(cfg.joined_tsv)
    selected = select_regions(joined, downstream_statuses, max_regions_per_dataset)
    region_index = build_region_index(cfg.round_root)

    training_rows: List[Dict[str, object]] = []
    manifest_rows: List[Dict[str, object]] = []
    missing_regions = 0

    for _, row in selected.iterrows():
        region_key = str(row.get("region_key", ""))
        source_scope = infer_source_scope(row)
        region_dir = region_index.get((source_scope or "", region_key))
        if region_dir is None:
            missing_regions += 1
            manifest_rows.append(
                {
                    "dataset_id": cfg.dataset_id,
                    "dataset_label": cfg.label,
                    "region_key": region_key,
                    "downstream_status": row.get("downstream_status", ""),
                    "source_scope": source_scope or "",
                    "region_found": False,
                    "region_dir": "",
                    "reads_exported": 0,
                    "num_cpg_total": 0,
                }
            )
            continue

        reads = load_reads(region_dir / "reads" / "reads.tsv")
        methyl, cpg_columns = load_methylation(region_dir / "methylation" / "methylation.csv")
        num_cpg_total = len(cpg_columns)

        merged = reads.merge(methyl, on="read_id", how="left", copy=False)
        base_payload = row_to_base_payload(cfg, row, region_dir)
        cpg_position_list = ",".join(cpg_columns) if include_vector else None

        for _, read_row in merged.iterrows():
            payload = dict(base_payload)
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
                    "phase1a_read_label": normalize_alt_support(read_row.get("alt_support", "")),
                }
            )

            methyl_payload = summarize_methylation_row(read_row[cpg_columns].tolist(), num_cpg_total, include_vector)
            payload.update(methyl_payload)
            if include_vector:
                payload["cpg_position_list"] = cpg_position_list
            training_rows.append(payload)

        manifest_rows.append(
            {
                "dataset_id": cfg.dataset_id,
                "dataset_label": cfg.label,
                "region_key": region_key,
                "downstream_status": row.get("downstream_status", ""),
                "source_scope": source_scope or "",
                "region_found": True,
                "region_dir": str(region_dir),
                "reads_exported": int(len(merged.index)),
                "num_cpg_total": num_cpg_total,
            }
        )

    summary = {
        "dataset_id": cfg.dataset_id,
        "dataset_label": cfg.label,
        "selected_regions": int(len(selected.index)),
        "regions_found": int(sum(1 for row in manifest_rows if row["region_found"])),
        "missing_regions": int(missing_regions),
        "read_rows_exported": int(len(training_rows)),
    }
    return training_rows, manifest_rows, summary


def main() -> None:
    args = parse_args()
    output_dir = ensure_dir(Path(args.output_dir).resolve())
    dataset_ids = args.dataset or list(DATASET_CONFIG.keys())
    downstream_statuses = args.downstream_status or ["caller_lost_tp", "caller_removed_fp"]

    read_rows: List[Dict[str, object]] = []
    manifest_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []

    read_fields = list(READ_LEVEL_FIELDS)
    if args.include_methyl_vector:
        read_fields.extend(["methyl_vector_compact", "cpg_position_list"])

    for dataset_id in dataset_ids:
        cfg = DATASET_CONFIG[dataset_id]
        dataset_read_rows, dataset_manifest_rows, dataset_summary = export_dataset(
            cfg=cfg,
            downstream_statuses=downstream_statuses,
            max_regions_per_dataset=args.max_regions_per_dataset,
            include_vector=args.include_methyl_vector,
        )
        read_rows.extend(dataset_read_rows)
        manifest_rows.extend(dataset_manifest_rows)
        summary_rows.append(dataset_summary)

    write_tsv_rows(output_dir / "phase1_read_training_table.tsv", read_fields, read_rows)
    write_tsv_rows(
        output_dir / "phase1_region_manifest.tsv",
        [
            "dataset_id",
            "dataset_label",
            "region_key",
            "downstream_status",
            "source_scope",
            "region_found",
            "region_dir",
            "reads_exported",
            "num_cpg_total",
        ],
        manifest_rows,
    )
    write_tsv_rows(
        output_dir / "phase1_export_summary.tsv",
        ["dataset_id", "dataset_label", "selected_regions", "regions_found", "missing_regions", "read_rows_exported"],
        summary_rows,
    )

    write_json(
        output_dir / "run_context.json",
        {
            "task": "Phase 1 pilot read-level training table export",
            "datasets": dataset_ids,
            "downstream_statuses": downstream_statuses,
            "max_regions_per_dataset": args.max_regions_per_dataset,
            "include_methyl_vector": args.include_methyl_vector,
            "output_dir": str(output_dir),
        },
    )

    notes_lines = [
        "# Phase 1 Pilot Read Training Export",
        "",
        f"- datasets: `{', '.join(dataset_ids)}`",
        f"- downstream_statuses: `{', '.join(downstream_statuses)}`",
        f"- max_regions_per_dataset: `{args.max_regions_per_dataset}`",
        f"- include_methyl_vector: `{args.include_methyl_vector}`",
        "",
        "## Outputs",
        "",
        "- `phase1_read_training_table.tsv`",
        "- `phase1_region_manifest.tsv`",
        "- `phase1_export_summary.tsv`",
        "- `run_context.json`",
        "",
        "## Notes",
        "",
        "- 本腳本目前以 candidate-specific rescue rounds 為資料來源。",
        "- 第一版訓練單位為 `read × region`。",
        "- `methyl_high_fraction` 與 `methyl_low_fraction` 目前用 `>=0.8` 與 `<=0.2` 做簡單切分。",
    ]
    (output_dir / "round_summary.md").write_text("\n".join(notes_lines) + "\n", encoding="utf-8")

    print(f"[phase1-export] wrote {output_dir / 'phase1_read_training_table.tsv'}")
    print(f"[phase1-export] wrote {output_dir / 'phase1_region_manifest.tsv'}")
    print(f"[phase1-export] wrote {output_dir / 'phase1_export_summary.tsv'}")


if __name__ == "__main__":
    main()
