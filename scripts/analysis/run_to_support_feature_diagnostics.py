#!/usr/bin/env python3
"""Select representative TO support-feature candidates and run diagnostics."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

LOCAL_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(LOCAL_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(LOCAL_REPO_ROOT))

from scripts.lib.verification_schema_contract import (
    VERIFICATION_PROVENANCE_COLUMNS,
    extract_provenance_frame,
)


REPO_ROOT = Path("/big8_disk/liaoyoyo2001/InterSubMod")
EXPORT_DIAG = REPO_ROOT / "scripts" / "analysis" / "export_region_diagnostics.py"
SAMTOOLS_SNAPSHOT = REPO_ROOT / "scripts" / "analysis" / "region_samtools_snapshot.sh"
SUMMARIZE_SNAPSHOTS = REPO_ROOT / "scripts" / "analysis" / "summarize_samtools_snapshots.py"


DATASET_CONFIG = {
    "hcc1395_5khz_to": {
        "label": "HCC1395 5kHz TO",
        "joined_tsv": Path(
            "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
            "20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_joined_features.tsv"
        ),
        "tp_root": Path(
            "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
            "20260308_hcc1395_to_candidate_rescue_methylation/intersubmod/intersubmod_tp"
        ),
        "fp_root": Path(
            "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
            "20260308_hcc1395_to_candidate_rescue_methylation/intersubmod/intersubmod_fp"
        ),
        "bam": Path(
            "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
            "20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam"
        ),
        "pairwise_mode": "high",
    },
    "hcc1395_dorado_to": {
        "label": "HCC1395 DORADO TO",
        "joined_tsv": Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260311_hcc1395_dorado_to_candidate_rescue/eval/rescue_joined_features.tsv"
        ),
        "tp_root": Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260311_hcc1395_dorado_to_candidate_rescue/intersubmod/intersubmod_tp"
        ),
        "fp_root": Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260311_hcc1395_dorado_to_candidate_rescue/intersubmod/intersubmod_fp"
        ),
        "bam": Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260311_hcc1395_dorado_to_candidate_rescue/step03_longphase_to/tumor_candidate_windows_tagged.bam"
        ),
        "pairwise_mode": "low",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        action="append",
        choices=sorted(DATASET_CONFIG.keys()),
        help="Dataset id to run. Default: all supported TO datasets.",
    )
    parser.add_argument("--top-k", type=int, default=3, help="Top candidates per feature/status.")
    parser.add_argument("--max-reads", type=int, default=120, help="Max reads for samtools snapshot.")
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    return parser.parse_args()


def attach_verification_provenance(frame: pd.DataFrame) -> pd.DataFrame:
    """Require complete, internally consistent v2 provenance before selection."""
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


def load_joined(path: Path) -> pd.DataFrame:
    df = attach_verification_provenance(pd.read_csv(path, sep="\t"))
    verification_provenance = dict(df.attrs["verification_provenance"])
    df = df[df["candidate_eligible"] == True].copy()
    df.attrs["verification_provenance"] = verification_provenance
    for col in ["Quality_Score", "PairwiseMedianDist", "hp_assign_rate", "af", "gq", "qual"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def normalize_region_key(key: str) -> str:
    return key.replace(":", "_")


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
    return info


def region_key_from_info(info: Dict[str, str]) -> str:
    chrom, pos = info["snv_position"].split(":")
    ref, alt = [part.strip() for part in info["snv_change"].split("->")]
    return f"{chrom}:{pos}:{ref}:{alt}"


def build_region_index(search_root: Path) -> Dict[str, Dict[str, str]]:
    index: Dict[str, Dict[str, str]] = {}
    for metadata in search_root.rglob("metadata.txt"):
        info = parse_metadata(metadata)
        if "snv_position" not in info or "snv_change" not in info:
            continue
        key = region_key_from_info(info)
        index[key] = {
            "region_dir": str(metadata.parent),
            "region": info.get("region", ""),
            "region_id": info.get("region_id", ""),
        }
    return index


def choose_rows(df: pd.DataFrame, feature: str, status: str, pairwise_mode: str, top_k: int) -> pd.DataFrame:
    sub = df[df["downstream_status"] == status].copy()
    sub = sub.dropna(subset=[feature])
    if feature == "Quality_Score":
        sub = sub[sub[feature] > 0]
        ascending = False
    elif feature == "hp_assign_rate":
        sub = sub[sub[feature] > 0]
        ascending = False
    elif feature == "PairwiseMedianDist":
        sub = sub[sub[feature] > 0]
        ascending = pairwise_mode == "low"
    else:
        ascending = False
    sub = sub.sort_values([feature, "gq", "qual", "af"], ascending=[ascending, False, False, False])
    return sub.head(top_k)


def feature_note(dataset_id: str, feature: str, row: pd.Series) -> str:
    if feature == "PairwiseMedianDist":
        if DATASET_CONFIG[dataset_id]["pairwise_mode"] == "high":
            return "high_pairwise_support"
        return "low_pairwise_support"
    if feature == "Quality_Score":
        return "high_quality_support"
    if feature == "hp_assign_rate":
        return "high_hp_assign_support"
    return feature


def run_command(cmd: List[str]) -> None:
    subprocess.run(cmd, check=True)


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: List[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_dir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    dataset_ids = args.dataset or sorted(DATASET_CONFIG.keys())

    manifest_rows: List[Dict[str, object]] = []
    selected_rows: Dict[Tuple[str, str], Dict[str, object]] = {}

    for dataset_id in dataset_ids:
        cfg = DATASET_CONFIG[dataset_id]
        df = load_joined(cfg["joined_tsv"])
        for status in ["caller_lost_tp", "caller_removed_fp"]:
            for feature in ["Quality_Score", "PairwiseMedianDist", "hp_assign_rate"]:
                chosen = choose_rows(df, feature, status, cfg["pairwise_mode"], args.top_k)
                for rank, (_, row) in enumerate(chosen.iterrows(), start=1):
                    selected_key = (dataset_id, row["region_key"])
                    note = feature_note(dataset_id, feature, row)
                    manifest_rows.append(
                        {
                            "dataset_id": dataset_id,
                            "dataset_label": cfg["label"],
                            "region_key": row["region_key"],
                            "downstream_status": row["downstream_status"],
                            "truth_status": row["truth_status"],
                            "feature": feature,
                            "feature_rank": rank,
                            "feature_direction_note": note,
                            "feature_value": row[feature],
                            **provenance_payload(row),
                            "agreement_type": row["agreement_type"],
                            "class_shift": row["class_shift"],
                            "af": row["af"],
                            "gq": row["gq"],
                            "qual": row["qual"],
                        }
                    )
                    entry = selected_rows.setdefault(
                        selected_key,
                        {
                            "dataset_id": dataset_id,
                            "dataset_label": cfg["label"],
                            "region_key": row["region_key"],
                            "downstream_status": row["downstream_status"],
                            "truth_status": row["truth_status"],
                            **provenance_payload(row),
                            "agreement_type": row["agreement_type"],
                            "class_shift": row["class_shift"],
                            "Quality_Score": row["Quality_Score"],
                            "PairwiseMedianDist": row["PairwiseMedianDist"],
                            "hp_assign_rate": row["hp_assign_rate"],
                            "af": row["af"],
                            "gq": row["gq"],
                            "qual": row["qual"],
                            "selected_features": [],
                        },
                    )
                    entry["selected_features"].append(f"{feature}:{note}:rank{rank}")

    manifest_fields = [
        "dataset_id",
        "dataset_label",
        "region_key",
        "downstream_status",
        "truth_status",
        "feature",
        "feature_rank",
        "feature_direction_note",
        "feature_value",
        *VERIFICATION_PROVENANCE_COLUMNS,
        "agreement_type",
        "class_shift",
        "af",
        "gq",
        "qual",
    ]
    write_tsv(output_root / "candidate_manifest.tsv", manifest_fields, manifest_rows)

    selected_fields = [
        "dataset_id",
        "dataset_label",
        "region_key",
        "downstream_status",
        "truth_status",
        *VERIFICATION_PROVENANCE_COLUMNS,
        "agreement_type",
        "class_shift",
        "Quality_Score",
        "PairwiseMedianDist",
        "hp_assign_rate",
        "af",
        "gq",
        "qual",
        "selected_features",
    ]
    dedup_rows = []
    for row in selected_rows.values():
        row["selected_features"] = ",".join(sorted(set(row["selected_features"])))
        dedup_rows.append(row)
    write_tsv(output_root / "selected_regions.tsv", selected_fields, dedup_rows)

    snapshot_roots: List[Path] = []
    summary_rows: List[Dict[str, object]] = []

    for dataset_id in dataset_ids:
        cfg = DATASET_CONFIG[dataset_id]
        dataset_dir = output_root / dataset_id
        dataset_dir.mkdir(parents=True, exist_ok=True)
        snapshot_roots.append(dataset_dir)
        tp_index = build_region_index(cfg["tp_root"])
        fp_index = build_region_index(cfg["fp_root"])

        rows = [row for row in dedup_rows if row["dataset_id"] == dataset_id]
        for row in rows:
            key = row["region_key"]
            index = tp_index if row["downstream_status"] == "caller_lost_tp" else fp_index
            if key not in index:
                continue
            info = index[key]
            region_safe = normalize_region_key(key)
            region_out = dataset_dir / "diagnostics" / row["downstream_status"] / region_safe
            region_out.mkdir(parents=True, exist_ok=True)
            run_command(
                [
                    "python3",
                    str(EXPORT_DIAG),
                    "--region-dir",
                    info["region_dir"],
                    "--output-dir",
                    str(region_out / "matrix_diagnostics"),
                ]
            )
            run_command(
                [
                    "bash",
                    str(SAMTOOLS_SNAPSHOT),
                    "--bam",
                    str(cfg["bam"]),
                    "--region",
                    info["region"],
                    "--output-dir",
                    str(region_out / "samtools_snapshot"),
                    "--max-reads",
                    str(args.max_reads),
                ]
            )
            row["region_dir"] = info["region_dir"]
            row["region"] = info["region"]
            row["matrix_notes"] = str(region_out / "matrix_diagnostics" / "region_notes.md")
            row["snapshot_dir"] = str(region_out / "samtools_snapshot")

        run_command(
            [
                "python3",
                str(SUMMARIZE_SNAPSHOTS),
                "--snapshot-root",
                str(dataset_dir),
                "--output-dir",
                str(dataset_dir / "snapshot_summary"),
            ]
        )

        snapshot_rows = read_tsv(dataset_dir / "snapshot_summary" / "snapshot_summary.tsv")
        by_region_id = {row["region_id"]: row for row in snapshot_rows}

        for row in rows:
            region_id = normalize_region_key(row["region_key"])
            snap = by_region_id.get(region_id, {})
            summary_rows.append(
                {
                    **row,
                    "target_alt_fraction": snap.get("target_alt_fraction", ""),
                    "top_hp_tag": snap.get("top_hp_tag", ""),
                    "top_hp_fraction": snap.get("top_hp_fraction", ""),
                    "collapsed_hp_balance_delta": snap.get("collapsed_hp_balance_delta", ""),
                    "na_hp_fraction": snap.get("na_hp_fraction", ""),
                    "snapshot_notes": snap.get("notes", ""),
                }
            )

    summary_fields = [
        "dataset_id",
        "dataset_label",
        "region_key",
        "downstream_status",
        "truth_status",
        *VERIFICATION_PROVENANCE_COLUMNS,
        "agreement_type",
        "class_shift",
        "Quality_Score",
        "PairwiseMedianDist",
        "hp_assign_rate",
        "af",
        "gq",
        "qual",
        "selected_features",
        "region",
        "region_dir",
        "matrix_notes",
        "snapshot_dir",
        "target_alt_fraction",
        "top_hp_tag",
        "top_hp_fraction",
        "collapsed_hp_balance_delta",
        "na_hp_fraction",
        "snapshot_notes",
    ]
    write_tsv(output_root / "diagnostic_summary.tsv", summary_fields, summary_rows)

    print(f"[run_to_support_feature_diagnostics] Wrote {output_root / 'candidate_manifest.tsv'}")
    print(f"[run_to_support_feature_diagnostics] Wrote {output_root / 'selected_regions.tsv'}")
    print(f"[run_to_support_feature_diagnostics] Wrote {output_root / 'diagnostic_summary.tsv'}")


if __name__ == "__main__":
    main()
