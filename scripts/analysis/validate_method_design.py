#!/usr/bin/env python3
"""Validate label-first and cluster-first design usability and agreement."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from research_common import as_bool, to_float, to_int, write_tsv_rows

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (
    SchemaContractError,
    read_evidence,
    select_current_view,
    select_legacy_view,
)


DESIGN_FIELDS = [
    "region_key",
    "region_id",
    "source_scope",
    "design_mode",
    "input_label_type",
    "cluster_count",
    "hp_assign_rate",
    "allele_assign_rate",
    "test_applicable",
    "failure_reason",
    "output_class",
    "interpretability_grade",
    "consistency_grade",
    "cluster_class",
    "cluster_class_source",
    "cluster_class_selection_field",
    "verification_schema_status",
    "dominant_label",
    "notes",
]

AGREEMENT_FIELDS = [
    "region_key",
    "region_id",
    "source_scope",
    "reads",
    "hp_assign_rate",
    "allele_assign_rate",
    "cluster_class",
    "cluster_class_source",
    "cluster_class_selection_field",
    "verification_schema_status",
    "label_class",
    "agreement_type",
    "cluster_permanova_p",
    "label_hp_permanova_p",
    "label_allele_permanova_p",
    "fisher_p",
    "effect_size",
    "class_shift",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-csv", action="append", required=True, help="Input significance_summary.csv")
    parser.add_argument("--region-root", action="append", default=[], help="Optional InterSubMod root for region reads.tsv")
    parser.add_argument("--haplotag-qc", default=None, help="Optional haplotag_qc.tsv")
    parser.add_argument("--sample", default="unknown", help="Sample name")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument(
        "--cluster-class-source",
        required=True,
        choices=("legacy", "evidence"),
        help=(
            "Explicit four-state comparison source: VerificationClass_Legacy or a "
            "four-state view derived from typed LabelFirstSupport/ClusterFirstSupport"
        ),
    )
    parser.add_argument("--min-reads", type=int, default=20, help="Minimum reads required for method validation")
    parser.add_argument(
        "--min-hp-assign-rate",
        type=float,
        default=0.05,
        help="Minimum HP assign rate to consider HP label-first usable",
    )
    return parser.parse_args()


def load_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def select_cluster_class_frame(frame: pd.DataFrame, source: str) -> pd.DataFrame:
    """Select the four-state comparison input without reading current names as legacy."""
    current = select_current_view(frame)
    if current.unknown_counts:
        raise SchemaContractError(
            "method design validation: unknown current VerificationClass values: "
            f"{current.unknown_counts}"
        )

    prepared = frame.copy()
    if source == "legacy":
        legacy = select_legacy_view(frame)
        prepared["_cluster_class_selected"] = legacy.values
        selection_field = legacy.field
        schema_status = f"{current.schema_status}+{legacy.schema_status}"
    elif source == "evidence":
        evidence = read_evidence(frame)
        label_support = evidence["LabelFirstSupport"]
        cluster_support = evidence["ClusterFirstSupport"]
        prepared["_cluster_class_selected"] = "Noise"
        prepared.loc[label_support & ~cluster_support, "_cluster_class_selected"] = "Weak"
        prepared.loc[~label_support & cluster_support, "_cluster_class_selected"] = "Subclone"
        prepared.loc[label_support & cluster_support, "_cluster_class_selected"] = "Strong"
        selection_field = "LabelFirstSupport+ClusterFirstSupport"
        schema_status = f"{current.schema_status}+TYPED_EVIDENCE"
    else:
        raise SchemaContractError(f"unsupported cluster class source: {source!r}")

    prepared["_cluster_class_source"] = source
    prepared["_cluster_class_selection_field"] = selection_field
    prepared["_verification_schema_status"] = schema_status
    return prepared


def load_selected_summary_rows(path: Path, source: str) -> List[Dict[str, object]]:
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    return select_cluster_class_frame(frame, source).to_dict(orient="records")


def parse_metadata(path: Path) -> Optional[Dict[str, str]]:
    region_dir = path.parent
    region_id = ""
    chrom = ""
    pos = ""
    ref = ""
    alt = ""
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith("Region ID:"):
                region_id = line.split(":", 1)[1].strip()
            elif line.startswith("SNV Position:"):
                chrom_pos = line.split(":", 1)[1].strip()
                if ":" in chrom_pos:
                    chrom, pos = chrom_pos.split(":", 1)
            elif line.startswith("SNV:"):
                parts = line.split(":", 1)[1].strip().split("->")
                if len(parts) == 2:
                    ref = parts[0].strip()
                    alt = parts[1].strip()
    if not chrom or not pos:
        return None
    return {
        "region_key": f"{chrom}:{pos}:{ref}:{alt}",
        "region_id": region_id,
        "region_dir": str(region_dir),
    }


def build_region_index(roots: List[Path]) -> Dict[str, Dict[str, str]]:
    index: Dict[str, Dict[str, str]] = {}
    for root in roots:
        if not root.exists():
            continue
        for metadata in root.rglob("metadata.txt"):
            entry = parse_metadata(metadata)
            if entry:
                index[entry["region_key"]] = entry
    return index


def compute_region_read_stats(region_dir: Path) -> Dict[str, float]:
    reads_path = region_dir / "reads" / "reads.tsv"
    stats = {
        "reads_total": 0,
        "hp_assign_rate": math.nan,
        "allele_assign_rate": math.nan,
    }
    if not reads_path.exists():
        return stats

    hp_assigned = 0
    allele_assigned = 0
    with reads_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    total = len(rows)
    if total == 0:
        return stats
    for row in rows:
        hp = str(row.get("hp", "")).strip()
        allele = str(row.get("alt_support", "")).strip().upper()
        if hp and hp not in {"0", "NA"}:
            hp_assigned += 1
        if allele and allele != "UNKNOWN":
            allele_assigned += 1
    stats["reads_total"] = total
    stats["hp_assign_rate"] = hp_assigned / total
    stats["allele_assign_rate"] = allele_assigned / total
    return stats


def infer_summary_read_stats(row: Dict[str, str]) -> Dict[str, float]:
    num_reads = to_int(row.get("NumReads"))
    hp1 = to_int(row.get("HP1FamilyN"))
    hp2 = to_int(row.get("HP2FamilyN"))
    stats = {
        "reads_total": num_reads,
        "hp_assign_rate": math.nan,
        "allele_assign_rate": math.nan,
    }
    if num_reads > 0 and (hp1 > 0 or hp2 > 0):
        assigned = min(max(hp1 + hp2, 0), num_reads)
        stats["hp_assign_rate"] = assigned / num_reads
    return stats


def read_sample_hp_rate(path: Optional[Path]) -> float:
    if path is None or not path.exists():
        return math.nan
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        return math.nan
    return to_float(rows[0].get("hp_assign_rate"))


def infer_hp_class(row: Dict[str, str]) -> str:
    hp_sig = (
        as_bool(row.get("HPMergedSig"))
        or as_bool(row.get("HPFineSig"))
        or (as_bool(row.get("LabelHPPermanovaValid")) and to_float(row.get("LabelHPPermanovaP")) < 0.05)
    )
    if hp_sig:
        return "Subclone"
    return "Noise"


def infer_allele_class(row: Dict[str, str]) -> str:
    allele_sig = as_bool(row.get("AlleleSig")) or (
        as_bool(row.get("LabelAllelePermanovaValid")) and to_float(row.get("LabelAllelePermanovaP")) < 0.05
    )
    if allele_sig:
        return "Weak"
    return "Noise"


def infer_combined_label_class(hp_class: str, allele_class: str) -> str:
    if hp_class == "Subclone" and allele_class == "Weak":
        return "Strong"
    if hp_class == "Subclone":
        return "Subclone"
    if allele_class == "Weak":
        return "Weak"
    return "Noise"


def interpretability_grade(output_class: str, applicable: bool, dominant_label: str) -> str:
    if not applicable:
        return "D"
    if output_class in {"Strong", "Subclone"} and dominant_label in {"hp", "allele"}:
        return "A"
    if output_class in {"Strong", "Subclone", "Weak"}:
        return "B"
    return "C"


def consistency_grade(cluster_class: str, design_class: str, applicable: bool) -> str:
    if not applicable:
        return "insufficient"
    if cluster_class == design_class:
        return "high"
    if cluster_class in {"Strong", "Subclone"} and design_class in {"Strong", "Subclone"}:
        return "medium"
    if cluster_class in {"Weak", "Noise"} and design_class in {"Weak", "Noise"}:
        return "medium"
    return "low"


def agreement_type(cluster_class: str, label_class: str, applicable: bool) -> str:
    if not applicable:
        return "insufficient_label"
    if cluster_class == "Strong" and label_class == "Strong":
        return "consistent_strong"
    if cluster_class == "Subclone" and label_class == "Subclone":
        return "consistent_subclone"
    if cluster_class in {"Weak", "Noise"} and label_class in {"Strong", "Subclone"}:
        return "label_upgrade"
    if cluster_class in {"Strong", "Subclone"} and label_class == "Weak":
        return "cluster_plus_weak_label"
    if cluster_class in {"Strong", "Subclone"} and label_class in {"Weak", "Noise"}:
        return "cluster_only"
    if cluster_class == label_class:
        return "consistent_weak_or_noise"
    return "conflict"


def format_ratio(value: float) -> str:
    if math.isnan(value):
        return ""
    return f"{value:.6f}"


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    region_roots = [Path(root) for root in args.region_root]
    region_index: Optional[Dict[str, Dict[str, str]]] = None
    region_stats_cache: Dict[str, Dict[str, float]] = {}
    sample_hp_rate = read_sample_hp_rate(Path(args.haplotag_qc) if args.haplotag_qc else None)

    design_rows: List[Dict[str, str]] = []
    agreement_rows: List[Dict[str, str]] = []

    for summary_csv in args.summary_csv:
        summary_path = Path(summary_csv)
        scope = summary_path.parent.name.replace("intersubmod_", "")
        for row in load_selected_summary_rows(summary_path, args.cluster_class_source):
            region_key = f"{row['Chr']}:{row['Pos']}:{row['Ref']}:{row['Alt']}"
            if region_key not in region_stats_cache:
                region_stats_cache[region_key] = infer_summary_read_stats(row)
            region_stats = region_stats_cache.get(region_key, {})
            region_dir = None
            if region_dir and math.isnan(region_stats.get("hp_assign_rate", math.nan)) and math.isnan(region_stats.get("allele_assign_rate", math.nan)):
                region_stats_cache[region_key] = compute_region_read_stats(region_dir)
                region_stats = region_stats_cache.get(region_key, {})
            elif math.isnan(region_stats.get("hp_assign_rate", math.nan)) and math.isnan(region_stats.get("allele_assign_rate", math.nan)):
                if region_index is None:
                    region_index = build_region_index(region_roots)
                region_info = region_index.get(region_key, {})
                region_dir = Path(region_info["region_dir"]) if "region_dir" in region_info else None
                if region_dir:
                    region_stats_cache[region_key] = compute_region_read_stats(region_dir)
                    region_stats = region_stats_cache.get(region_key, {})

            hp_assign_rate = region_stats.get("hp_assign_rate", sample_hp_rate)
            allele_assign_rate = region_stats.get("allele_assign_rate", math.nan)
            num_reads = to_int(row.get("NumReads"))
            cluster_class = row["_cluster_class_selected"]
            selection_metadata = {
                "cluster_class_source": row["_cluster_class_source"],
                "cluster_class_selection_field": row[
                    "_cluster_class_selection_field"
                ],
                "verification_schema_status": row["_verification_schema_status"],
            }
            dominant_label = row.get("DominantLabel", "none") or "none"
            cluster_count = row.get("LocalBestCluster", "")

            hp_applicable = num_reads >= args.min_reads and not math.isnan(hp_assign_rate) and hp_assign_rate >= args.min_hp_assign_rate
            allele_applicable = num_reads >= args.min_reads and (
                (not math.isnan(allele_assign_rate) and allele_assign_rate > 0.0)
                or to_float(row.get("AlleleDelta"), 0.0) != 0.0
                or as_bool(row.get("AlleleSig"))
            )

            hp_failure = ""
            if num_reads < args.min_reads:
                hp_failure = "insufficient_reads"
            elif math.isnan(hp_assign_rate):
                hp_failure = "missing_hp_rate"
            elif hp_assign_rate < args.min_hp_assign_rate:
                hp_failure = "low_hp_assign_rate"

            allele_failure = ""
            if num_reads < args.min_reads:
                allele_failure = "insufficient_reads"
            elif math.isnan(allele_assign_rate) and to_float(row.get("AlleleDelta"), 0.0) == 0.0 and not as_bool(row.get("AlleleSig")):
                allele_failure = "missing_allele_signal"

            hp_class = infer_hp_class(row) if hp_applicable else "NA"
            allele_class = infer_allele_class(row) if allele_applicable else "NA"
            combined_label_class = infer_combined_label_class(
                hp_class if hp_class != "NA" else "Noise",
                allele_class if allele_class != "NA" else "Noise",
            )
            combined_applicable = hp_applicable or allele_applicable

            design_rows.append(
                {
                    "region_key": region_key,
                    "region_id": row.get("RegionID", ""),
                    "source_scope": scope,
                    "design_mode": "cluster_first",
                    "input_label_type": "cluster",
                    "cluster_count": cluster_count,
                    "hp_assign_rate": format_ratio(hp_assign_rate),
                    "allele_assign_rate": format_ratio(allele_assign_rate),
                    "test_applicable": "true" if num_reads >= args.min_reads else "false",
                    "failure_reason": "" if num_reads >= args.min_reads else "insufficient_reads",
                    "output_class": cluster_class,
                    "interpretability_grade": interpretability_grade(cluster_class, num_reads >= args.min_reads, dominant_label),
                    "consistency_grade": "self",
                    "cluster_class": cluster_class,
                    **selection_metadata,
                    "dominant_label": dominant_label,
                    "notes": f"sample={args.sample}",
                }
            )

            design_rows.append(
                {
                    "region_key": region_key,
                    "region_id": row.get("RegionID", ""),
                    "source_scope": scope,
                    "design_mode": "label_first_hp",
                    "input_label_type": "hp",
                    "cluster_count": cluster_count,
                    "hp_assign_rate": format_ratio(hp_assign_rate),
                    "allele_assign_rate": format_ratio(allele_assign_rate),
                    "test_applicable": "true" if hp_applicable else "false",
                    "failure_reason": hp_failure,
                    "output_class": hp_class,
                    "interpretability_grade": interpretability_grade(hp_class, hp_applicable, dominant_label),
                    "consistency_grade": consistency_grade(cluster_class, hp_class, hp_applicable),
                    "cluster_class": cluster_class,
                    **selection_metadata,
                    "dominant_label": dominant_label,
                    "notes": f"HPMergedSig={row.get('HPMergedSig', '')}; HPFineSig={row.get('HPFineSig', '')}",
                }
            )

            design_rows.append(
                {
                    "region_key": region_key,
                    "region_id": row.get("RegionID", ""),
                    "source_scope": scope,
                    "design_mode": "label_first_allele",
                    "input_label_type": "allele",
                    "cluster_count": cluster_count,
                    "hp_assign_rate": format_ratio(hp_assign_rate),
                    "allele_assign_rate": format_ratio(allele_assign_rate),
                    "test_applicable": "true" if allele_applicable else "false",
                    "failure_reason": allele_failure,
                    "output_class": allele_class,
                    "interpretability_grade": interpretability_grade(allele_class, allele_applicable, dominant_label),
                    "consistency_grade": consistency_grade(cluster_class, allele_class, allele_applicable),
                    "cluster_class": cluster_class,
                    **selection_metadata,
                    "dominant_label": dominant_label,
                    "notes": f"AlleleSig={row.get('AlleleSig', '')}",
                }
            )

            effect_size = max(
                abs(to_float(row.get("CramersV"), 0.0)),
                abs(to_float(row.get("HPMergedDelta"), 0.0)),
                abs(to_float(row.get("AlleleDelta"), 0.0)),
            )
            agreement_rows.append(
                {
                    "region_key": region_key,
                    "region_id": row.get("RegionID", ""),
                    "source_scope": scope,
                    "reads": str(num_reads),
                    "hp_assign_rate": format_ratio(hp_assign_rate),
                    "allele_assign_rate": format_ratio(allele_assign_rate),
                    "cluster_class": cluster_class,
                    **selection_metadata,
                    "label_class": combined_label_class,
                    "agreement_type": agreement_type(cluster_class, combined_label_class, combined_applicable),
                    "cluster_permanova_p": row.get("ClusterPermanovaP", ""),
                    "label_hp_permanova_p": row.get("LabelHPPermanovaP", ""),
                    "label_allele_permanova_p": row.get("LabelAllelePermanovaP", ""),
                    "fisher_p": row.get("GlobalP", ""),
                    "effect_size": f"{effect_size:.6f}",
                    "class_shift": "unchanged" if cluster_class == combined_label_class else f"{cluster_class}->{combined_label_class}",
                }
            )

    write_tsv_rows(output_dir / "method_design_validation.tsv", DESIGN_FIELDS, design_rows)
    write_tsv_rows(output_dir / "label_cluster_agreement.tsv", AGREEMENT_FIELDS, agreement_rows)

    print(f"[validate_method_design] Wrote {output_dir / 'method_design_validation.tsv'}")
    print(f"[validate_method_design] Wrote {output_dir / 'label_cluster_agreement.tsv'}")


if __name__ == "__main__":
    main()
