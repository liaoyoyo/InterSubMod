#!/usr/bin/env python3
"""Build the Round 1 LOH cross-sample audit workspace."""

from __future__ import annotations

import argparse
import bisect
import json
import math
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

from research_common import ensure_dir, markdown_table, read_json, write_json, write_tsv_rows

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

PROVENANCE_METADATA_FIELDS = (
    "VerificationProvenanceStatus",
    "VerificationProvenanceSourceField",
    "VerificationProvenanceWarnings",
    "VerificationUnknownCurrentCounts",
    "LOHProvenanceSourceField",
)


OUTPUT_ROOT_DEFAULT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260327_loh_round1_cross_sample_audit"
)


@dataclass(frozen=True)
class DatasetConfig:
    dataset_id: str
    dataset_label: str
    sample: str
    sample_label: str
    platform: str
    mode: str
    mode_label: str
    source_kind: str
    source_priority: str
    base_dir: Path
    tp_summary: Path
    fp_summary: Path
    tp_region_root: Path
    fp_region_root: Path
    tp_vcf: Path
    fp_vcf: Path
    benchmark_tsv: Path
    context_json: Path
    tagged_bam: Path
    phased_vcf: Optional[Path] = None
    loh_bed: Optional[Path] = None


def paired_config(sample: str, sample_label: str, platform: str, run_id: str) -> DatasetConfig:
    base = Path("/big7_disk/liaoyoyo2001/big7_disk_output/canonical") / sample / "paired_full" / run_id
    bam_name = f"{sample}_tagged.bam"
    return DatasetConfig(
        dataset_id=f"{sample_label.lower()}_paired".replace("-", "_"),
        dataset_label=f"{sample_label} paired_full",
        sample=sample,
        sample_label=sample_label,
        platform=platform,
        mode="paired",
        mode_label="paired-pure",
        source_kind="canonical",
        source_priority="canonical_primary",
        base_dir=base,
        tp_summary=base / "intersubmod_tp" / "significance_summary.csv",
        fp_summary=base / "intersubmod_fp" / "significance_summary.csv",
        tp_region_root=base / "intersubmod_tp" / "filtered_snv_tp",
        fp_region_root=base / "intersubmod_fp" / "filtered_snv_fp",
        tp_vcf=base / "longphase_s" / "filtered_snv_tp.vcf.gz",
        fp_vcf=base / "longphase_s" / "filtered_snv_fp.vcf.gz",
        benchmark_tsv=base / "benchmark_comparison.tsv",
        context_json=base / "run_context.json",
        tagged_bam=base / "longphase_s" / bam_name,
    )


def to_config(sample: str, sample_label: str, platform: str, round_id: str) -> DatasetConfig:
    base = Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds") / round_id
    return DatasetConfig(
        dataset_id=f"{sample_label.lower()}_to".replace("-", "_"),
        dataset_label=f"{sample_label} tumor_only",
        sample=sample,
        sample_label=sample_label,
        platform=platform,
        mode="to",
        mode_label="to-pure",
        source_kind="research_round",
        source_priority="research_round_primary",
        base_dir=base,
        tp_summary=base / "step05_intersubmod" / "intersubmod_tp" / "significance_summary.csv",
        fp_summary=base / "step05_intersubmod" / "intersubmod_fp" / "significance_summary.csv",
        tp_region_root=base / "step05_intersubmod" / "intersubmod_tp" / "filtered_snv_tp",
        fp_region_root=base / "step05_intersubmod" / "intersubmod_fp" / "filtered_snv_fp",
        tp_vcf=base / "step04_benchmark_longphase_to" / "filtered_snv_tp.vcf.gz",
        fp_vcf=base / "step04_benchmark_longphase_to" / "filtered_snv_fp.vcf.gz",
        benchmark_tsv=base / "benchmark_comparison.tsv",
        context_json=base / "round_context.json",
        tagged_bam=base / "step03_longphase_to" / "tumor_tagged.bam",
        phased_vcf=base / "step03_longphase_to" / "tumor_phased.vcf",
        loh_bed=base / "step03_longphase_to" / "tumor_phased_LOH.bed",
    )


DATASET_CONFIGS: Sequence[DatasetConfig] = [
    paired_config("HCC1395", "HCC1395_HKU_5kHz", "ONT_5kHz", "20260314_HCC1395_paired_full_full_complete_matrix"),
    to_config("HCC1395", "HCC1395_HKU_5kHz", "ONT_5kHz", "20260315_hcc1395_to_pilot"),
    paired_config("HCC1395_DORADO", "HCC1395_DORADO", "ONT_Dorado", "20260315_HCC1395_DORADO_paired_full_full_complete_matrix"),
    to_config("HCC1395_DORADO", "HCC1395_DORADO", "ONT_Dorado", "20260315_hcc1395_dorado_to_pilot"),
    paired_config("COLO829", "COLO829", "ONT", "20260315_COLO829_paired_full_full_complete_matrix"),
    to_config("COLO829", "COLO829", "ONT", "20260317_colo829_to_pilot_1"),
    paired_config("H1437", "H1437", "ONT", "20260315_H1437_paired_full_full_complete_matrix"),
    to_config("H1437", "H1437", "ONT", "20260318_h1437_to_pilot_fastresume"),
    paired_config("H2009", "H2009", "ONT", "20260315_H2009_paired_full_full_complete_matrix"),
    to_config("H2009", "H2009", "ONT", "20260318_h2009_to_pilot_fastresume"),
    paired_config("HCC1937", "HCC1937", "ONT", "20260315_HCC1937_paired_full_full_complete_matrix"),
    to_config("HCC1937", "HCC1937", "ONT", "20260318_hcc1937_to_pilot_fastresume"),
    paired_config("HCC1954", "HCC1954", "ONT", "20260315_HCC1954_paired_full_full_complete_matrix"),
    to_config("HCC1954", "HCC1954", "ONT", "20260318_hcc1954_to_pilot_fastresume"),
]


MANIFEST_FIELDS = [
    "dataset_id",
    "dataset_label",
    "sample",
    "sample_label",
    "platform",
    "mode",
    "mode_label",
    "source_kind",
    "source_priority",
    "availability_status",
    "missing_components",
    "base_dir",
    "tp_summary",
    "fp_summary",
    "tp_vcf",
    "fp_vcf",
    "benchmark_tsv",
    "context_json",
    "tagged_bam",
    "phased_vcf",
    "loh_bed",
    "notes",
]

DECISION_FIELDS = [
    "decision_id",
    "decision_scope_round1",
    "decision_topic",
    "decision",
    "status",
    "rationale",
    "evidence_table",
    "evidence_figure",
    "generated_at",
]

CASE_REGISTRY_FIELDS = [
    "case_id",
    "case_category",
    "sample",
    "sample_label",
    "variant_key",
    "chrom",
    "pos",
    "ref",
    "alt",
    "locus_overlap_class",
    "paired_status",
    "to_status",
    "paired_verification_class",
    "to_verification_class",
    "paired_core_loh_like",
    "to_core_loh_like",
    "paired_hp_ratio_core",
    "to_hp_ratio_core",
    "paired_effective_hp_reads",
    "to_effective_hp_reads",
    "case_dir",
]

SUMMARY_METADATA_COLS = [
    "sample",
    "sample_label",
    "mode",
    "mode_label",
    "tool_role",
    "source_file",
    "source_record_count",
    "analysis_subset_rule",
    "generated_at",
]

HP_RATIO_BINS = [
    (-math.inf, math.inf, "no_effective_hp"),
    (-math.inf, 0.05, "<0.05"),
    (0.05, 0.10, "0.05-0.10"),
    (0.10, 0.20, "0.10-0.20"),
    (0.20, 0.80, "0.20-0.80"),
    (0.80, 0.90, "0.80-0.90"),
    (0.90, 0.95, "0.90-0.95"),
    (0.95, math.inf, ">0.95"),
]

AF_BINS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01]
AF_BIN_LABELS = [
    "0.0-0.1",
    "0.1-0.2",
    "0.2-0.3",
    "0.3-0.4",
    "0.4-0.5",
    "0.5-0.6",
    "0.6-0.7",
    "0.7-0.8",
    "0.8-0.9",
    ">=0.9",
]
EFFECTIVE_HP_BINS = [0, 1, 5, 10, 20, 40, 80, math.inf]
EFFECTIVE_HP_BIN_LABELS = ["0", "1-4", "5-9", "10-19", "20-39", "40-79", ">=80"]

CASE_CATEGORY_ORDER = [
    "both_tp_loh_like",
    "both_fp_loh_like",
    "paired_tp_to_fp",
    "paired_fp_to_tp",
    "paired_only_loh_like",
    "to_only_loh_like",
    "high_hp0_or_hp3",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(OUTPUT_ROOT_DEFAULT), help="Observation workspace output directory")
    parser.add_argument("--sample", action="append", default=None, help="Optional sample filter")
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help="Explicitly authorize historical H1 summaries without VerificationSchemaVersion.",
    )
    return parser.parse_args()


def generated_at_string() -> str:
    return datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")


def variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{chrom}:{pos}:{ref}:{alt}"


def safe_div(num: float, den: float) -> float:
    if den in (0, 0.0) or den is None or (isinstance(den, float) and math.isnan(den)):
        return math.nan
    return num / den


def first_scalar(value):
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return None
        return first_scalar(value.tolist())
    if isinstance(value, (list, tuple)):
        if not value:
            return None
        return first_scalar(value[0])
    return value


def to_float_or_nan(value) -> float:
    value = first_scalar(value)
    if value is None:
        return math.nan
    try:
        text = str(value).strip()
        if not text or text == ".":
            return math.nan
        return float(text)
    except (TypeError, ValueError):
        return math.nan


def to_int_or_zero(value) -> int:
    value = first_scalar(value)
    if value is None:
        return 0
    try:
        text = str(value).strip()
        if not text or text == ".":
            return 0
        return int(float(text))
    except (TypeError, ValueError):
        return 0


def bool_from_value(value) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def format_metric(value, digits: int = 4, integer: bool = False) -> str:
    if value is None:
        return "NA"
    try:
        number = float(value)
    except (TypeError, ValueError):
        text = str(value).strip()
        return text or "NA"
    if math.isnan(number):
        return "NA"
    if integer:
        return str(int(round(number)))
    return f"{number:.{digits}f}"


def format_gt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        parts = []
        for item in value:
            if item is None:
                parts.append(".")
            else:
                parts.append(str(item))
        return "/".join(parts)
    return str(value)


def normalize_hp_tag(value) -> str:
    if value is None:
        return "HP0"
    text = str(value).strip()
    if not text or text.lower() == "none":
        return "HP0"
    if text == "0":
        return "HP0"
    return f"HP{text}"


def hp_family(tag: str) -> str:
    if tag in {"HP1", "HP1-1"}:
        return "HP1_family"
    if tag in {"HP2", "HP2-1"}:
        return "HP2_family"
    if tag == "HP3":
        return "HP3"
    if tag == "HP0":
        return "HP0"
    return "OTHER"


def hp_ratio_bin(ratio: float, effective_hp_reads: float) -> str:
    if effective_hp_reads <= 0 or math.isnan(ratio):
        return "no_effective_hp"
    if ratio < 0.05:
        return "<0.05"
    if ratio < 0.10:
        return "0.05-0.10"
    if ratio < 0.20:
        return "0.10-0.20"
    if ratio < 0.80:
        return "0.20-0.80"
    if ratio < 0.90:
        return "0.80-0.90"
    if ratio < 0.95:
        return "0.90-0.95"
    return ">0.95"


def af_bin(af: float) -> str:
    if math.isnan(af):
        return "unknown"
    idx = np.digitize([af], AF_BINS, right=False)[0] - 1
    idx = max(0, min(idx, len(AF_BIN_LABELS) - 1))
    return AF_BIN_LABELS[idx]


def effective_hp_bin(value: float) -> str:
    if math.isnan(value):
        return "unknown"
    if value == 0:
        return "0"
    if value < 5:
        return "1-4"
    if value < 10:
        return "5-9"
    if value < 20:
        return "10-19"
    if value < 40:
        return "20-39"
    if value < 80:
        return "40-79"
    return ">=80"


def sample_mode_label(sample_label: str, mode: str) -> str:
    return f"{sample_label}|{mode}"


def parse_filter_terms(record: pysam.VariantRecord) -> str:
    terms = list(record.filter.keys())
    return ";".join(terms or ["PASS"])


def parse_ad_field(value) -> tuple[int, int]:
    if value is None:
        return 0, 0
    if isinstance(value, (list, tuple)):
        parts = list(value)
    else:
        parts = str(value).split(",")
    if len(parts) < 2:
        return to_int_or_zero(parts[0] if parts else 0), 0
    return to_int_or_zero(parts[0]), to_int_or_zero(parts[1])


class BedIndex:
    def __init__(self, path: Optional[Path]) -> None:
        self.path = path
        self.by_chrom: Dict[str, tuple[List[int], List[int]]] = {}
        if path is None or not path.exists():
            return
        starts_by_chrom: Dict[str, List[int]] = {}
        ends_by_chrom: Dict[str, List[int]] = {}
        with path.open("r", encoding="utf-8") as handle:
            for raw in handle:
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                chrom = fields[0]
                starts_by_chrom.setdefault(chrom, []).append(int(fields[1]))
                ends_by_chrom.setdefault(chrom, []).append(int(fields[2]))
        for chrom in starts_by_chrom:
            pairs = sorted(zip(starts_by_chrom[chrom], ends_by_chrom[chrom]))
            self.by_chrom[chrom] = ([x[0] for x in pairs], [x[1] for x in pairs])

    def contains(self, chrom: str, pos_1based: int) -> bool:
        if chrom not in self.by_chrom:
            return False
        starts, ends = self.by_chrom[chrom]
        query = pos_1based - 1
        idx = bisect.bisect_right(starts, query) - 1
        if idx < 0:
            return False
        return query < ends[idx]


def load_split_vcf_features(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["variant_key"])
    rows: List[Dict] = []
    vf = pysam.VariantFile(str(path))
    for record in vf.fetch():
        if len(record.ref) != 1 or not record.alts or len(record.alts) != 1 or len(record.alts[0]) != 1:
            continue
        sample_data = next(iter(record.samples.values()), None)
        af = math.nan
        dp = 0
        gq = 0
        ad_ref = 0
        ad_alt = 0
        naf = math.nan
        ndp = 0
        nad_ref = 0
        nad_alt = 0
        if sample_data is not None:
            af = to_float_or_nan(sample_data.get("AF"))
            dp = to_int_or_zero(sample_data.get("DP"))
            gq = to_int_or_zero(sample_data.get("GQ"))
            ad_ref, ad_alt = parse_ad_field(sample_data.get("AD"))
            naf = to_float_or_nan(sample_data.get("NAF"))
            ndp = to_int_or_zero(sample_data.get("NDP"))
            nad_ref, nad_alt = parse_ad_field(sample_data.get("NAD"))
        rows.append(
            {
                "variant_key": variant_key(record.chrom, record.pos, record.ref, record.alts[0]),
                "caller_filter": parse_filter_terms(record),
                "caller_af": af,
                "caller_dp": dp,
                "caller_gq": gq,
                "caller_ad_ref": ad_ref,
                "caller_ad_alt": ad_alt,
                "caller_naf": naf,
                "caller_ndp": ndp,
                "caller_nad_ref": nad_ref,
                "caller_nad_alt": nad_alt,
            }
        )
    vf.close()
    return pd.DataFrame.from_records(rows)


def load_to_phased_features(path: Optional[Path]) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame(columns=["variant_key"])
    rows: List[Dict] = []
    vf = pysam.VariantFile(str(path))
    for record in vf.fetch():
        if len(record.ref) != 1 or not record.alts or len(record.alts) != 1 or len(record.alts[0]) != 1:
            continue
        sample_data = next(iter(record.samples.values()), None)
        if sample_data is None:
            continue
        rows.append(
            {
                "variant_key": variant_key(record.chrom, record.pos, record.ref, record.alts[0]),
                "to_phase_filter": parse_filter_terms(record),
                "to_phase_ps": sample_data.get("PS") if sample_data.get("PS") is not None else "",
                "to_phase_gt": format_gt(sample_data.get("GT")),
                "to_phase_gt2": format_gt(sample_data.get("GT2")),
                "to_phase_ps2": sample_data.get("PS2") if sample_data.get("PS2") is not None else "",
                "to_phase_gt3": format_gt(sample_data.get("GT3")),
                "to_phase_ps3": sample_data.get("PS3") if sample_data.get("PS3") is not None else "",
            }
        )
    vf.close()
    return pd.DataFrame.from_records(rows)


def availability_status(config: DatasetConfig) -> Dict[str, str]:
    required = {
        "tp_summary": config.tp_summary,
        "fp_summary": config.fp_summary,
        "tp_vcf": config.tp_vcf,
        "fp_vcf": config.fp_vcf,
        "benchmark_tsv": config.benchmark_tsv,
        "context_json": config.context_json,
        "tagged_bam": config.tagged_bam,
    }
    if config.mode == "to":
        required["phased_vcf"] = config.phased_vcf
        required["loh_bed"] = config.loh_bed
    missing = [name for name, path in required.items() if path is None or not Path(path).exists()]
    status = "ready" if not missing else "missing"
    return {
        "availability_status": status,
        "missing_components": ",".join(missing),
    }


def build_manifest_rows(configs: Sequence[DatasetConfig]) -> List[Dict]:
    rows: List[Dict] = []
    for config in configs:
        availability = availability_status(config)
        notes = ""
        if config.mode == "to":
            notes = "TO uses research_round pilot/fastresume outputs; PS only available via tumor_phased.vcf"
        elif config.mode == "paired":
            notes = "paired uses canonical paired_full complete_matrix outputs; PS only available from tagged BAM reads"
        row = {
            "dataset_id": config.dataset_id,
            "dataset_label": config.dataset_label,
            "sample": config.sample,
            "sample_label": config.sample_label,
            "platform": config.platform,
            "mode": config.mode,
            "mode_label": config.mode_label,
            "source_kind": config.source_kind,
            "source_priority": config.source_priority,
            "availability_status": availability["availability_status"],
            "missing_components": availability["missing_components"],
            "base_dir": str(config.base_dir),
            "tp_summary": str(config.tp_summary),
            "fp_summary": str(config.fp_summary),
            "tp_vcf": str(config.tp_vcf),
            "fp_vcf": str(config.fp_vcf),
            "benchmark_tsv": str(config.benchmark_tsv),
            "context_json": str(config.context_json),
            "tagged_bam": str(config.tagged_bam),
            "phased_vcf": str(config.phased_vcf) if config.phased_vcf else "",
            "loh_bed": str(config.loh_bed) if config.loh_bed else "",
            "notes": notes,
        }
        rows.append(row)
    return rows


def attach_verification_contract(
    df: pd.DataFrame,
    allow_unversioned_v1: bool = False,
) -> pd.DataFrame:
    """Attach normalized C2/LOH-L views while preserving every raw provenance field."""

    selected = df.copy()
    if "VerificationSchemaVersion" in selected.columns:
        current = select_current_view(selected)
        select_legacy_view(selected)
        read_evidence(selected)
        loh = select_loh_legacy(selected)
        missing = [field for field in VERIFICATION_PROVENANCE_COLUMNS if field not in selected.columns]
        if missing:
            raise SchemaContractError(
                "LOH round provenance: missing required fields: " + ", ".join(missing)
            )
        known = current.values != UNKNOWN_CURRENT_CLASS
        if known.any():
            extract_provenance_frame(selected.loc[known])
        status = "V2"
        warning_payload = {"unknown_current_counts": current.unknown_counts}
    else:
        if not allow_unversioned_v1:
            raise SchemaContractError(
                "LOH round provenance: VerificationSchemaVersion missing; explicit "
                "--allow-unversioned-v1 authorization is required for historical H1 input"
            )
        current = select_legacy_view(selected, allow_unversioned_v1=True)
        loh = select_loh_legacy(selected, allow_unversioned_v1=True)
        provenance = extract_provenance_frame(selected, allow_unversioned_v1=True)
        status = "UNVERSIONED_V1"
        warning_payload = {
            "authorization": "--allow-unversioned-v1",
            "messages": list(provenance.attrs.get("warnings", [])),
        }
        for field in VERIFICATION_PROVENANCE_COLUMNS:
            if field not in selected.columns:
                selected[field] = ""

    selected["verification_class"] = current.values
    selected["loh_subtype"] = loh.values
    selected["VerificationProvenanceStatus"] = status
    selected["VerificationProvenanceSourceField"] = current.field
    selected["LOHProvenanceSourceField"] = loh.field
    selected["VerificationUnknownCurrentCounts"] = json.dumps(
        current.unknown_counts, ensure_ascii=False, sort_keys=True
    )
    selected["VerificationProvenanceWarnings"] = json.dumps(
        warning_payload, ensure_ascii=False, sort_keys=True
    )
    return selected


def load_dataset_rows(
    config: DatasetConfig,
    allow_unversioned_v1: bool = False,
) -> pd.DataFrame:
    vcf_maps = {
        "TP": load_split_vcf_features(config.tp_vcf),
        "FP": load_split_vcf_features(config.fp_vcf),
    }
    to_phase_df = load_to_phased_features(config.phased_vcf)
    bed_index = BedIndex(config.loh_bed)
    context = read_json(config.context_json) if config.context_json.exists() else {}
    frames: List[pd.DataFrame] = []
    for truth_label, summary_path in [("TP", config.tp_summary), ("FP", config.fp_summary)]:
        df = attach_verification_contract(
            pd.read_csv(summary_path, keep_default_na=False, low_memory=False),
            allow_unversioned_v1=allow_unversioned_v1,
        )
        if df.empty:
            continue
        df = df.copy()
        numeric_cols = [
            "Pos",
            "NumReads",
            "NumCpGs",
            "HP1FamilyN",
            "HP2FamilyN",
            "NHP0",
            "NHP3",
            "AlleleDelta",
            "PairwiseMedianDist",
            "Quality_Score",
            "CramersV",
            "HP_Ratio",
        ]
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        df["truth_label"] = truth_label
        df["variant_key"] = (
            df["Chr"].astype(str)
            + ":"
            + df["Pos"].astype("Int64").astype(str)
            + ":"
            + df["Ref"].astype(str)
            + ":"
            + df["Alt"].astype(str)
        )
        df["sample"] = config.sample
        df["sample_label"] = config.sample_label
        df["platform"] = config.platform
        df["mode"] = config.mode
        df["mode_label"] = config.mode_label
        df["dataset_id"] = config.dataset_id
        df["dataset_label"] = config.dataset_label
        df["source_kind"] = config.source_kind
        df["source_priority"] = config.source_priority
        df["source_base_dir"] = str(config.base_dir)
        df["source_summary_file"] = str(summary_path)
        df["source_region_root"] = str(config.tp_region_root if truth_label == "TP" else config.fp_region_root)
        df["source_vcf_file"] = str(config.tp_vcf if truth_label == "TP" else config.fp_vcf)
        df["source_tagged_bam"] = str(config.tagged_bam)
        df["source_phased_vcf"] = str(config.phased_vcf) if config.phased_vcf else ""
        df["source_loh_bed"] = str(config.loh_bed) if config.loh_bed else ""
        df["context_truth_total"] = context.get("truth_total", "")
        df["effective_hp_reads"] = df["HP1FamilyN"].fillna(0) + df["HP2FamilyN"].fillna(0)
        df["hp_ratio_core"] = np.where(
            df["effective_hp_reads"] > 0,
            df["HP1FamilyN"].fillna(0) / df["effective_hp_reads"],
            np.nan,
        )
        df["hp_assign_rate"] = np.where(
            df["NumReads"] > 0,
            df["effective_hp_reads"] / df["NumReads"],
            np.nan,
        )
        df["hp0_ratio"] = np.where(df["NumReads"] > 0, df["NHP0"].fillna(0) / df["NumReads"], np.nan)
        df["hp3_ratio"] = np.where(df["NumReads"] > 0, df["NHP3"].fillna(0) / df["NumReads"], np.nan)
        df["tool_potential_loh"] = df["Potential_LOH"].map(bool_from_value)
        df["core_loh_like"] = (df["effective_hp_reads"] > 0) & (
            (df["hp_ratio_core"] < 0.1) | (df["hp_ratio_core"] > 0.9)
        )
        df["tool_core_loh_match"] = np.where(
            df["effective_hp_reads"] > 0,
            df["tool_potential_loh"] == df["core_loh_like"],
            False,
        )
        df["hp_ratio_bin"] = [
            hp_ratio_bin(ratio, eff)
            for ratio, eff in zip(df["hp_ratio_core"].tolist(), df["effective_hp_reads"].tolist())
        ]
        df["allele_delta"] = df["AlleleDelta"]
        df["pairwise_median_dist"] = df["PairwiseMedianDist"]
        df["quality_score"] = df["Quality_Score"]
        df["cramers_v"] = df["CramersV"]
        df["num_reads"] = df["NumReads"]
        df["num_cpgs"] = df["NumCpGs"]
        vcf_df = vcf_maps[truth_label]
        if not vcf_df.empty:
            df = df.merge(vcf_df, on="variant_key", how="left")
        else:
            for col in [
                "caller_filter",
                "caller_af",
                "caller_dp",
                "caller_gq",
                "caller_ad_ref",
                "caller_ad_alt",
                "caller_naf",
                "caller_ndp",
                "caller_nad_ref",
                "caller_nad_alt",
            ]:
                df[col] = np.nan if col != "caller_filter" else ""
        df["af_bin"] = [af_bin(x) for x in df["caller_af"].tolist()]
        df["effective_hp_bin"] = [effective_hp_bin(x) for x in df["effective_hp_reads"].tolist()]
        if config.mode == "to":
            if not to_phase_df.empty:
                df = df.merge(to_phase_df, on="variant_key", how="left")
            else:
                for col in [
                    "to_phase_filter",
                    "to_phase_ps",
                    "to_phase_gt",
                    "to_phase_gt2",
                    "to_phase_ps2",
                    "to_phase_gt3",
                    "to_phase_ps3",
                ]:
                    df[col] = ""
            df["to_loh_bed_hit"] = [
                bed_index.contains(chrom, int(pos))
                for chrom, pos in zip(df["Chr"].astype(str).tolist(), df["Pos"].fillna(0).astype(int).tolist())
            ]
            df["phase_block_status"] = np.where(
                df["to_phase_ps"].fillna("").astype(str) != "",
                "to_ps_available",
                "to_ps_missing",
            )
        else:
            for col in [
                "to_phase_filter",
                "to_phase_ps",
                "to_phase_gt",
                "to_phase_gt2",
                "to_phase_ps2",
                "to_phase_gt3",
                "to_phase_ps3",
            ]:
                df[col] = ""
            df["to_loh_bed_hit"] = False
            df["phase_block_status"] = "paired_vcf_ps_missing"
        frames.append(df.copy())
    if not frames:
        return pd.DataFrame()
    result = pd.concat(frames, ignore_index=True)
    statuses = sorted(result["VerificationProvenanceStatus"].dropna().unique().tolist())
    if len(statuses) != 1:
        raise SchemaContractError(
            f"LOH round provenance: mixed TP/FP schema statuses for {config.dataset_id}: {statuses}"
        )
    result = result.sort_values(["sample_label", "mode", "truth_label", "Chr", "Pos"]).reset_index(drop=True)
    return result


def load_benchmark_summary(configs: Sequence[DatasetConfig]) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    for config in configs:
        if not config.benchmark_tsv.exists():
            continue
        df = pd.read_csv(config.benchmark_tsv, sep="\t")
        df["dataset_id"] = config.dataset_id
        df["dataset_label"] = config.dataset_label
        df["sample_label"] = config.sample_label
        df["mode_group"] = config.mode
        df["mode_label"] = config.mode_label
        df["source_kind"] = config.source_kind
        df["benchmark_tsv"] = str(config.benchmark_tsv)
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def attach_summary_metadata(df: pd.DataFrame, source_file: Path, analysis_rule: str, tool_role: str, generated_at: str) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    df["tool_role"] = tool_role
    df["source_file"] = str(source_file)
    if "source_record_count" not in df.columns:
        df["source_record_count"] = df.get("count", 0)
    df["analysis_subset_rule"] = analysis_rule
    df["generated_at"] = generated_at
    return df


def summarize_loh_distribution(all_df: pd.DataFrame, source_file: Path, generated_at: str) -> pd.DataFrame:
    totals = (
        all_df.groupby(["sample", "sample_label", "mode", "mode_label", "truth_label"], dropna=False)
        .size()
        .rename("truth_group_total")
        .reset_index()
    )
    grouped = (
        all_df.groupby(["sample", "sample_label", "mode", "mode_label", "truth_label", "hp_ratio_bin"], dropna=False)
        .agg(
            count=("variant_key", "size"),
            core_loh_like_count=("core_loh_like", "sum"),
            tool_potential_loh_count=("tool_potential_loh", "sum"),
            tool_core_mismatch_count=("tool_core_loh_match", lambda s: int((~s).sum())),
            median_effective_hp_reads=("effective_hp_reads", "median"),
            median_num_reads=("num_reads", "median"),
            median_hp_assign_rate=("hp_assign_rate", "median"),
        )
        .reset_index()
    )
    grouped = grouped.merge(totals, on=["sample", "sample_label", "mode", "mode_label", "truth_label"], how="left")
    grouped["fraction_of_truth_group"] = grouped["count"] / grouped["truth_group_total"]
    grouped["source_record_count"] = grouped["count"]
    grouped = attach_summary_metadata(
        grouped,
        source_file,
        "group all region rows by hp_ratio_bin within sample/mode/truth",
        "loh_distribution_summary",
        generated_at,
    )
    return grouped


def summarize_loh_enrichment(all_df: pd.DataFrame, source_file: Path, generated_at: str) -> pd.DataFrame:
    rows: List[Dict] = []
    for (sample, sample_label, mode, mode_label), group in all_df.groupby(["sample", "sample_label", "mode", "mode_label"], dropna=False):
        tp_group = group[group["truth_label"] == "TP"]
        fp_group = group[group["truth_label"] == "FP"]
        tp_total = len(tp_group)
        fp_total = len(fp_group)
        tp_loh = int(tp_group["core_loh_like"].sum())
        fp_loh = int(fp_group["core_loh_like"].sum())
        tp_tool = int(tp_group["tool_potential_loh"].sum())
        fp_tool = int(fp_group["tool_potential_loh"].sum())
        tp_frac = tp_loh / tp_total if tp_total else math.nan
        fp_frac = fp_loh / fp_total if fp_total else math.nan
        rows.append(
            {
                "sample": sample,
                "sample_label": sample_label,
                "mode": mode,
                "mode_label": mode_label,
                "tp_total": tp_total,
                "fp_total": fp_total,
                "tp_core_loh_like_n": tp_loh,
                "fp_core_loh_like_n": fp_loh,
                "tp_core_loh_like_frac": tp_frac,
                "fp_core_loh_like_frac": fp_frac,
                "core_loh_fp_vs_tp_enrichment": fp_frac / tp_frac if tp_total and tp_loh else math.nan,
                "tp_tool_potential_loh_n": tp_tool,
                "fp_tool_potential_loh_n": fp_tool,
                "tp_tool_potential_loh_frac": tp_tool / tp_total if tp_total else math.nan,
                "fp_tool_potential_loh_frac": fp_tool / fp_total if fp_total else math.nan,
                "median_tp_effective_hp_reads": tp_group["effective_hp_reads"].median() if tp_total else math.nan,
                "median_fp_effective_hp_reads": fp_group["effective_hp_reads"].median() if fp_total else math.nan,
                "source_record_count": len(group),
            }
        )
    df = pd.DataFrame.from_records(rows)
    return attach_summary_metadata(
        df,
        source_file,
        "compare TP vs FP LOH-like proportions within sample/mode",
        "loh_enrichment_summary",
        generated_at,
    )


def summarize_hp_balance(all_df: pd.DataFrame, source_file: Path, generated_at: str) -> pd.DataFrame:
    grouped = (
        all_df.groupby(["sample", "sample_label", "mode", "mode_label", "truth_label"], dropna=False)
        .agg(
            count=("variant_key", "size"),
            median_hp_ratio_core=("hp_ratio_core", "median"),
            q10_hp_ratio_core=("hp_ratio_core", lambda s: s.quantile(0.10)),
            q90_hp_ratio_core=("hp_ratio_core", lambda s: s.quantile(0.90)),
            balanced_fraction=("hp_ratio_core", lambda s: float(((s >= 0.4) & (s <= 0.6)).mean())),
            extreme_fraction=("core_loh_like", "mean"),
            no_effective_hp_fraction=("effective_hp_reads", lambda s: float((s <= 0).mean())),
            median_effective_hp_reads=("effective_hp_reads", "median"),
            median_num_reads=("num_reads", "median"),
        )
        .reset_index()
    )
    grouped["source_record_count"] = grouped["count"]
    return attach_summary_metadata(
        grouped,
        source_file,
        "summarize hp balance and effective HP support within sample/mode/truth",
        "hp_family_balance_summary",
        generated_at,
    )


def summarize_hp_qc(all_df: pd.DataFrame, source_file: Path, generated_at: str) -> pd.DataFrame:
    grouped = (
        all_df.groupby(["sample", "sample_label", "mode", "mode_label", "truth_label"], dropna=False)
        .agg(
            count=("variant_key", "size"),
            median_hp_assign_rate=("hp_assign_rate", "median"),
            median_hp0_ratio=("hp0_ratio", "median"),
            median_hp3_ratio=("hp3_ratio", "median"),
            hp0_ge_0_5_fraction=("hp0_ratio", lambda s: float((s >= 0.5).mean())),
            hp3_gt_0_fraction=("hp3_ratio", lambda s: float((s > 0).mean())),
            effective_hp_lt10_fraction=("effective_hp_reads", lambda s: float((s < 10).mean())),
            effective_hp_eq0_fraction=("effective_hp_reads", lambda s: float((s == 0).mean())),
            median_quality_score=("quality_score", "median"),
        )
        .reset_index()
    )
    grouped["source_record_count"] = grouped["count"]
    return attach_summary_metadata(
        grouped,
        source_file,
        "summarize phase/QC related ratios within sample/mode/truth",
        "hp_coverage_qc_summary",
        generated_at,
    )


def summarize_loh_vs_features(all_df: pd.DataFrame, source_file: Path, generated_at: str) -> pd.DataFrame:
    grouped = (
        all_df.groupby(["sample", "sample_label", "mode", "mode_label", "truth_label", "core_loh_like"], dropna=False)
        .agg(
            count=("variant_key", "size"),
            median_allele_delta=("allele_delta", "median"),
            median_pairwise_median_dist=("pairwise_median_dist", "median"),
            median_quality_score=("quality_score", "median"),
            median_hp_assign_rate=("hp_assign_rate", "median"),
            median_hp0_ratio=("hp0_ratio", "median"),
            median_hp3_ratio=("hp3_ratio", "median"),
            median_caller_af=("caller_af", "median"),
        )
        .reset_index()
    )
    grouped["core_loh_label"] = np.where(grouped["core_loh_like"], "LOH-like", "non-LOH-like")
    grouped["source_record_count"] = grouped["count"]
    return attach_summary_metadata(
        grouped,
        source_file,
        "compare methylation and QC features between LOH-like and non-LOH-like groups",
        "loh_vs_methyl_feature_summary",
        generated_at,
    )


def summarize_verification_by_loh(all_df: pd.DataFrame, source_file: Path, generated_at: str) -> pd.DataFrame:
    provenance_group_fields = [
        *VERIFICATION_PROVENANCE_COLUMNS,
        "VerificationProvenanceStatus",
        "VerificationProvenanceSourceField",
        "VerificationProvenanceWarnings",
        "LOHProvenanceSourceField",
        "VerificationUnknownCurrentCounts",
    ]
    grouped = (
        all_df.groupby(
            [
                "sample",
                "sample_label",
                "mode",
                "mode_label",
                "truth_label",
                "verification_class",
                "loh_subtype",
                *provenance_group_fields,
            ],
            dropna=False,
        )
        .size()
        .rename("count")
        .reset_index()
    )
    totals = (
        all_df.groupby(
            [
                "sample",
                "sample_label",
                "mode",
                "mode_label",
                "truth_label",
                "verification_class",
                *provenance_group_fields,
            ],
            dropna=False,
        )
        .size()
        .rename("class_total")
        .reset_index()
    )
    grouped = grouped.merge(
        totals,
        on=[
            "sample",
            "sample_label",
            "mode",
            "mode_label",
            "truth_label",
            "verification_class",
            *provenance_group_fields,
        ],
        how="left",
    )
    grouped["fraction_within_verification_class"] = grouped["count"] / grouped["class_total"]
    grouped["source_record_count"] = grouped["count"]
    return attach_summary_metadata(
        grouped,
        source_file,
        "cross-tab verification_class and loh_subtype within sample/mode/truth",
        "verificationclass_by_loh_subtype",
        generated_at,
    )


def build_heatmap_table(all_df: pd.DataFrame, source_file: Path, generated_at: str) -> pd.DataFrame:
    rows: List[Dict] = []
    for metric, bin_col in [("hp_ratio_core", "hp_ratio_bin"), ("caller_af", "af_bin"), ("effective_hp_reads", "effective_hp_bin")]:
        grouped = (
            all_df.groupby(["sample", "sample_label", "mode", "mode_label", bin_col], dropna=False)
            .size()
            .rename("count")
            .reset_index()
        )
        totals = (
            all_df.groupby(["sample", "sample_label", "mode", "mode_label"], dropna=False)
            .size()
            .rename("total")
            .reset_index()
        )
        grouped = grouped.merge(totals, on=["sample", "sample_label", "mode", "mode_label"], how="left")
        for row in grouped.to_dict("records"):
            rows.append(
                {
                    "sample": row["sample"],
                    "sample_label": row["sample_label"],
                    "mode": row["mode"],
                    "mode_label": row["mode_label"],
                    "metric": metric,
                    "bin_label": row[bin_col],
                    "count": row["count"],
                    "fraction": row["count"] / row["total"] if row["total"] else math.nan,
                }
            )
    df = pd.DataFrame.from_records(rows)
    df["source_record_count"] = df["count"]
    return attach_summary_metadata(
        df,
        source_file,
        "sample/mode bin fractions for hp_ratio_core, caller_af, effective_hp_reads",
        "sample_bin_heatmap",
        generated_at,
    )


def compare_status(paired_status: str, to_status: str, overlap_class: str) -> str:
    if overlap_class == "both":
        if paired_status == "TP" and to_status == "TP":
            return "both_tp"
        if paired_status == "FP" and to_status == "FP":
            return "both_fp"
        if paired_status == "TP" and to_status == "FP":
            return "paired_tp_to_fp"
        if paired_status == "FP" and to_status == "TP":
            return "paired_fp_to_tp"
        return "both_other"
    if overlap_class == "paired_only":
        return f"paired_only_{paired_status.lower()}" if paired_status else "paired_only"
    if overlap_class == "to_only":
        return f"to_only_{to_status.lower()}" if to_status else "to_only"
    return "unknown"


def build_same_locus_compare(all_df: pd.DataFrame) -> pd.DataFrame:
    keep_cols = [
        "sample",
        "sample_label",
        "platform",
        "variant_key",
        "Chr",
        "Pos",
        "Ref",
        "Alt",
        "truth_label",
        "verification_class",
        "loh_subtype",
        "core_loh_like",
        "tool_potential_loh",
        "hp_ratio_core",
        "effective_hp_reads",
        "hp_assign_rate",
        "hp0_ratio",
        "hp3_ratio",
        "allele_delta",
        "pairwise_median_dist",
        "quality_score",
        "caller_af",
        "caller_dp",
        "caller_gq",
        "caller_filter",
        "to_phase_filter",
        "to_phase_ps",
        "to_phase_gt",
        "to_phase_gt2",
        "to_phase_ps2",
        "to_phase_gt3",
        "to_phase_ps3",
        "to_loh_bed_hit",
        "source_summary_file",
        "source_vcf_file",
        "source_tagged_bam",
        "phase_block_status",
        *VERIFICATION_PROVENANCE_COLUMNS,
        *PROVENANCE_METADATA_FIELDS,
    ]
    rows: List[pd.DataFrame] = []
    for sample in sorted(all_df["sample"].dropna().unique().tolist()):
        paired = all_df[(all_df["sample"] == sample) & (all_df["mode"] == "paired")][keep_cols].copy()
        to = all_df[(all_df["sample"] == sample) & (all_df["mode"] == "to")][keep_cols].copy()
        if paired.empty and to.empty:
            continue
        merged = paired.merge(
            to,
            on=["sample", "variant_key"],
            how="outer",
            suffixes=("_paired", "_to"),
        )
        merged["sample_label"] = merged["sample_label_paired"].combine_first(merged["sample_label_to"])
        merged["platform"] = merged["platform_paired"].combine_first(merged["platform_to"])
        merged["chrom"] = merged["Chr_paired"].combine_first(merged["Chr_to"])
        merged["pos"] = merged["Pos_paired"].combine_first(merged["Pos_to"])
        merged["ref"] = merged["Ref_paired"].combine_first(merged["Ref_to"])
        merged["alt"] = merged["Alt_paired"].combine_first(merged["Alt_to"])
        merged["paired_status"] = merged["truth_label_paired"].fillna("")
        merged["to_status"] = merged["truth_label_to"].fillna("")
        merged["paired_present"] = merged["paired_status"] != ""
        merged["to_present"] = merged["to_status"] != ""
        merged["locus_overlap_class"] = np.select(
            [merged["paired_present"] & merged["to_present"], merged["paired_present"], merged["to_present"]],
            ["both", "paired_only", "to_only"],
            default="none",
        )
        merged["paired_to_concordance"] = [
            compare_status(paired_status, to_status, overlap)
            for paired_status, to_status, overlap in zip(
                merged["paired_status"].tolist(),
                merged["to_status"].tolist(),
                merged["locus_overlap_class"].tolist(),
            )
        ]
        merged["paired_verification_class"] = merged["verification_class_paired"].fillna("")
        merged["to_verification_class"] = merged["verification_class_to"].fillna("")
        merged["paired_loh_subtype"] = merged["loh_subtype_paired"].fillna("")
        merged["to_loh_subtype"] = merged["loh_subtype_to"].fillna("")
        merged["max_effective_hp_reads"] = merged[["effective_hp_reads_paired", "effective_hp_reads_to"]].max(axis=1, skipna=True)
        rows.append(merged)
    if not rows:
        return pd.DataFrame()
    df = pd.concat(rows, ignore_index=True)
    order_cols = [
        "sample",
        "sample_label",
        "platform",
        "variant_key",
        "chrom",
        "pos",
        "ref",
        "alt",
        "locus_overlap_class",
        "paired_to_concordance",
        "paired_status",
        "to_status",
        "paired_verification_class",
        "to_verification_class",
        "paired_loh_subtype",
        "to_loh_subtype",
        "core_loh_like_paired",
        "core_loh_like_to",
        "tool_potential_loh_paired",
        "tool_potential_loh_to",
        "hp_ratio_core_paired",
        "hp_ratio_core_to",
        "effective_hp_reads_paired",
        "effective_hp_reads_to",
        "hp_assign_rate_paired",
        "hp_assign_rate_to",
        "hp0_ratio_paired",
        "hp0_ratio_to",
        "hp3_ratio_paired",
        "hp3_ratio_to",
        "allele_delta_paired",
        "allele_delta_to",
        "pairwise_median_dist_paired",
        "pairwise_median_dist_to",
        "quality_score_paired",
        "quality_score_to",
        "caller_af_paired",
        "caller_af_to",
        "caller_filter_paired",
        "caller_filter_to",
        "to_phase_ps_to",
        "to_phase_gt_to",
        "to_phase_gt2_to",
        "to_phase_ps2_to",
        "to_phase_gt3_to",
        "to_phase_ps3_to",
        "to_loh_bed_hit_to",
        "source_summary_file_paired",
        "source_summary_file_to",
        "source_vcf_file_paired",
        "source_vcf_file_to",
        "source_tagged_bam_paired",
        "source_tagged_bam_to",
    ]
    existing_cols = [col for col in order_cols if col in df.columns]
    other_cols = [col for col in df.columns if col not in existing_cols]
    return df[existing_cols + other_cols].sort_values(["sample_label", "pos", "variant_key"]).reset_index(drop=True)


def locate_region_dir(region_root: Path, chrom: str, pos: int) -> Optional[Path]:
    base = region_root / chrom / f"{chrom}_{pos}"
    if not base.exists():
        return None
    children = [child for child in base.iterdir() if child.is_dir()]
    if children:
        return sorted(children)[0]
    return base


def pileup_hp_counts(bam_path: Path, chrom: str, pos: int, ref: str, alt: str) -> tuple[List[Dict], List[Dict]]:
    if not bam_path.exists():
        return [], []
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    per_read: Dict[str, Dict] = {}
    for column in bam.pileup(chrom, pos - 1, pos, truncate=True, stepper="samtools", max_depth=100000):
        if column.reference_pos != pos - 1:
            continue
        for pileupread in column.pileups:
            aln = pileupread.alignment
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            hp_tag = normalize_hp_tag(aln.get_tag("HP") if aln.has_tag("HP") else None)
            ps_value = aln.get_tag("PS") if aln.has_tag("PS") else ""
            if pileupread.is_del or pileupread.is_refskip:
                base = "-"
                base_category = "deletion_or_skip"
            else:
                base = aln.query_sequence[pileupread.query_position].upper()
                if base == ref.upper():
                    base_category = "ref"
                elif base == alt.upper():
                    base_category = "alt"
                else:
                    base_category = "other"
            per_read[aln.query_name] = {
                "read_name": aln.query_name,
                "hp_tag": hp_tag,
                "hp_family": hp_family(hp_tag),
                "ps": ps_value,
                "base": base,
                "base_category": base_category,
                "mapq": aln.mapping_quality,
                "read_start": aln.reference_start + 1,
                "read_end": aln.reference_end,
            }
    bam.close()
    detailed_rows = sorted(per_read.values(), key=lambda row: (str(row["ps"]), row["hp_tag"], row["read_name"]))
    counter = Counter((row["hp_tag"], row["hp_family"], str(row["ps"]), row["base_category"]) for row in detailed_rows)
    summary_rows = [
        {
            "hp_tag": hp_tag,
            "hp_family": hp_family_name,
            "ps": ps_value,
            "base_category": base_category,
            "read_count": count,
        }
        for (hp_tag, hp_family_name, ps_value, base_category), count in sorted(counter.items())
    ]
    return detailed_rows, summary_rows


def build_case_feature_plot(case_dir: Path, rows: List[Dict]) -> str:
    data_rows = [row for row in rows if row.get("present")]
    if not data_rows:
        return ""
    feature_names = [
        ("hp_ratio_core", "hp_ratio_core"),
        ("effective_hp_reads", "effective_hp_reads"),
        ("hp_assign_rate", "hp_assign_rate"),
        ("hp0_ratio", "hp0_ratio"),
        ("hp3_ratio", "hp3_ratio"),
        ("allele_delta", "allele_delta"),
        ("pairwise_median_dist", "pairwise_median_dist"),
        ("quality_score", "quality_score"),
        ("caller_af", "caller_af"),
    ]
    raw_values = np.array(
        [[to_float_or_nan(row.get(key)) for _, key in feature_names] for row in data_rows],
        dtype=float,
    )
    normalized = raw_values.copy()
    for col_idx in range(normalized.shape[1]):
        col = normalized[:, col_idx]
        finite = col[np.isfinite(col)]
        if finite.size <= 1:
            normalized[:, col_idx] = 0.5 if finite.size == 1 else np.nan
            continue
        cmin = finite.min()
        cmax = finite.max()
        if cmax == cmin:
            normalized[:, col_idx] = 0.5
        else:
            normalized[:, col_idx] = (col - cmin) / (cmax - cmin)
    fig, ax = plt.subplots(figsize=(10, max(2.5, 0.9 * len(data_rows))))
    im = ax.imshow(normalized, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(feature_names)))
    ax.set_xticklabels([name for name, _ in feature_names], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(data_rows)))
    ax.set_yticklabels([row["mode"] for row in data_rows])
    for i in range(len(data_rows)):
        for j in range(len(feature_names)):
            value = raw_values[i, j]
            text = "NA" if math.isnan(value) else f"{value:.3f}"
            ax.text(j, i, text, ha="center", va="center", color="white", fontsize=8)
    ax.set_title("Case feature comparison")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="pairwise normalized value")
    fig.tight_layout()
    out_path = case_dir / "feature_comparison.png"
    fig.savefig(out_path, dpi=220)
    plt.close(fig)
    return str(out_path)


def build_case_hp_count_plot(case_dir: Path, case_id: str, summary_rows: Dict[str, List[Dict]]) -> str:
    mode_names = [mode for mode, rows in summary_rows.items() if rows]
    if not mode_names:
        return ""
    hp_family_order = ["HP1_family", "HP2_family", "HP0", "HP3", "OTHER"]
    fig, axes = plt.subplots(1, len(mode_names), figsize=(5 * len(mode_names), 4), squeeze=False)
    for idx, mode in enumerate(mode_names):
        ax = axes[0][idx]
        counts = Counter()
        for row in summary_rows[mode]:
            counts[row["hp_family"]] += int(row["read_count"])
        values = [counts.get(name, 0) for name in hp_family_order]
        ax.bar(hp_family_order, values, color=["#1f77b4", "#ff7f0e", "#7f7f7f", "#2ca02c", "#d62728"])
        ax.set_title(f"{case_id} {mode}")
        ax.set_ylabel("read count")
        ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    out_path = case_dir / "hp_tag_counts.png"
    fig.savefig(out_path, dpi=220)
    plt.close(fig)
    return str(out_path)


def describe_phase_status(row: Dict) -> str:
    status_bits: List[str] = []
    if row.get("to_phase_ps"):
        status_bits.append("TO phased VCF PS available")
    else:
        status_bits.append("TO phased VCF PS missing")
    status_bits.append("paired VCF PS unavailable; paired PS only inspectable from tagged BAM read tags")
    status_bits.append("different PS blocks cannot be treated as parent/maternal labels across blocks")
    return "; ".join(status_bits)


def build_case_panel(
    all_df: pd.DataFrame,
    same_df: pd.DataFrame,
    output_dir: Path,
) -> pd.DataFrame:
    cases_dir = ensure_dir(output_dir / "cases")
    selected_rows: List[Dict] = []
    used_keys: set[tuple[str, str]] = set()

    def add_case(case_category: str, candidate_df: pd.DataFrame) -> None:
        if candidate_df.empty:
            return
        sorted_df = candidate_df.sort_values(
            ["max_effective_hp_reads", "pos"],
            ascending=[False, True],
        )
        for row in sorted_df.to_dict("records"):
            key = (row["sample"], row["variant_key"])
            if key in used_keys:
                continue
            row["case_category"] = case_category
            selected_rows.append(row)
            used_keys.add(key)
            return

    add_case(
        "both_tp_loh_like",
        same_df[
            (same_df["paired_status"] == "TP")
            & (same_df["to_status"] == "TP")
            & ((same_df["core_loh_like_paired"]) | (same_df["core_loh_like_to"]))
        ],
    )
    add_case(
        "both_fp_loh_like",
        same_df[
            (same_df["paired_status"] == "FP")
            & (same_df["to_status"] == "FP")
            & ((same_df["core_loh_like_paired"]) | (same_df["core_loh_like_to"]))
        ],
    )
    add_case(
        "paired_tp_to_fp",
        same_df[(same_df["paired_status"] == "TP") & (same_df["to_status"] == "FP")],
    )
    add_case(
        "paired_fp_to_tp",
        same_df[(same_df["paired_status"] == "FP") & (same_df["to_status"] == "TP")],
    )
    add_case(
        "paired_only_loh_like",
        same_df[(same_df["locus_overlap_class"] == "paired_only") & (same_df["core_loh_like_paired"])],
    )
    add_case(
        "to_only_loh_like",
        same_df[(same_df["locus_overlap_class"] == "to_only") & (same_df["core_loh_like_to"])],
    )
    if len(selected_rows) < 8:
        high_noise = all_df[
            ((all_df["hp0_ratio"] >= 0.5) | (all_df["hp3_ratio"] > 0))
            & (all_df["num_reads"] >= 10)
        ].sort_values(["hp0_ratio", "hp3_ratio", "effective_hp_reads"], ascending=[False, False, False])
        for raw in high_noise.to_dict("records"):
            key = (raw["sample"], raw["variant_key"])
            if key in used_keys:
                continue
            selected_rows.append(
                {
                    "sample": raw["sample"],
                    "sample_label": raw["sample_label"],
                    "variant_key": raw["variant_key"],
                    "chrom": raw["Chr"],
                    "pos": int(raw["Pos"]),
                    "ref": raw["Ref"],
                    "alt": raw["Alt"],
                    "locus_overlap_class": raw["mode"],
                    "paired_status": raw["truth_label"] if raw["mode"] == "paired" else "",
                    "to_status": raw["truth_label"] if raw["mode"] == "to" else "",
                    "paired_verification_class": raw["verification_class"] if raw["mode"] == "paired" else "",
                    "to_verification_class": raw["verification_class"] if raw["mode"] == "to" else "",
                    "core_loh_like_paired": bool(raw["core_loh_like"]) if raw["mode"] == "paired" else False,
                    "core_loh_like_to": bool(raw["core_loh_like"]) if raw["mode"] == "to" else False,
                    "hp_ratio_core_paired": raw["hp_ratio_core"] if raw["mode"] == "paired" else math.nan,
                    "hp_ratio_core_to": raw["hp_ratio_core"] if raw["mode"] == "to" else math.nan,
                    "effective_hp_reads_paired": raw["effective_hp_reads"] if raw["mode"] == "paired" else math.nan,
                    "effective_hp_reads_to": raw["effective_hp_reads"] if raw["mode"] == "to" else math.nan,
                    "max_effective_hp_reads": raw["effective_hp_reads"],
                    "case_category": "high_hp0_or_hp3",
                }
            )
            used_keys.add(key)
            break

    registry_rows: List[Dict] = []
    for idx, base_row in enumerate(selected_rows, start=1):
        sample = base_row["sample"]
        variant = base_row["variant_key"]
        chrom = str(base_row["chrom"])
        pos = int(base_row["pos"])
        ref = str(base_row["ref"])
        alt = str(base_row["alt"])
        case_id = f"case_{idx:02d}_{base_row['case_category']}_{sample}_{chrom}_{pos}"
        case_dir = ensure_dir(cases_dir / case_id)
        mode_rows: List[Dict] = []
        same_rows = same_df[(same_df["sample"] == sample) & (same_df["variant_key"] == variant)]
        if not same_rows.empty:
            same_row = same_rows.iloc[0].to_dict()
        else:
            same_row = {}
        for mode in ["paired", "to"]:
            subset = all_df[(all_df["sample"] == sample) & (all_df["mode"] == mode) & (all_df["variant_key"] == variant)]
            if subset.empty:
                mode_rows.append({"mode": mode, "present": False})
                continue
            row = subset.iloc[0].to_dict()
            mode_rows.append(
                {
                    "mode": mode,
                    "present": True,
                    "truth_label": row["truth_label"],
                    "verification_class": row["verification_class"],
                    "loh_subtype": row["loh_subtype"],
                    "core_loh_like": bool(row["core_loh_like"]),
                    "tool_potential_loh": bool(row["tool_potential_loh"]),
                    "hp_ratio_core": row["hp_ratio_core"],
                    "effective_hp_reads": row["effective_hp_reads"],
                    "hp_assign_rate": row["hp_assign_rate"],
                    "hp0_ratio": row["hp0_ratio"],
                    "hp3_ratio": row["hp3_ratio"],
                    "allele_delta": row["allele_delta"],
                    "pairwise_median_dist": row["pairwise_median_dist"],
                    "quality_score": row["quality_score"],
                    "caller_af": row["caller_af"],
                    "source_summary_file": row["source_summary_file"],
                    "source_vcf_file": row["source_vcf_file"],
                    "source_tagged_bam": row["source_tagged_bam"],
                    "source_region_root": row["source_region_root"],
                    "source_phased_vcf": row["source_phased_vcf"],
                    "source_loh_bed": row["source_loh_bed"],
                    "to_phase_ps": row["to_phase_ps"],
                    "to_phase_gt": row["to_phase_gt"],
                    "to_phase_gt2": row["to_phase_gt2"],
                    "to_phase_ps2": row["to_phase_ps2"],
                    "to_phase_gt3": row["to_phase_gt3"],
                    "to_phase_ps3": row["to_phase_ps3"],
                    "to_loh_bed_hit": bool(row["to_loh_bed_hit"]),
                    "phase_block_status": row["phase_block_status"],
                    "region_dir": str(locate_region_dir(Path(row["source_region_root"]), chrom, pos) or ""),
                }
            )
        mode_rows_by_name = {row["mode"]: row for row in mode_rows}
        evidence_rows: List[Dict] = []
        bam_summary_by_mode: Dict[str, List[Dict]] = {}
        for mode in ["paired", "to"]:
            row = mode_rows_by_name[mode]
            if not row.get("present"):
                continue
            bam_path = Path(str(row["source_tagged_bam"]))
            detailed_rows, summary_rows = pileup_hp_counts(bam_path, chrom, pos, ref, alt)
            bam_summary_by_mode[mode] = summary_rows
            if summary_rows:
                write_tsv_rows(case_dir / f"{mode}_tag_bam_counts.tsv", ["hp_tag", "hp_family", "ps", "base_category", "read_count"], summary_rows)
            else:
                write_tsv_rows(case_dir / f"{mode}_tag_bam_counts.tsv", ["hp_tag", "hp_family", "ps", "base_category", "read_count"], [])
            if detailed_rows:
                pd.DataFrame.from_records(detailed_rows).to_csv(
                    case_dir / f"{mode}_tag_bam_reads.tsv.gz",
                    sep="\t",
                    index=False,
                    compression="gzip",
                )
            region_dir = row.get("region_dir", "")
            if region_dir:
                evidence_rows.extend(
                    [
                        {"artifact_type": f"{mode}_region_dir", "path": region_dir, "exists": Path(region_dir).exists(), "notes": "region root"},
                        {"artifact_type": f"{mode}_clustering_dir", "path": str(Path(region_dir) / "clustering"), "exists": (Path(region_dir) / "clustering").exists(), "notes": "existing clustering outputs"},
                        {"artifact_type": f"{mode}_distance_dir", "path": str(Path(region_dir) / "distance"), "exists": (Path(region_dir) / "distance").exists(), "notes": "existing distance outputs"},
                        {"artifact_type": f"{mode}_methylation_dir", "path": str(Path(region_dir) / "methylation"), "exists": (Path(region_dir) / "methylation").exists(), "notes": "existing methylation outputs"},
                    ]
                )
            evidence_rows.extend(
                [
                    {"artifact_type": f"{mode}_summary", "path": str(row["source_summary_file"]), "exists": Path(str(row["source_summary_file"])).exists(), "notes": "significance summary source"},
                    {"artifact_type": f"{mode}_vcf", "path": str(row["source_vcf_file"]), "exists": Path(str(row["source_vcf_file"])).exists(), "notes": "benchmark split VCF source"},
                    {"artifact_type": f"{mode}_tagged_bam", "path": str(row["source_tagged_bam"]), "exists": bam_path.exists(), "notes": "tagged BAM used for read-level counts"},
                ]
            )
            if row.get("source_phased_vcf"):
                evidence_rows.append(
                    {
                        "artifact_type": f"{mode}_phased_vcf",
                        "path": str(row["source_phased_vcf"]),
                        "exists": Path(str(row["source_phased_vcf"])).exists(),
                        "notes": "TO phased VCF for PS lookup",
                    }
                )
            if row.get("source_loh_bed"):
                evidence_rows.append(
                    {
                        "artifact_type": f"{mode}_loh_bed",
                        "path": str(row["source_loh_bed"]),
                        "exists": Path(str(row["source_loh_bed"])).exists(),
                        "notes": "TO LOH bed reference",
                    }
                )
        feature_plot = build_case_feature_plot(case_dir, mode_rows)
        hp_plot = build_case_hp_count_plot(case_dir, case_id, bam_summary_by_mode)
        write_tsv_rows(case_dir / "evidence_links.tsv", ["artifact_type", "path", "exists", "notes"], evidence_rows)
        context_payload = {
            "case_id": case_id,
            "case_category": base_row["case_category"],
            "sample": sample,
            "sample_label": base_row["sample_label"],
            "variant_key": variant,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "locus_overlap_class": base_row.get("locus_overlap_class", ""),
            "paired_to_concordance": same_row.get("paired_to_concordance", ""),
            "phase_block_status": describe_phase_status(mode_rows_by_name.get("to", {})),
            "figure_paths": {
                "feature_comparison": feature_plot,
                "hp_tag_counts": hp_plot,
            },
            "source_table_rows": {
                "same_locus_compare_variant_key": variant,
                "all_region_sample": sample,
            },
        }
        write_json(case_dir / "case_context.json", context_payload)
        report_lines = [
            f"# {case_id}",
            "",
            "## Basic Information",
            "",
            f"- Sample: `{sample}` ({base_row['sample_label']})",
            f"- Variant: `{variant}`",
            f"- Category: `{base_row['case_category']}`",
            f"- Overlap class: `{base_row.get('locus_overlap_class', '')}`",
            f"- Paired vs TO concordance: `{same_row.get('paired_to_concordance', '')}`",
            "",
            "## Direct Observations",
            "",
        ]
        for row in mode_rows:
            if not row.get("present"):
                report_lines.append(f"- `{row['mode']}`: locus not present in this mode")
                continue
            report_lines.append(
                f"- `{row['mode']}`: truth=`{row['truth_label']}`, VC=`{row['verification_class']}`, "
                f"LOH_Subtype=`{row['loh_subtype']}`, core_LOH_like=`{row['core_loh_like']}`, "
                f"hp_ratio_core={format_metric(row['hp_ratio_core'])}, "
                f"effective_hp_reads={format_metric(row['effective_hp_reads'], integer=True)}, "
                f"hp0_ratio={format_metric(row['hp0_ratio'])}, "
                f"hp3_ratio={format_metric(row['hp3_ratio'])}, "
                f"AF={format_metric(row['caller_af'])}"
            )
        report_lines.extend(
            [
                "",
                "## Inference",
                "",
                "- 這裡只把 `hp_ratio_core`、`LOH_Subtype`、`VerificationClass`、`PS` 當作 evidence panel；不直接宣稱正式 LOH 結論。",
                "- 若 paired 與 TO 同位點結論一致，可視為跨 mode 支持；若不一致，優先檢查 AF、有效 HP read 數與 `HP0/HP3`。",
                "",
                "## Phase Block Limitation",
                "",
                f"- {describe_phase_status(mode_rows_by_name.get('to', {}))}",
                "",
                "## Figures",
                "",
                f"- Feature comparison: `{feature_plot}`" if feature_plot else "- Feature comparison: not generated",
                f"- HP tag counts: `{hp_plot}`" if hp_plot else "- HP tag counts: not generated",
                "",
                "## Evidence",
                "",
                f"- evidence_links.tsv: `{case_dir / 'evidence_links.tsv'}`",
            ]
        )
        (case_dir / "case_report.md").write_text("\n".join(report_lines), encoding="utf-8")
        registry_rows.append(
            {
                "case_id": case_id,
                "case_category": base_row["case_category"],
                "sample": sample,
                "sample_label": base_row["sample_label"],
                "variant_key": variant,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "locus_overlap_class": base_row.get("locus_overlap_class", ""),
                "paired_status": base_row.get("paired_status", ""),
                "to_status": base_row.get("to_status", ""),
                "paired_verification_class": base_row.get("paired_verification_class", ""),
                "to_verification_class": base_row.get("to_verification_class", ""),
                "paired_core_loh_like": base_row.get("core_loh_like_paired", False),
                "to_core_loh_like": base_row.get("core_loh_like_to", False),
                "paired_hp_ratio_core": base_row.get("hp_ratio_core_paired", math.nan),
                "to_hp_ratio_core": base_row.get("hp_ratio_core_to", math.nan),
                "paired_effective_hp_reads": base_row.get("effective_hp_reads_paired", math.nan),
                "to_effective_hp_reads": base_row.get("effective_hp_reads_to", math.nan),
                "case_dir": str(case_dir),
            }
        )
    registry_df = pd.DataFrame.from_records(registry_rows)
    write_tsv_rows(cases_dir / "case_registry.tsv", CASE_REGISTRY_FIELDS, registry_rows)
    return registry_df


def save_png(fig: plt.Figure, path: Path) -> None:
    ensure_dir(path.parent)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_fig01_loh_rates(summary_df: pd.DataFrame, out_path: Path) -> None:
    if summary_df.empty:
        return
    df = summary_df.copy()
    df["sample_mode"] = [sample_mode_label(s, m) for s, m in zip(df["sample_label"], df["mode"])]
    order = df.sort_values(["sample_label", "mode"])["sample_mode"].tolist()
    order = list(dict.fromkeys(order))
    x = np.arange(len(order))
    width = 0.38
    tp = df.set_index("sample_mode").reindex(order)["tp_core_loh_like_frac"].fillna(0).values
    fp = df.set_index("sample_mode").reindex(order)["fp_core_loh_like_frac"].fillna(0).values
    tp_n = df.set_index("sample_mode").reindex(order)["tp_core_loh_like_n"].fillna(0).astype(int).values
    fp_n = df.set_index("sample_mode").reindex(order)["fp_core_loh_like_n"].fillna(0).astype(int).values
    fig, ax = plt.subplots(figsize=(max(10, len(order) * 0.9), 5.2))
    ax.bar(x - width / 2, tp, width, label="TP", color="#1f77b4")
    ax.bar(x + width / 2, fp, width, label="FP", color="#d62728")
    for idx, (y1, y2, n1, n2) in enumerate(zip(tp, fp, tp_n, fp_n)):
        ax.text(idx - width / 2, y1 + 0.005, f"{n1}\n{y1:.2%}", ha="center", va="bottom", fontsize=7)
        ax.text(idx + width / 2, y2 + 0.005, f"{n2}\n{y2:.2%}", ha="center", va="bottom", fontsize=7)
    ax.set_ylabel("core LOH-like fraction")
    ax.set_title("Fig01 TP/FP LOH-like fractions by sample and mode")
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=45, ha="right")
    ax.legend()
    ax.set_ylim(0, max(0.25, np.nanmax(np.concatenate([tp, fp])) * 1.25))
    fig.tight_layout()
    save_png(fig, out_path)


def plot_fig02_hp_ratio_hist(all_df: pd.DataFrame, out_path: Path) -> None:
    if all_df.empty:
        return
    pairs = [(sample, mode) for sample in sorted(all_df["sample_label"].unique()) for mode in ["paired", "to"]]
    valid_pairs = [(s, m) for s, m in pairs if not all_df[(all_df["sample_label"] == s) & (all_df["mode"] == m)].empty]
    n = len(valid_pairs)
    cols = 2
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(14, max(3 * rows, 4)), squeeze=False)
    bins = np.linspace(0, 1, 21)
    for ax in axes.flatten():
        ax.axis("off")
    for idx, (sample_label, mode) in enumerate(valid_pairs):
        ax = axes[idx // cols][idx % cols]
        ax.axis("on")
        subset = all_df[(all_df["sample_label"] == sample_label) & (all_df["mode"] == mode) & (all_df["effective_hp_reads"] > 0)]
        tp = subset[subset["truth_label"] == "TP"]["hp_ratio_core"].dropna().values
        fp = subset[subset["truth_label"] == "FP"]["hp_ratio_core"].dropna().values
        ax.hist(tp, bins=bins, alpha=0.55, label=f"TP n={len(tp)}", color="#1f77b4", density=True)
        ax.hist(fp, bins=bins, alpha=0.55, label=f"FP n={len(fp)}", color="#d62728", density=True)
        ax.set_title(f"{sample_label} {mode}")
        ax.set_xlabel("hp_ratio_core")
        ax.set_ylabel("density")
        ax.legend(fontsize=7)
    fig.suptitle("Fig02 hp_ratio_core distributions (effective HP reads > 0)", y=1.01)
    fig.tight_layout()
    save_png(fig, out_path)


def plot_fig03_effective_hp_scatter(all_df: pd.DataFrame, out_path: Path) -> None:
    if all_df.empty:
        return
    pairs = [(sample, mode) for sample in sorted(all_df["sample_label"].unique()) for mode in ["paired", "to"]]
    valid_pairs = [(s, m) for s, m in pairs if not all_df[(all_df["sample_label"] == s) & (all_df["mode"] == m)].empty]
    cols = 2
    rows = math.ceil(len(valid_pairs) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(14, max(3.2 * rows, 4)), squeeze=False)
    for ax in axes.flatten():
        ax.axis("off")
    colors = {"TP": "#1f77b4", "FP": "#d62728"}
    for idx, (sample_label, mode) in enumerate(valid_pairs):
        ax = axes[idx // cols][idx % cols]
        ax.axis("on")
        subset = all_df[(all_df["sample_label"] == sample_label) & (all_df["mode"] == mode)]
        for truth_label in ["TP", "FP"]:
            part = subset[subset["truth_label"] == truth_label]
            if part.empty:
                continue
            sizes = np.clip(part["num_reads"].fillna(0).values, 5, 120)
            ax.scatter(
                part["effective_hp_reads"],
                part["hp_ratio_core"],
                s=sizes,
                alpha=0.25,
                color=colors[truth_label],
                label=truth_label,
                rasterized=True,
            )
        ax.set_title(f"{sample_label} {mode}")
        ax.set_xlabel("effective_hp_reads")
        ax.set_ylabel("hp_ratio_core")
        ax.set_ylim(-0.02, 1.02)
        ax.legend(fontsize=7)
    fig.suptitle("Fig03 effective HP reads vs hp_ratio_core", y=1.01)
    fig.tight_layout()
    save_png(fig, out_path)


def plot_fig04_verification_loh(verification_df: pd.DataFrame, out_path: Path) -> None:
    if verification_df.empty:
        return
    df = verification_df.copy()
    df["sample_mode_truth"] = [
        f"{sample_mode_label(sample_label, mode)}|{truth}"
        for sample_label, mode, truth in zip(df["sample_label"], df["mode"], df["truth_label"])
    ]
    df["combo"] = df["verification_class"] + "|" + df["loh_subtype"]
    pivot = (
        df.pivot_table(
            index="sample_mode_truth",
            columns="combo",
            values="count",
            aggfunc="sum",
            fill_value=0,
        )
        .sort_index()
    )
    row_sums = pivot.sum(axis=1).replace(0, np.nan)
    frac = pivot.div(row_sums, axis=0)
    fig, ax = plt.subplots(figsize=(max(12, 0.55 * len(frac.index)), 6))
    bottom = np.zeros(len(frac.index))
    colors = plt.cm.tab20(np.linspace(0, 1, max(1, len(frac.columns))))
    for color, column in zip(colors, frac.columns):
        values = frac[column].fillna(0).values
        ax.bar(np.arange(len(frac.index)), values, bottom=bottom, label=column, color=color)
        bottom += values
    ax.set_xticks(np.arange(len(frac.index)))
    ax.set_xticklabels(frac.index.tolist(), rotation=60, ha="right")
    ax.set_ylabel("fraction within sample/mode/truth")
    ax.set_title("Fig04 VerificationClass x LOH_Subtype structure")
    ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    save_png(fig, out_path)


def plot_fig05_feature_boxplots(all_df: pd.DataFrame, out_path: Path) -> None:
    if all_df.empty:
        return
    df = all_df.copy()
    df["group_label"] = [
        f"{mode}|{truth}|{'LOH' if loh else 'non-LOH'}"
        for mode, truth, loh in zip(df["mode"], df["truth_label"], df["core_loh_like"])
    ]
    group_order = [
        "paired|TP|LOH",
        "paired|TP|non-LOH",
        "paired|FP|LOH",
        "paired|FP|non-LOH",
        "to|TP|LOH",
        "to|TP|non-LOH",
        "to|FP|LOH",
        "to|FP|non-LOH",
    ]
    features = [
        ("allele_delta", "AlleleDelta"),
        ("pairwise_median_dist", "PairwiseMedianDist"),
        ("quality_score", "Quality_Score"),
        ("hp_assign_rate", "hp_assign_rate"),
        ("hp0_ratio", "hp0_ratio"),
        ("hp3_ratio", "hp3_ratio"),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(16, 11), squeeze=False)
    for ax, (column, title) in zip(axes.flatten(), features):
        data = [df[df["group_label"] == group][column].dropna().values for group in group_order]
        ax.boxplot(data, labels=group_order, showfliers=False)
        ax.set_title(title)
        ax.tick_params(axis="x", rotation=40)
    fig.suptitle("Fig05 LOH / non-LOH feature distributions", y=1.01)
    fig.tight_layout()
    save_png(fig, out_path)


def plot_fig06_heatmap(heatmap_df: pd.DataFrame, out_path: Path) -> None:
    if heatmap_df.empty:
        return
    metrics = ["hp_ratio_core", "caller_af", "effective_hp_reads"]
    fig, axes = plt.subplots(1, 3, figsize=(18, 7), squeeze=False)
    for idx, metric in enumerate(metrics):
        subset = heatmap_df[heatmap_df["metric"] == metric].copy()
        subset["sample_mode"] = [
            sample_mode_label(sample_label, mode)
            for sample_label, mode in zip(subset["sample_label"], subset["mode"])
        ]
        pivot = subset.pivot_table(index="sample_mode", columns="bin_label", values="fraction", aggfunc="sum", fill_value=0)
        pivot = pivot.sort_index()
        ax = axes[0][idx]
        im = ax.imshow(pivot.values, aspect="auto", cmap="magma", vmin=0, vmax=max(0.2, float(np.nanmax(pivot.values))))
        ax.set_title(metric)
        ax.set_xticks(np.arange(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns.tolist(), rotation=45, ha="right")
        ax.set_yticks(np.arange(len(pivot.index)))
        ax.set_yticklabels(pivot.index.tolist())
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("Fig06 sample x bin heatmaps", y=1.02)
    fig.tight_layout()
    save_png(fig, out_path)


def top_rows(df: pd.DataFrame, columns: Sequence[str], n: int = 6) -> str:
    if df.empty:
        return "_No rows_"
    display = df.loc[:, list(columns)].head(n).copy()
    return markdown_table(display.columns.tolist(), display.to_dict("records"))


def build_decision_ledger(output_dir: Path, generated_at: str, availability_missing: int, ps_fraction: float, loh_enrichment_df: pd.DataFrame) -> pd.DataFrame:
    to_rows = loh_enrichment_df[loh_enrichment_df["mode"] == "to"]
    paired_rows = loh_enrichment_df[loh_enrichment_df["mode"] == "paired"]
    max_to = to_rows["fp_core_loh_like_frac"].max() if not to_rows.empty else math.nan
    max_paired = paired_rows["fp_core_loh_like_frac"].max() if not paired_rows.empty else math.nan
    rows = [
        {
            "decision_id": "A",
            "decision_scope_round1": "sample_scope",
            "decision_topic": "Round 1 sample scope",
            "decision": "納入 7 個 paired 樣本與目前所有可用 TO 甲基 pilot",
            "status": "implemented",
            "rationale": f"manifest missing datasets={availability_missing}",
            "evidence_table": str(output_dir / "input_manifest.tsv"),
            "evidence_figure": "",
            "generated_at": generated_at,
        },
        {
            "decision_id": "B",
            "decision_scope_round1": "threshold",
            "decision_topic": "LOH-like threshold",
            "decision": "先保留 core hp_ratio full bins；<0.1/>0.9 只作 legacy reference",
            "status": "implemented",
            "rationale": "Round 1 先盤清分佈，再決定是否有更合理 cut",
            "evidence_table": str(output_dir / "loh_distribution_summary.tsv"),
            "evidence_figure": str(output_dir / "figures" / "fig02_hp_ratio_core_distribution.png"),
            "generated_at": generated_at,
        },
        {
            "decision_id": "C",
            "decision_scope_round1": "ps_strategy",
            "decision_topic": "PS handling",
            "decision": "paired 先不補 C++ export；TO 用 phased VCF 補查 PS；paired 僅在 tag.bam case audit 看 read-level PS",
            "status": "implemented",
            "rationale": f"TO phased-VCF PS non-empty fraction={ps_fraction:.3f}",
            "evidence_table": str(output_dir / "same_locus_compare.tsv"),
            "evidence_figure": "",
            "generated_at": generated_at,
        },
        {
            "decision_id": "D",
            "decision_scope_round1": "same_locus",
            "decision_topic": "paired-vs-TO compare",
            "decision": "same-locus key 使用 chr:pos:ref:alt，並固定輸出 overlap / concordance",
            "status": "implemented",
            "rationale": "先用 canonical/research_round 既有 split VCF 正規化後 key 做對齊",
            "evidence_table": str(output_dir / "same_locus_compare.tsv"),
            "evidence_figure": "",
            "generated_at": generated_at,
        },
        {
            "decision_id": "E",
            "decision_scope_round1": "first_observation",
            "decision_topic": "LOH observation strength",
            "decision": "LOH 更適合作為 observation / diagnostics layer，不直接當主證明",
            "status": "observed",
            "rationale": f"max TO FP LOH-like frac={max_to:.3f} ; max paired FP LOH-like frac={max_paired:.3f}",
            "evidence_table": str(output_dir / "loh_enrichment_by_sample_mode.tsv"),
            "evidence_figure": str(output_dir / "figures" / "fig01_loh_like_fraction_overview.png"),
            "generated_at": generated_at,
        },
    ]
    df = pd.DataFrame.from_records(rows)
    write_tsv_rows(output_dir / "decision_ledger.tsv", DECISION_FIELDS, rows)
    return df


def write_decision_notes(output_dir: Path, generated_at: str) -> None:
    notes = [
        f"# Round 1 decision notes ({generated_at})",
        "",
        "1. 這輪先盤點所有樣本 paired / TO 的分佈，而不是先追求單一 cut 或單一主張。",
        "2. `HP_Ratio` 不直接拿來當主判斷，因為 summary 中存在 `effective_hp_reads=0` 仍給 `HP_Ratio=0.5` 的情況。",
        "3. `hp_ratio_core = HP1FamilyN / (HP1FamilyN + HP2FamilyN)` 為本輪核心定義；`Potential_LOH` 只當 tool reference。",
        "4. TO 的 `PS` 可從 phased VCF 補查；paired 沒有 variant-level PS，僅能從 tag.bam read tags 看 local block evidence。",
        "5. 不同 PS block 的 HP1/HP2 不能外推出父母源，也不能直接跨 block 串聯成同一條長程 haplotype 故事。",
    ]
    (output_dir / "decision_notes.md").write_text("\n".join(notes), encoding="utf-8")


def write_open_questions(output_dir: Path, manifest_df: pd.DataFrame, all_df: pd.DataFrame, same_df: pd.DataFrame) -> None:
    missing_rows = manifest_df[manifest_df["availability_status"] != "ready"]
    ps_nonempty_frac = float(
        (
            all_df[(all_df["mode"] == "to")]["to_phase_ps"].fillna("").astype(str).str.len() > 0
        ).mean()
    ) if not all_df.empty else math.nan
    threshold_note = "LOH-like bins 在多數樣本呈現連續分佈，尚未支持單一 universal cut。" \
        if not all_df.empty else "尚未取得資料"
    overlap_note = ""
    if not same_df.empty:
        overlap_frac = float((same_df["locus_overlap_class"] == "both").mean())
        overlap_note = f"- same-locus overlap fraction = {overlap_frac:.3f}; 需注意 paired-only / TO-only 也可能主導觀察。"
    lines = [
        "# Open Questions",
        "",
        "## blocked by data availability",
        "",
    ]
    if missing_rows.empty:
        lines.append("- 目前 Round 1 manifest 內資料都可讀；沒有 sample 因缺檔被排除。")
    else:
        for row in missing_rows.to_dict("records"):
            lines.append(f"- {row['dataset_id']}: missing `{row['missing_components']}`")
    lines.extend(
        [
            "",
            "## blocked by PS",
            "",
            f"- TO phased VCF 的 non-empty PS fraction = {ps_nonempty_frac:.3f}",
            "- paired 目前沒有 variant-level PS export；只能用 tag.bam read-level PS 輔助解釋。",
            "- 不同 PS block 的 HP1/HP2 無法直接映射成父母源，也不能跨 block 直接串聯。",
            "",
            "## blocked by unclear threshold",
            "",
            f"- {threshold_note}",
            "- `Potential_LOH` 與自算 `core_loh_like` 是否要長期並存，需等更多 cross-sample round 再決定。",
            overlap_note,
            "",
            "## needs user confirmation",
            "",
            "- 若後續要把 phase-aware 敘事升級為正式論證，是否要先補 paired VCF / summary 的 PS export。",
            "- 若後續要做 gene-level second-hit 整合，是否優先把 paired-vs-TO overlap loci 擴展到 CNV/LOH 區段層。",
        ]
    )
    (output_dir / "open_questions.md").write_text("\n".join(lines), encoding="utf-8")


def write_round_summary(
    output_dir: Path,
    manifest_df: pd.DataFrame,
    benchmark_df: pd.DataFrame,
    loh_enrichment_df: pd.DataFrame,
    hp_qc_df: pd.DataFrame,
    same_df: pd.DataFrame,
    case_registry_df: pd.DataFrame,
    all_region_rows_count: int,
    generated_at: str,
) -> None:
    intersubmod_rows = benchmark_df[benchmark_df["method"] == "InterSubMod"].copy()
    mode_f1 = (
        intersubmod_rows.groupby("mode_group")
        .agg(
            mean_f1=("f1", "mean"),
            min_f1=("f1", "min"),
            max_f1=("f1", "max"),
            n=("sample_label", "size"),
        )
        .reset_index()
    )
    overlap_summary = (
        same_df.groupby(["sample_label", "paired_to_concordance"])
        .size()
        .rename("count")
        .reset_index()
        .sort_values(["sample_label", "count"], ascending=[True, False])
    ) if not same_df.empty else pd.DataFrame()
    lines = [
        f"# Round 1 LOH cross-sample audit summary ({generated_at})",
        "",
        "## Scope",
        "",
        f"- datasets in manifest: {len(manifest_df)}",
        f"- ready datasets: {(manifest_df['availability_status'] == 'ready').sum()}",
        f"- all region rows: {all_region_rows_count}",
        f"- same-locus compare rows: {len(same_df)}",
        f"- cases built: {len(case_registry_df)}",
        "",
        "## Key observations",
        "",
    ]
    if not loh_enrichment_df.empty:
        top_fp = loh_enrichment_df.sort_values("fp_core_loh_like_frac", ascending=False).head(5)
        lines.extend(
            [
                "- `core_loh_like` 在不同 sample/mode 間差異很大，不能先假設單一 cut 具有全域一致性。",
                "- paired 與 TO 都有 LOH-like 位點，但 TP/FP 的占比與 enrich 程度不一致，需拆樣本與 mode 解讀。",
                "- `Potential_LOH` 仍適合當 tool reference；Round 1 主分析使用自算 `hp_ratio_core`。",
                "",
                "### Highest FP LOH-like fractions",
                "",
                top_rows(top_fp, ["sample_label", "mode", "fp_core_loh_like_frac", "tp_core_loh_like_frac", "core_loh_fp_vs_tp_enrichment"]),
                "",
            ]
        )
    if not hp_qc_df.empty:
        lines.extend(
            [
                "### HP/phase QC snapshot",
                "",
                top_rows(
                    hp_qc_df.sort_values(["mode", "hp0_ge_0_5_fraction"], ascending=[True, False]),
                    ["sample_label", "mode", "truth_label", "median_hp_assign_rate", "median_hp0_ratio", "median_hp3_ratio", "effective_hp_eq0_fraction"],
                ),
                "",
            ]
        )
    if not overlap_summary.empty:
        lines.extend(
            [
                "### same-locus paired vs TO",
                "",
                top_rows(overlap_summary, ["sample_label", "paired_to_concordance", "count"], n=12),
                "",
            ]
        )
    if not mode_f1.empty:
        lines.extend(
            [
                "### InterSubMod F1 snapshot from benchmark_comparison.tsv",
                "",
                markdown_table(mode_f1.columns.tolist(), mode_f1.to_dict("records")),
                "",
            ]
        )
    lines.extend(
        [
            "## Deliverables",
            "",
            f"- input manifest: `{output_dir / 'input_manifest.tsv'}`",
            f"- all region table: `{output_dir / 'all_region_rows.tsv.gz'}`",
            f"- same locus compare: `{output_dir / 'same_locus_compare.tsv'}`",
            f"- figures: `{output_dir / 'figures'}`",
            f"- cases: `{output_dir / 'cases'}`",
        ]
    )
    (output_dir / "round_summary.md").write_text("\n".join(lines), encoding="utf-8")


def write_round_context(
    output_dir: Path,
    configs: Sequence[DatasetConfig],
    generated_at: str,
    verification_receipt: Optional[Dict[str, object]] = None,
) -> None:
    payload = {
        "generated_at": generated_at,
        "workspace_type": "observation_workspace",
        "workspace_goal": "Round 1 LOH cross-sample audit",
        "script": str(Path(__file__).resolve()),
        "dataset_ids": [config.dataset_id for config in configs],
        "samples": sorted({config.sample for config in configs}),
        "modes": sorted({config.mode for config in configs}),
        "verification_contract": verification_receipt or {"status": "PENDING_INPUT_VALIDATION"},
        "notes": [
            "LOH round is observation-first, not intervention-first.",
            "core hp ratio is recomputed from HP1FamilyN and HP2FamilyN.",
            "PS-aware narrative is partial because paired lacks variant-level PS export.",
        ],
    }
    write_json(output_dir / "round_context.json", payload)


def build_workspace(
    output_dir: Path,
    configs: Sequence[DatasetConfig],
    allow_unversioned_v1: bool = False,
) -> None:
    generated_at = generated_at_string()
    ensure_dir(output_dir)
    ensure_dir(output_dir / "figures")
    print(f"[1/9] initialize workspace: {output_dir}")
    write_round_context(output_dir, configs, generated_at)

    manifest_rows = build_manifest_rows(configs)
    manifest_df = pd.DataFrame.from_records(manifest_rows)
    write_tsv_rows(output_dir / "input_manifest.tsv", MANIFEST_FIELDS, manifest_rows)
    print(f"[2/9] manifest written: datasets={len(manifest_df)} ready={(manifest_df['availability_status'] == 'ready').sum()}")

    ready_configs = [config for config in configs if availability_status(config)["availability_status"] == "ready"]
    all_frames = []
    for idx, config in enumerate(ready_configs, start=1):
        print(f"[3/9] load dataset {idx}/{len(ready_configs)}: {config.dataset_id}")
        all_frames.append(
            load_dataset_rows(config, allow_unversioned_v1=allow_unversioned_v1)
        )
    all_frames = [frame for frame in all_frames if not frame.empty]
    if not all_frames:
        raise RuntimeError("No ready datasets produced region rows.")
    all_df = pd.concat(all_frames, ignore_index=True)
    provenance_statuses = sorted(
        all_df["VerificationProvenanceStatus"].dropna().astype(str).unique().tolist()
    )
    if len(provenance_statuses) != 1:
        raise SchemaContractError(
            f"LOH round provenance: mixed schema statuses across datasets: {provenance_statuses}"
        )
    verification_receipt = {
        "schema_status": provenance_statuses[0],
        "schema_versions": sorted(
            {
                value
                for value in all_df["VerificationSchemaVersion"].dropna().astype(str)
                if value.strip()
            }
        ),
        "current_selection_field": sorted(
            all_df["VerificationProvenanceSourceField"].dropna().astype(str).unique().tolist()
        ),
        "loh_selection_field": sorted(
            all_df["LOHProvenanceSourceField"].dropna().astype(str).unique().tolist()
        ),
        "allow_unversioned_v1": allow_unversioned_v1,
        "unknown_current_counts": sorted(
            all_df["VerificationUnknownCurrentCounts"].dropna().astype(str).unique().tolist()
        ),
        "provenance_fields": list(VERIFICATION_PROVENANCE_COLUMNS),
    }
    write_json(output_dir / "verification_schema_contract.json", verification_receipt)
    write_round_context(output_dir, configs, generated_at, verification_receipt)
    all_df = all_df.sort_values(["sample_label", "mode", "truth_label", "Chr", "Pos"]).reset_index(drop=True)
    all_region_path = output_dir / "all_region_rows.tsv.gz"
    all_df.to_csv(all_region_path, sep="\t", index=False, compression="gzip")
    print(f"[4/9] all region rows written: {len(all_df)}")

    benchmark_df = load_benchmark_summary(ready_configs)
    benchmark_df.to_csv(output_dir / "dataset_benchmark_summary.tsv", sep="\t", index=False)

    print("[5/9] build summary tables")
    loh_distribution_df = summarize_loh_distribution(all_df, all_region_path, generated_at)
    loh_distribution_df.to_csv(output_dir / "loh_distribution_summary.tsv", sep="\t", index=False)

    loh_enrichment_df = summarize_loh_enrichment(all_df, all_region_path, generated_at)
    loh_enrichment_df.to_csv(output_dir / "loh_enrichment_by_sample_mode.tsv", sep="\t", index=False)

    hp_balance_df = summarize_hp_balance(all_df, all_region_path, generated_at)
    hp_balance_df.to_csv(output_dir / "hp_family_balance_summary.tsv", sep="\t", index=False)

    hp_qc_df = summarize_hp_qc(all_df, all_region_path, generated_at)
    hp_qc_df.to_csv(output_dir / "hp_coverage_qc_summary.tsv", sep="\t", index=False)

    loh_feature_df = summarize_loh_vs_features(all_df, all_region_path, generated_at)
    loh_feature_df.to_csv(output_dir / "loh_vs_methyl_feature_summary.tsv", sep="\t", index=False)

    verification_df = summarize_verification_by_loh(all_df, all_region_path, generated_at)
    verification_df.to_csv(output_dir / "verificationclass_by_loh_subtype.tsv", sep="\t", index=False)

    heatmap_df = build_heatmap_table(all_df, all_region_path, generated_at)
    heatmap_df.to_csv(output_dir / "sample_bin_heatmap.tsv", sep="\t", index=False)

    print("[6/9] build same-locus compare")
    same_df = build_same_locus_compare(all_df)
    same_df.to_csv(output_dir / "same_locus_compare.tsv", sep="\t", index=False)

    print("[7/9] build case panel")
    case_registry_df = build_case_panel(all_df, same_df, output_dir)

    print("[8/9] render figures")
    plot_fig01_loh_rates(loh_enrichment_df, output_dir / "figures" / "fig01_loh_like_fraction_overview.png")
    plot_fig02_hp_ratio_hist(all_df, output_dir / "figures" / "fig02_hp_ratio_core_distribution.png")
    plot_fig03_effective_hp_scatter(all_df, output_dir / "figures" / "fig03_effective_hp_vs_hp_ratio_scatter.png")
    plot_fig04_verification_loh(verification_df, output_dir / "figures" / "fig04_verificationclass_lohsubtype_structure.png")
    plot_fig05_feature_boxplots(all_df, output_dir / "figures" / "fig05_loh_nonloh_feature_boxplots.png")
    plot_fig06_heatmap(heatmap_df, output_dir / "figures" / "fig06_sample_bin_heatmap.png")

    ps_fraction = float(
        (
            all_df[(all_df["mode"] == "to")]["to_phase_ps"].fillna("").astype(str).str.len() > 0
        ).mean()
    ) if not all_df.empty else math.nan
    build_decision_ledger(
        output_dir,
        generated_at,
        int((manifest_df["availability_status"] != "ready").sum()),
        ps_fraction,
        loh_enrichment_df,
    )
    write_decision_notes(output_dir, generated_at)
    write_open_questions(output_dir, manifest_df, all_df, same_df)
    write_round_summary(
        output_dir,
        manifest_df,
        benchmark_df,
        loh_enrichment_df,
        hp_qc_df,
        same_df,
        case_registry_df,
        len(all_df),
        generated_at,
    )
    print("[9/9] workspace complete")


def main() -> None:
    args = parse_args()
    selected_samples = set(args.sample or [])
    configs = [config for config in DATASET_CONFIGS if not selected_samples or config.sample in selected_samples or config.sample_label in selected_samples]
    output_dir = Path(args.output_dir)
    build_workspace(
        output_dir,
        configs,
        allow_unversioned_v1=args.allow_unversioned_v1,
    )


if __name__ == "__main__":
    main()
