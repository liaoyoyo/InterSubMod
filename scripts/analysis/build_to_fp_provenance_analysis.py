#!/usr/bin/env python3
"""Build TO FP provenance analysis workspace for HCC1395 5kHz / DORADO."""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import pandas as pd
import pysam

from research_common import compute_metrics, ensure_dir, infer_platform, load_tsv_rows, read_json, to_float, to_int, write_tsv_rows

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (
    SchemaContractError,
    read_evidence,
    select_current_view,
)


DEFAULT_SAMPLES = ["HCC1395", "HCC1395_DORADO"]

SAMPLE_CONFIGS = {
    "HCC1395": {
        "to_round_dir": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_to_pilot",
        "paired_dir": "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix",
    },
    "HCC1395_DORADO": {
        "to_round_dir": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_dorado_to_pilot",
        "paired_dir": "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395_DORADO/paired_full/20260315_HCC1395_DORADO_paired_full_full_complete_matrix",
    },
}

MASTER_FIELDS = [
    "sample",
    "platform",
    "variant_key",
    "chrom",
    "pos",
    "ref",
    "alt",
    "truth_status",
    "to_caller_filter",
    "to_caller_status",
    "to_pon_hit",
    "to_pon_hit_count",
    "to_verdict_somatic",
    "to_verdict_subclonal",
    "to_verdict_germline",
    "to_has_h_flag",
    "to_multihap_flag",
    "to_noancestry_flag",
    "to_qual",
    "to_gq",
    "to_dp",
    "to_af",
    "to_ad_ref",
    "to_ad_alt",
    "primary_class",
    "primary_detail",
    "to_baseline_pass",
    "to_longphase_status",
    "to_postprocess_status",
    "paired_resolution_stage",
    "paired_resolution_reason",
    "paired_normal_resolvable",
    "paired_raw_found",
    "paired_raw_filter",
    "paired_raw_qual",
    "paired_raw_gq",
    "paired_raw_dp",
    "paired_raw_af",
    "paired_raw_ad_ref",
    "paired_raw_ad_alt",
    "paired_raw_naf",
    "paired_raw_ndp",
    "paired_raw_nad_ref",
    "paired_raw_nad_alt",
    "paired_longphase_status",
    "paired_postprocess_status",
    "paired_final_status",
]

SAMPLE_SUMMARY_FIELDS = [
    "sample",
    "platform",
    "truth_total",
    "to_raw_fp_count",
    "caller_pon_filtered",
    "caller_nonpon_filtered",
    "longphase_to_removed",
    "to_postprocess_removed",
    "to_residual_final_fp",
    "paired_resolved_count",
    "paired_persistent_count",
    "to_baseline_tp",
    "to_baseline_fp",
    "to_baseline_fn",
    "to_baseline_precision",
    "to_baseline_recall",
    "to_baseline_f1",
    "to_final_tp",
    "to_final_fp",
    "to_final_fn",
    "to_final_precision",
    "to_final_recall",
    "to_final_f1",
    "paired_longphase_tp",
    "paired_longphase_fp",
    "paired_longphase_fn",
    "paired_longphase_precision",
    "paired_longphase_recall",
    "paired_longphase_f1",
    "paired_final_tp",
    "paired_final_fp",
    "paired_final_fn",
    "paired_final_precision",
    "paired_final_recall",
    "paired_final_f1",
    "to_fp_rule_trigger_count",
    "to_tp_rule_trigger_count",
    "paired_fp_rule_trigger_count",
    "paired_tp_rule_trigger_count",
    "to_round_dir",
    "paired_dir",
    "to_raw_vcf",
    "paired_raw_vcf",
]

CLASS_SUMMARY_FIELDS = [
    "sample",
    "primary_class",
    "primary_detail",
    "count",
    "fraction_of_to_raw_fp",
]

PAIRED_REASON_FIELDS = [
    "sample",
    "paired_resolution_stage",
    "paired_resolution_reason",
    "count",
    "fraction_of_to_raw_fp",
]

CLASS_PAIRED_CROSSTAB_FIELDS = [
    "sample",
    "primary_class",
    "paired_resolution_stage",
    "paired_resolution_reason",
    "count",
]

FEATURE_GROUP_FIELDS = [
    "sample",
    "group_name",
    "count",
    "median_af",
    "median_allele_delta",
    "median_pairwise_median_dist",
    "median_cramers_v",
    "median_quality_score",
    "median_hp_assign_rate",
    "median_allele_assign_rate",
    "top_verification_class",
    "verification_class_field",
    "verification_schema_status",
    "verification_rule_source",
    "verification_rule_selection_field",
    "top_agreement_type",
    "top_class_shift",
]

RULE_SWEEP_FIELDS = [
    "sample",
    "phase",
    "verification_rule_source",
    "rule_rank",
    "rule_id",
    "rule_label",
    "tp_removed",
    "fp_removed",
    "tp_total_before",
    "fp_total_before",
    "tp_after",
    "fp_after",
    "precision_after",
    "recall_after",
    "f1_after",
    "delta_f1_vs_final",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", action="append", default=None, help="Sample name, repeatable")
    parser.add_argument("--output-dir", required=True, help="Workspace output directory")
    parser.add_argument(
        "--verification-rule-source",
        required=True,
        choices=("current-bidirectional", "cluster-first-evidence"),
        help=(
            "Explicit meaning of the historical Strong-based candidate rule: "
            "canonical Strong_Bidirectional class or ClusterFirstSupport evidence"
        ),
    )
    return parser.parse_args()


def region_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{chrom}:{pos}:{ref}:{alt}"


def median_or_blank(values: Iterable[float]) -> str:
    clean = sorted(value for value in values if value is not None and not math.isnan(value))
    if not clean:
        return ""
    mid = len(clean) // 2
    if len(clean) % 2:
        return f"{clean[mid]:.6f}"
    return f"{(clean[mid - 1] + clean[mid]) / 2.0:.6f}"


def top_or_blank(values: Iterable[str]) -> str:
    counts = Counter(value for value in values if str(value).strip())
    if not counts:
        return ""
    value, _ = counts.most_common(1)[0]
    return value


def load_truth_keys(path: Path) -> set[str]:
    keys: set[str] = set()
    vf = pysam.VariantFile(str(path))
    for record in vf.fetch():
        if len(record.ref) != 1:
            continue
        if not record.alts or len(record.alts) != 1 or len(record.alts[0]) != 1:
            continue
        keys.add(region_key(record.chrom, record.pos, record.ref, record.alts[0]))
    vf.close()
    return keys


def load_bed_intervals(path: Path) -> dict[str, list[tuple[int, int]]]:
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            by_chrom[fields[0]].append((int(fields[1]), int(fields[2])))
    return by_chrom


def parse_filter_terms(record: pysam.VariantRecord) -> list[str]:
    terms = list(record.filter.keys())
    return terms or ["PASS"]


def parse_ad_field(value) -> tuple[int, int]:
    if value is None:
        return 0, 0
    if isinstance(value, (list, tuple)):
        parts = list(value)
    else:
        parts = str(value).split(",")
    if len(parts) < 2:
        return to_int(parts[0] if parts else 0), 0
    return to_int(parts[0]), to_int(parts[1])


def first_sample(record: pysam.VariantRecord):
    samples = list(record.samples.values())
    return samples[0] if samples else {}


def variant_payload_from_record(record: pysam.VariantRecord, paired: bool = False) -> dict[str, object]:
    sample = first_sample(record)
    filters = parse_filter_terms(record)
    ad_ref, ad_alt = parse_ad_field(sample.get("AD"))
    payload = {
        "filter_terms": filters,
        "filter": ";".join(filters),
        "qual": to_float(record.qual, 0.0),
        "gq": to_int(sample.get("GQ")),
        "dp": to_int(sample.get("DP")),
        "af": to_float(sample.get("AF")),
        "ad_ref": ad_ref,
        "ad_alt": ad_alt,
        "has_h_flag": "H" in record.info,
        "verdict_somatic": "Verdict_Somatic" in record.info,
        "verdict_subclonal": "Verdict_SubclonalSomatic" in record.info,
        "verdict_germline": "Verdict_Germline" in record.info,
    }
    if paired:
        nad_ref, nad_alt = parse_ad_field(sample.get("NAD"))
        payload.update(
            {
                "naf": to_float(sample.get("NAF")),
                "ndp": to_int(sample.get("NDP")),
                "nad_ref": nad_ref,
                "nad_alt": nad_alt,
            }
        )
    return payload


def iter_truth_scoped_records(vcf_path: Path, bed_intervals: dict[str, list[tuple[int, int]]]) -> Iterator[pysam.VariantRecord]:
    vf = pysam.VariantFile(str(vcf_path))
    try:
        for chrom in sorted(bed_intervals):
            if chrom not in vf.header.contigs:
                continue
            for start, end in bed_intervals[chrom]:
                for record in vf.fetch(chrom, start, end):
                    if len(record.ref) != 1:
                        continue
                    if not record.alts or len(record.alts) != 1 or len(record.alts[0]) != 1:
                        continue
                    yield record
    finally:
        vf.close()


def parse_vcf_keys(path: Path) -> set[str]:
    keys: set[str] = set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            keys.add(region_key(fields[0], int(fields[1]), fields[3], fields[4]))
    return keys


def parse_feature_vcf_records(path: Path, paired: bool = False) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            chrom, pos, _vid, ref, alt, qual, filt, _info = fields[:8]
            fmt = fields[8].split(":") if len(fields) >= 10 else []
            sample_vals = fields[9].split(":") if len(fields) >= 10 else []
            fmt_map = dict(zip(fmt, sample_vals))
            ad_ref, ad_alt = parse_ad_field(fmt_map.get("AD"))
            row = {
                "variant_key": region_key(chrom, int(pos), ref, alt),
                "qual": to_float(qual, 0.0),
                "filter": filt or "PASS",
                "gq": to_int(fmt_map.get("GQ")),
                "dp": to_int(fmt_map.get("DP")),
                "af": to_float(fmt_map.get("AF")),
                "ad_ref": ad_ref,
                "ad_alt": ad_alt,
            }
            if paired:
                nad_ref, nad_alt = parse_ad_field(fmt_map.get("NAD"))
                row.update(
                    {
                        "naf": to_float(fmt_map.get("NAF")),
                        "ndp": to_int(fmt_map.get("NDP")),
                        "nad_ref": nad_ref,
                        "nad_alt": nad_alt,
                    }
                )
            records[row["variant_key"]] = row
    return records


def attach_verification_rule_view(
    frame: pd.DataFrame,
    verification_rule_source: str,
) -> pd.DataFrame:
    """Attach the explicitly selected v2 class/evidence rule input."""
    current_view = select_current_view(frame)
    if current_view.unknown_counts:
        raise SchemaContractError(
            "TO FP provenance rule selection: unknown current VerificationClass values: "
            f"{current_view.unknown_counts}"
        )

    prepared = frame.copy()
    prepared["_verification_class_v2"] = current_view.values
    if verification_rule_source == "current-bidirectional":
        prepared["_verification_rule_match"] = current_view.values.eq("Strong_Bidirectional")
        selection_field = "VerificationClass"
    elif verification_rule_source == "cluster-first-evidence":
        evidence = read_evidence(frame)
        prepared["_verification_rule_match"] = evidence["ClusterFirstSupport"]
        prepared["_verification_evidence_path"] = evidence["EvidencePath"]
        selection_field = "ClusterFirstSupport"
    else:
        raise SchemaContractError(
            f"unsupported verification rule source: {verification_rule_source!r}"
        )

    prepared["_verification_rule_source"] = verification_rule_source
    prepared["_verification_rule_selection_field"] = selection_field
    prepared["_verification_schema_status"] = current_view.schema_status
    return prepared


def load_significance_map(
    path: Path,
    verification_rule_source: str,
) -> dict[str, dict[str, object]]:
    if not path.exists():
        return {}
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    prepared = attach_verification_rule_view(frame, verification_rule_source)
    rows: dict[str, dict[str, object]] = {}
    for row in prepared.to_dict(orient="records"):
        if row.get("RegionKey"):
            key = str(row["RegionKey"])
        elif row.get("region_key"):
            key = str(row["region_key"])
        else:
            chrom = row.get("Chr") or row.get("chrom")
            pos = row.get("Pos") or row.get("pos")
            ref = row.get("Ref") or row.get("ref")
            alt = row.get("Alt") or row.get("alt")
            if not all([chrom, pos, ref, alt]):
                continue
            key = region_key(str(chrom), int(str(pos)), str(ref), str(alt))
        rows[key] = row
    return rows


def load_label_agreement_map(path: Path, source_scope: str) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    if not path.exists():
        return rows
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("source_scope") != source_scope:
                continue
            key = row.get("region_key") or row.get("RegionKey")
            if not key:
                continue
            rows[key] = row
    return rows


def build_feature_map(
    vcf_records: dict[str, dict[str, object]],
    summary_map: dict[str, dict[str, object]],
    label_map: dict[str, dict[str, str]],
) -> dict[str, dict[str, object]]:
    rows: dict[str, dict[str, object]] = {}
    for key, vcf_row in vcf_records.items():
        summary = summary_map.get(key, {})
        label = label_map.get(key, {})
        row = {
            "variant_key": key,
            "af": to_float(vcf_row.get("af")),
            "qual": to_float(vcf_row.get("qual")),
            "gq": to_int(vcf_row.get("gq")),
            "dp": to_int(vcf_row.get("dp")),
            "ad_alt": to_int(vcf_row.get("ad_alt")),
            "ad_ref": to_int(vcf_row.get("ad_ref")),
            "pairwise_median_dist": to_float(summary.get("PairwiseMedianDist")),
            "cramers_v": to_float(summary.get("CramersV")),
            "quality_score": to_float(summary.get("Quality_Score")),
            "allele_delta": to_float(summary.get("AlleleDelta") or summary.get("AlleleDelta_abs")),
            "verification_class": summary.get("_verification_class_v2", ""),
            "verification_class_field": "VerificationClass" if summary else "",
            "verification_rule_match": summary.get("_verification_rule_match"),
            "verification_rule_source": summary.get("_verification_rule_source", ""),
            "verification_rule_selection_field": summary.get(
                "_verification_rule_selection_field", ""
            ),
            "verification_schema_status": summary.get("_verification_schema_status", ""),
            "quality_tier": summary.get("Quality_Tier", ""),
            "potential_loh": summary.get("Potential_LOH", ""),
            "coverage_category": summary.get("Coverage_Category", ""),
            "suggest_filter": summary.get("SuggestFilter", ""),
            "hp_assign_rate": to_float(label.get("hp_assign_rate")),
            "allele_assign_rate": to_float(label.get("allele_assign_rate")),
            "agreement_type": label.get("agreement_type", ""),
            "class_shift": label.get("class_shift", ""),
            "effect_size": to_float(label.get("effect_size")),
            "cluster_class": label.get("cluster_class", ""),
            "label_class": label.get("label_class", ""),
        }
        rows[key] = row
    return rows


def load_stage_feature_maps(
    base_dir: Path,
    tp_scope: str,
    fp_scope: str,
    verification_rule_source: str,
) -> tuple[dict[str, dict[str, object]], dict[str, dict[str, object]]]:
    label_path = base_dir / "label_cluster_agreement.tsv"
    tp_vcf = parse_feature_vcf_records(base_dir / tp_scope / "tp.vcf" if (base_dir / tp_scope / "tp.vcf").exists() else base_dir / tp_scope / "filtered_snv_tp.vcf.gz")
    fp_vcf = parse_feature_vcf_records(base_dir / fp_scope / "fp.vcf" if (base_dir / fp_scope / "fp.vcf").exists() else base_dir / fp_scope / "filtered_snv_fp.vcf.gz")
    tp_summary = load_significance_map(base_dir / "step05_intersubmod" / "intersubmod_tp" / "significance_summary.csv", verification_rule_source) if (base_dir / "step05_intersubmod").exists() else load_significance_map(base_dir / "intersubmod_tp" / "significance_summary.csv", verification_rule_source)
    fp_summary = load_significance_map(base_dir / "step05_intersubmod" / "intersubmod_fp" / "significance_summary.csv", verification_rule_source) if (base_dir / "step05_intersubmod").exists() else load_significance_map(base_dir / "intersubmod_fp" / "significance_summary.csv", verification_rule_source)
    tp_label = load_label_agreement_map(label_path, "tp")
    fp_label = load_label_agreement_map(label_path, "fp")
    return build_feature_map(tp_vcf, tp_summary, tp_label), build_feature_map(fp_vcf, fp_summary, fp_label)


def load_to_feature_maps(
    to_round_dir: Path,
    verification_rule_source: str,
) -> tuple[dict[str, dict[str, object]], dict[str, dict[str, object]]]:
    tp_vcf = parse_feature_vcf_records(to_round_dir / "step04_benchmark_longphase_to" / "tp.vcf")
    fp_vcf = parse_feature_vcf_records(to_round_dir / "step04_benchmark_longphase_to" / "fp.vcf")
    tp_summary = load_significance_map(to_round_dir / "step05_intersubmod" / "intersubmod_tp" / "significance_summary.csv", verification_rule_source)
    fp_summary = load_significance_map(to_round_dir / "step05_intersubmod" / "intersubmod_fp" / "significance_summary.csv", verification_rule_source)
    tp_label = load_label_agreement_map(to_round_dir / "label_cluster_agreement.tsv", "tp")
    fp_label = load_label_agreement_map(to_round_dir / "label_cluster_agreement.tsv", "fp")
    return build_feature_map(tp_vcf, tp_summary, tp_label), build_feature_map(fp_vcf, fp_summary, fp_label)


def load_paired_feature_maps(
    paired_dir: Path,
    verification_rule_source: str,
) -> tuple[dict[str, dict[str, object]], dict[str, dict[str, object]]]:
    tp_vcf = parse_feature_vcf_records(paired_dir / "longphase_s" / "tp.vcf", paired=True)
    fp_vcf = parse_feature_vcf_records(paired_dir / "longphase_s" / "fp.vcf", paired=True)
    tp_summary = load_significance_map(paired_dir / "intersubmod_tp" / "significance_summary.csv", verification_rule_source)
    fp_summary = load_significance_map(paired_dir / "intersubmod_fp" / "significance_summary.csv", verification_rule_source)
    tp_label = load_label_agreement_map(paired_dir / "label_cluster_agreement.tsv", "tp")
    fp_label = load_label_agreement_map(paired_dir / "label_cluster_agreement.tsv", "fp")
    return build_feature_map(tp_vcf, tp_summary, tp_label), build_feature_map(fp_vcf, fp_summary, fp_label)


def rule_trigger(feature: dict[str, object], rule: dict[str, object]) -> bool:
    if not rule.get("apply_filter", False):
        return False
    af = to_float(feature.get("af"))
    allele_delta = to_float(feature.get("allele_delta"))
    if math.isnan(af) or math.isnan(allele_delta):
        return False
    if af > to_float(rule.get("vaf_max"), default=math.inf):
        return False
    if allele_delta < to_float(rule.get("allele_delta_min"), default=0.0):
        return False
    if rule.get("require_cv_support", False):
        cramers_v = to_float(feature.get("cramers_v"))
        if math.isnan(cramers_v):
            return False
        if cramers_v > to_float(rule.get("cv_support_max"), default=math.inf):
            return False
    return True


def load_configs(sample: str) -> dict[str, object]:
    config = SAMPLE_CONFIGS[sample].copy()
    to_round_dir = Path(config["to_round_dir"]).resolve()
    paired_dir = Path(config["paired_dir"]).resolve()
    to_context = read_json(to_round_dir / "round_context.json")
    paired_context = read_json(paired_dir / "run_context.json")
    metrics_to = read_json(to_round_dir / "metrics.json")
    metrics_paired = read_json(paired_dir / "metrics.json")
    config.update(
        {
            "sample": sample,
            "platform": infer_platform(sample),
            "to_round_dir": to_round_dir,
            "paired_dir": paired_dir,
            "to_raw_vcf": Path(to_context["somatic_vcf"]).resolve(),
            "paired_raw_vcf": Path(paired_context["somatic_vcf"]).resolve(),
            "truth_vcf": Path(to_context["truth_vcf"]).resolve(),
            "truth_bed": Path(to_context["truth_bed"]).resolve(),
            "truth_total": int(to_context["truth_total"]),
            "metrics_to": metrics_to,
            "metrics_paired": metrics_paired,
        }
    )
    return config


def classify_paired_resolution(
    key: str,
    paired_raw_map: dict[str, dict[str, object]],
    paired_longphase_fp_keys: set[str],
    paired_post_removed_fp_keys: set[str],
) -> tuple[str, str, bool, str, str, str]:
    raw = paired_raw_map.get(key)
    if raw is None:
        return ("raw_absent", "paired_raw_absent", True, "absent", "not_applicable", "not_present_in_paired_raw")

    filters = set(raw["filter_terms"])
    germline_filter = "Germline" in filters
    normal_alt_support = to_float(raw.get("naf")) >= 0.05 and to_int(raw.get("nad_alt")) >= 3

    if germline_filter:
        return ("raw_filter", "germline_filter", True, "filtered", "not_applicable", "removed_in_paired_raw")
    if normal_alt_support:
        return ("raw_filter", "normal_alt_support", True, "filtered", "not_applicable", "removed_in_paired_raw")
    if "PASS" not in filters or len(filters) != 1:
        return ("raw_filter", "paired_raw_artifact_filtered", True, "filtered", "not_applicable", "removed_in_paired_raw")
    if key not in paired_longphase_fp_keys:
        return ("longphase_s", "paired_longphase_s_removed", True, "pass", "removed", "removed_before_paired_final")
    if key in paired_post_removed_fp_keys:
        return ("postprocess", "paired_postprocess_removed", True, "pass", "kept", "removed_in_paired_postprocess")
    return ("persistent", "paired_persistent_final_fp", False, "pass", "kept", "kept_in_paired_final")


def build_feature_group_summary(sample: str, group_name: str, rows: Sequence[dict[str, object]]) -> dict[str, object]:
    def unique_metadata(field: str) -> str:
        values = sorted({str(row.get(field)) for row in rows if str(row.get(field, "")).strip()})
        if len(values) > 1:
            raise SchemaContractError(
                f"feature group {group_name}: mixed {field} values: {values}"
            )
        return values[0] if values else ""

    return {
        "sample": sample,
        "group_name": group_name,
        "count": len(rows),
        "median_af": median_or_blank(to_float(row.get("af")) for row in rows),
        "median_allele_delta": median_or_blank(to_float(row.get("allele_delta")) for row in rows),
        "median_pairwise_median_dist": median_or_blank(to_float(row.get("pairwise_median_dist")) for row in rows),
        "median_cramers_v": median_or_blank(to_float(row.get("cramers_v")) for row in rows),
        "median_quality_score": median_or_blank(to_float(row.get("quality_score")) for row in rows),
        "median_hp_assign_rate": median_or_blank(to_float(row.get("hp_assign_rate")) for row in rows),
        "median_allele_assign_rate": median_or_blank(to_float(row.get("allele_assign_rate")) for row in rows),
        "top_verification_class": top_or_blank(row.get("verification_class", "") for row in rows),
        "verification_class_field": unique_metadata("verification_class_field"),
        "verification_schema_status": unique_metadata("verification_schema_status"),
        "verification_rule_source": unique_metadata("verification_rule_source"),
        "verification_rule_selection_field": unique_metadata(
            "verification_rule_selection_field"
        ),
        "top_agreement_type": top_or_blank(row.get("agreement_type", "") for row in rows),
        "top_class_shift": top_or_blank(row.get("class_shift", "") for row in rows),
    }


def apply_rule(rows: Sequence[dict[str, object]], predicate) -> tuple[int, list[str]]:
    matched = [row["variant_key"] for row in rows if predicate(row)]
    return len(matched), matched


def make_numeric_rule(rule_id: str, label: str, predicate):
    return {"rule_id": rule_id, "rule_label": label, "predicate": predicate}


def require_verification_rule_match(row: dict[str, object]) -> bool:
    value = row.get("verification_rule_match")
    if not isinstance(value, bool):
        raise SchemaContractError(
            "candidate rule requires a matched, validated significance-summary row; "
            f"variant={row.get('variant_key', '<UNKNOWN>')} value={value!r}"
        )
    return value


def build_candidate_rules(
    sample: str,
    verification_rule_source: str,
) -> list[dict[str, object]]:
    rules: list[dict[str, object]] = []
    for af_max, ad_min in itertools.product([0.03, 0.05, 0.08, 0.10, 0.12], [0.15, 0.20, 0.25, 0.30]):
        rules.append(
            make_numeric_rule(
                f"af_le_{af_max:.2f}_ad_ge_{ad_min:.2f}",
                f"AF<={af_max:.2f} and AlleleDelta>={ad_min:.2f}",
                lambda row, af_max=af_max, ad_min=ad_min: to_float(row.get("af")) <= af_max and to_float(row.get("allele_delta")) >= ad_min,
            )
        )
    for af_max, ad_min, cv_max in itertools.product([0.05, 0.08, 0.10], [0.15, 0.20, 0.25], [0.03, 0.05, 0.10]):
        rules.append(
            make_numeric_rule(
                f"af_le_{af_max:.2f}_ad_ge_{ad_min:.2f}_cv_le_{cv_max:.2f}",
                f"AF<={af_max:.2f} and AlleleDelta>={ad_min:.2f} and CramersV<={cv_max:.2f}",
                lambda row, af_max=af_max, ad_min=ad_min, cv_max=cv_max: to_float(row.get("af")) <= af_max and to_float(row.get("allele_delta")) >= ad_min and to_float(row.get("cramers_v")) <= cv_max,
            )
        )
    for af_max, ad_min, pairwise_max in itertools.product([0.05, 0.08, 0.10], [0.15, 0.20, 0.25], [0.05, 0.10, 0.15, 0.20]):
        rules.append(
            make_numeric_rule(
                f"af_le_{af_max:.2f}_ad_ge_{ad_min:.2f}_pairwise_le_{pairwise_max:.2f}",
                f"AF<={af_max:.2f} and AlleleDelta>={ad_min:.2f} and PairwiseMedianDist<={pairwise_max:.2f}",
                lambda row, af_max=af_max, ad_min=ad_min, pairwise_max=pairwise_max: to_float(row.get("af")) <= af_max and to_float(row.get("allele_delta")) >= ad_min and to_float(row.get("pairwise_median_dist")) <= pairwise_max,
            )
        )
    if verification_rule_source == "current-bidirectional":
        verification_rule = make_numeric_rule(
            "strong_bidirectional_lowaf_highad",
            "VerificationClass(v2)=Strong_Bidirectional and AF<=0.10 and AlleleDelta>=0.15",
            lambda row: require_verification_rule_match(row)
            and to_float(row.get("af")) <= 0.10
            and to_float(row.get("allele_delta")) >= 0.15,
        )
    elif verification_rule_source == "cluster-first-evidence":
        verification_rule = make_numeric_rule(
            "cluster_first_support_lowaf_highad",
            "ClusterFirstSupport=true and AF<=0.10 and AlleleDelta>=0.15",
            lambda row: require_verification_rule_match(row)
            and to_float(row.get("af")) <= 0.10
            and to_float(row.get("allele_delta")) >= 0.15,
        )
    else:
        raise SchemaContractError(
            f"unsupported verification rule source: {verification_rule_source!r}"
        )

    rules.extend(
        [
            verification_rule,
            make_numeric_rule(
                "cluster_plus_weak_lowaf_highad",
                "agreement=cluster_plus_weak_label and AF<=0.10 and AlleleDelta>=0.15",
                lambda row: row.get("agreement_type") == "cluster_plus_weak_label" and to_float(row.get("af")) <= 0.10 and to_float(row.get("allele_delta")) >= 0.15,
            ),
            make_numeric_rule(
                "strong_to_weak_lowaf_highad",
                "class_shift=Strong->Weak and AF<=0.10 and AlleleDelta>=0.15",
                lambda row: row.get("class_shift") == "Strong->Weak" and to_float(row.get("af")) <= 0.10 and to_float(row.get("allele_delta")) >= 0.15,
            ),
            make_numeric_rule(
                "quality_low_and_lowaf",
                "QualityScore<=20 and AF<=0.10",
                lambda row: to_float(row.get("quality_score")) <= 20 and to_float(row.get("af")) <= 0.10,
            ),
        ]
    )
    return rules


def run_rule_sweep(
    sample: str,
    tp_rows: Sequence[dict[str, object]],
    fp_rows: Sequence[dict[str, object]],
    final_metrics: dict[str, object],
    candidate_rules: Sequence[dict[str, object]],
    phase: str,
    verification_rule_source: str,
    rule_limit: Optional[int] = None,
) -> list[dict[str, object]]:
    total_tp = int(final_metrics["filtered"]["tp"])
    total_fp = int(final_metrics["filtered"]["fp"])
    truth_total = int(final_metrics["truth_total"])
    baseline_f1 = to_float(final_metrics["filtered"]["f1"])
    results: list[dict[str, object]] = []
    for idx, rule in enumerate(candidate_rules, start=1):
        tp_removed, _ = apply_rule(tp_rows, rule["predicate"])
        fp_removed, _ = apply_rule(fp_rows, rule["predicate"])
        metrics = compute_metrics(total_tp - tp_removed, total_fp - fp_removed, truth_total)
        results.append(
            {
                "sample": sample,
                "phase": phase,
                "verification_rule_source": verification_rule_source,
                "rule_rank": idx,
                "rule_id": rule["rule_id"],
                "rule_label": rule["rule_label"],
                "tp_removed": tp_removed,
                "fp_removed": fp_removed,
                "tp_total_before": total_tp,
                "fp_total_before": total_fp,
                "tp_after": metrics["tp"],
                "fp_after": metrics["fp"],
                "precision_after": f"{metrics['precision']:.6f}",
                "recall_after": f"{metrics['recall']:.6f}",
                "f1_after": f"{metrics['f1']:.6f}",
                "delta_f1_vs_final": f"{metrics['f1'] - baseline_f1:.6f}",
            }
        )
    results.sort(key=lambda row: (float(row["delta_f1_vs_final"]), row["fp_removed"] - row["tp_removed"]), reverse=True)
    if rule_limit is not None:
        return results[:rule_limit]
    return results


def write_gzip_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    ensure_dir(path.parent)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def analyze_sample(config: dict[str, object], output_dir: Path) -> dict[str, object]:
    sample = config["sample"]
    platform = config["platform"]
    verification_rule_source = str(config["verification_rule_source"])
    truth_keys = load_truth_keys(config["truth_vcf"])
    bed_intervals = load_bed_intervals(config["truth_bed"])

    to_baseline_tp_keys = parse_vcf_keys(config["to_round_dir"] / "step02_benchmark_clairs_to" / "tp.vcf")
    to_baseline_fp_keys = parse_vcf_keys(config["to_round_dir"] / "step02_benchmark_clairs_to" / "fp.vcf")
    to_longphase_tp_keys = parse_vcf_keys(config["to_round_dir"] / "step04_benchmark_longphase_to" / "tp.vcf")
    to_longphase_fp_keys = parse_vcf_keys(config["to_round_dir"] / "step04_benchmark_longphase_to" / "fp.vcf")

    paired_longphase_tp_keys = parse_vcf_keys(config["paired_dir"] / "longphase_s" / "tp.vcf")
    paired_longphase_fp_keys = parse_vcf_keys(config["paired_dir"] / "longphase_s" / "fp.vcf")

    to_tp_feature_map, to_fp_feature_map = load_to_feature_maps(
        config["to_round_dir"], verification_rule_source
    )
    paired_tp_feature_map, paired_fp_feature_map = load_paired_feature_maps(
        config["paired_dir"], verification_rule_source
    )

    to_rule = config["metrics_to"]["rule"]
    paired_rule = config["metrics_paired"]["rule"]

    to_post_removed_fp_keys = {key for key, row in to_fp_feature_map.items() if rule_trigger(row, to_rule)}
    to_post_removed_tp_keys = {key for key, row in to_tp_feature_map.items() if rule_trigger(row, to_rule)}
    paired_post_removed_fp_keys = {key for key, row in paired_fp_feature_map.items() if rule_trigger(row, paired_rule)}
    paired_post_removed_tp_keys = {key for key, row in paired_tp_feature_map.items() if rule_trigger(row, paired_rule)}

    paired_raw_map: dict[str, dict[str, object]] = {}
    for record in iter_truth_scoped_records(config["paired_raw_vcf"], bed_intervals):
        key = region_key(record.chrom, record.pos, record.ref, record.alts[0])
        row = variant_payload_from_record(record, paired=True)
        paired_raw_map[key] = row

    sample_rows: list[dict[str, object]] = []
    primary_counter: Counter[tuple[str, str]] = Counter()
    paired_reason_counter: Counter[tuple[str, str]] = Counter()
    class_paired_counter: Counter[tuple[str, str, str]] = Counter()

    master_path = output_dir / f"{sample.lower()}_to_fp_provenance_master.tsv.gz"
    with gzip.open(master_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MASTER_FIELDS, delimiter="\t")
        writer.writeheader()
        for record in iter_truth_scoped_records(config["to_raw_vcf"], bed_intervals):
            key = region_key(record.chrom, record.pos, record.ref, record.alts[0])
            truth_status = "TP" if key in truth_keys else "FP"
            if truth_status != "FP":
                continue
            to_row = variant_payload_from_record(record, paired=False)
            filters = set(to_row["filter_terms"])
            to_caller_pass = "PASS" in filters and len(filters) == 1
            pon_hit_count = sum(1 for tag in ("PoN_1", "PoN_2", "PoN_3", "PoN_4") if tag in record.info)
            pon_hit = pon_hit_count > 0 or "NonSomatic" in filters
            to_baseline_pass = key in to_baseline_fp_keys

            if not to_caller_pass:
                primary_class = "caller_pon_filtered" if pon_hit else "caller_nonpon_filtered"
                primary_detail = to_row["filter"]
                to_longphase_status = "not_applicable"
                to_postprocess_status = "not_applicable"
            else:
                if key not in to_longphase_fp_keys:
                    primary_class = "longphase_to_removed"
                    primary_detail = "removed_before_longphase_fp"
                    to_longphase_status = "removed"
                    to_postprocess_status = "not_applicable"
                else:
                    to_longphase_status = "kept"
                    if key in to_post_removed_fp_keys:
                        primary_class = "to_postprocess_removed"
                        primary_detail = to_rule["rule_id"]
                        to_postprocess_status = "removed"
                    else:
                        primary_class = "to_residual_final_fp"
                        primary_detail = "kept_after_to_final_rule"
                        to_postprocess_status = "kept"

            (
                paired_stage,
                paired_reason,
                paired_resolvable,
                paired_longphase_status,
                paired_postprocess_status,
                paired_final_status,
            ) = classify_paired_resolution(
                key,
                paired_raw_map,
                paired_longphase_fp_keys,
                paired_post_removed_fp_keys,
            )
            paired_raw = paired_raw_map.get(key, {})

            out_row = {
                "sample": sample,
                "platform": platform,
                "variant_key": key,
                "chrom": record.chrom,
                "pos": record.pos,
                "ref": record.ref,
                "alt": record.alts[0],
                "truth_status": truth_status,
                "to_caller_filter": to_row["filter"],
                "to_caller_status": "pass" if to_caller_pass else "filtered",
                "to_pon_hit": pon_hit,
                "to_pon_hit_count": pon_hit_count,
                "to_verdict_somatic": to_row["verdict_somatic"],
                "to_verdict_subclonal": to_row["verdict_subclonal"],
                "to_verdict_germline": to_row["verdict_germline"],
                "to_has_h_flag": to_row["has_h_flag"],
                "to_multihap_flag": "MultiHap" in filters,
                "to_noancestry_flag": "NoAncestry" in filters,
                "to_qual": f"{to_row['qual']:.6f}",
                "to_gq": to_row["gq"],
                "to_dp": to_row["dp"],
                "to_af": f"{to_row['af']:.6f}" if not math.isnan(to_row["af"]) else "",
                "to_ad_ref": to_row["ad_ref"],
                "to_ad_alt": to_row["ad_alt"],
                "primary_class": primary_class,
                "primary_detail": primary_detail,
                "to_baseline_pass": to_baseline_pass,
                "to_longphase_status": to_longphase_status,
                "to_postprocess_status": to_postprocess_status,
                "paired_resolution_stage": paired_stage,
                "paired_resolution_reason": paired_reason,
                "paired_normal_resolvable": paired_resolvable,
                "paired_raw_found": key in paired_raw_map,
                "paired_raw_filter": paired_raw.get("filter", ""),
                "paired_raw_qual": f"{to_float(paired_raw.get('qual')):.6f}" if paired_raw else "",
                "paired_raw_gq": paired_raw.get("gq", ""),
                "paired_raw_dp": paired_raw.get("dp", ""),
                "paired_raw_af": f"{to_float(paired_raw.get('af')):.6f}" if paired_raw and not math.isnan(to_float(paired_raw.get('af'))) else "",
                "paired_raw_ad_ref": paired_raw.get("ad_ref", ""),
                "paired_raw_ad_alt": paired_raw.get("ad_alt", ""),
                "paired_raw_naf": f"{to_float(paired_raw.get('naf')):.6f}" if paired_raw and not math.isnan(to_float(paired_raw.get('naf'))) else "",
                "paired_raw_ndp": paired_raw.get("ndp", ""),
                "paired_raw_nad_ref": paired_raw.get("nad_ref", ""),
                "paired_raw_nad_alt": paired_raw.get("nad_alt", ""),
                "paired_longphase_status": paired_longphase_status,
                "paired_postprocess_status": paired_postprocess_status,
                "paired_final_status": paired_final_status,
            }
            writer.writerow(out_row)
            sample_rows.append(out_row)
            primary_counter[(primary_class, primary_detail)] += 1
            paired_reason_counter[(paired_stage, paired_reason)] += 1
            class_paired_counter[(primary_class, paired_stage, paired_reason)] += 1

    to_fp_total = len(sample_rows)
    to_final_metrics = config["metrics_to"]["filtered"]
    to_baseline_metrics = config["metrics_to"]["baseline"]
    paired_longphase_metrics = config["metrics_paired"]["baseline"]
    paired_final_metrics = config["metrics_paired"]["filtered"]

    class_summary_rows = []
    for (primary_class, primary_detail), count in sorted(primary_counter.items()):
        class_summary_rows.append(
            {
                "sample": sample,
                "primary_class": primary_class,
                "primary_detail": primary_detail,
                "count": count,
                "fraction_of_to_raw_fp": f"{count / to_fp_total:.6f}",
            }
        )

    paired_summary_rows = []
    for (stage, reason), count in sorted(paired_reason_counter.items()):
        paired_summary_rows.append(
            {
                "sample": sample,
                "paired_resolution_stage": stage,
                "paired_resolution_reason": reason,
                "count": count,
                "fraction_of_to_raw_fp": f"{count / to_fp_total:.6f}",
            }
        )

    crosstab_rows = []
    for (primary_class, stage, reason), count in sorted(class_paired_counter.items()):
        crosstab_rows.append(
            {
                "sample": sample,
                "primary_class": primary_class,
                "paired_resolution_stage": stage,
                "paired_resolution_reason": reason,
                "count": count,
            }
        )

    final_tp_rows = [row for key, row in to_tp_feature_map.items() if key not in to_post_removed_tp_keys]
    final_fp_rows = [row for key, row in to_fp_feature_map.items() if key not in to_post_removed_fp_keys]
    residual_fp_paired_resolved = [row for row in final_fp_rows if classify_paired_resolution(row["variant_key"], paired_raw_map, paired_longphase_fp_keys, paired_post_removed_fp_keys)[2]]
    residual_fp_paired_persistent = [row for row in final_fp_rows if not classify_paired_resolution(row["variant_key"], paired_raw_map, paired_longphase_fp_keys, paired_post_removed_fp_keys)[2]]
    to_post_removed_fp_rows = [row for key, row in to_fp_feature_map.items() if key in to_post_removed_fp_keys]

    feature_group_rows = [
        build_feature_group_summary(sample, "to_final_tp", final_tp_rows),
        build_feature_group_summary(sample, "to_postprocess_removed_fp", to_post_removed_fp_rows),
        build_feature_group_summary(sample, "to_residual_fp_paired_resolved", residual_fp_paired_resolved),
        build_feature_group_summary(sample, "to_residual_fp_paired_persistent", residual_fp_paired_persistent),
    ]

    discovery_rules = build_candidate_rules(sample, verification_rule_source)
    sweep_rows = run_rule_sweep(
        sample,
        final_tp_rows,
        final_fp_rows,
        config["metrics_to"],
        discovery_rules,
        phase="discovery",
        verification_rule_source=verification_rule_source,
    )

    summary_row = {
        "sample": sample,
        "platform": platform,
        "truth_total": config["truth_total"],
        "to_raw_fp_count": to_fp_total,
        "caller_pon_filtered": sum(1 for row in sample_rows if row["primary_class"] == "caller_pon_filtered"),
        "caller_nonpon_filtered": sum(1 for row in sample_rows if row["primary_class"] == "caller_nonpon_filtered"),
        "longphase_to_removed": sum(1 for row in sample_rows if row["primary_class"] == "longphase_to_removed"),
        "to_postprocess_removed": sum(1 for row in sample_rows if row["primary_class"] == "to_postprocess_removed"),
        "to_residual_final_fp": sum(1 for row in sample_rows if row["primary_class"] == "to_residual_final_fp"),
        "paired_resolved_count": sum(1 for row in sample_rows if row["paired_normal_resolvable"]),
        "paired_persistent_count": sum(1 for row in sample_rows if not row["paired_normal_resolvable"]),
        "to_baseline_tp": to_baseline_metrics["tp"],
        "to_baseline_fp": to_baseline_metrics["fp"],
        "to_baseline_fn": config["truth_total"] - int(to_baseline_metrics["tp"]),
        "to_baseline_precision": f"{to_float(to_baseline_metrics['precision']):.6f}",
        "to_baseline_recall": f"{to_float(to_baseline_metrics['recall']):.6f}",
        "to_baseline_f1": f"{to_float(to_baseline_metrics['f1']):.6f}",
        "to_final_tp": to_final_metrics["tp"],
        "to_final_fp": to_final_metrics["fp"],
        "to_final_fn": config["truth_total"] - int(to_final_metrics["tp"]),
        "to_final_precision": f"{to_float(to_final_metrics['precision']):.6f}",
        "to_final_recall": f"{to_float(to_final_metrics['recall']):.6f}",
        "to_final_f1": f"{to_float(to_final_metrics['f1']):.6f}",
        "paired_longphase_tp": paired_longphase_metrics["tp"],
        "paired_longphase_fp": paired_longphase_metrics["fp"],
        "paired_longphase_fn": config["truth_total"] - int(paired_longphase_metrics["tp"]),
        "paired_longphase_precision": f"{to_float(paired_longphase_metrics['precision']):.6f}",
        "paired_longphase_recall": f"{to_float(paired_longphase_metrics['recall']):.6f}",
        "paired_longphase_f1": f"{to_float(paired_longphase_metrics['f1']):.6f}",
        "paired_final_tp": paired_final_metrics["tp"],
        "paired_final_fp": paired_final_metrics["fp"],
        "paired_final_fn": config["truth_total"] - int(paired_final_metrics["tp"]),
        "paired_final_precision": f"{to_float(paired_final_metrics['precision']):.6f}",
        "paired_final_recall": f"{to_float(paired_final_metrics['recall']):.6f}",
        "paired_final_f1": f"{to_float(paired_final_metrics['f1']):.6f}",
        "to_fp_rule_trigger_count": len(to_post_removed_fp_keys),
        "to_tp_rule_trigger_count": len(to_post_removed_tp_keys),
        "paired_fp_rule_trigger_count": len(paired_post_removed_fp_keys),
        "paired_tp_rule_trigger_count": len(paired_post_removed_tp_keys),
        "to_round_dir": str(config["to_round_dir"]),
        "paired_dir": str(config["paired_dir"]),
        "to_raw_vcf": str(config["to_raw_vcf"]),
        "paired_raw_vcf": str(config["paired_raw_vcf"]),
    }

    return {
        "summary_row": summary_row,
        "class_summary_rows": class_summary_rows,
        "paired_summary_rows": paired_summary_rows,
        "crosstab_rows": crosstab_rows,
        "feature_group_rows": feature_group_rows,
        "sweep_rows": sweep_rows,
        "final_tp_rows": final_tp_rows,
        "final_fp_rows": final_fp_rows,
    }


def build_workspace_report(output_dir: Path, sample_summaries: list[dict[str, object]], best_rules: list[dict[str, object]]) -> None:
    lines = [
        "# TO FP provenance analysis workspace",
        "",
        "## Sample summary",
        "",
        "| sample | to_raw_fp_count | caller_pon_filtered | caller_nonpon_filtered | longphase_to_removed | to_postprocess_removed | to_residual_final_fp | paired_resolved_count | paired_persistent_count | to_final_f1 | paired_final_f1 |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in sample_summaries:
        lines.append(
            f"| {row['sample']} | {row['to_raw_fp_count']} | {row['caller_pon_filtered']} | {row['caller_nonpon_filtered']} | {row['longphase_to_removed']} | {row['to_postprocess_removed']} | {row['to_residual_final_fp']} | {row['paired_resolved_count']} | {row['paired_persistent_count']} | {row['to_final_f1']} | {row['paired_final_f1']} |"
        )
    lines.extend(
        [
            "",
            "## Discovery top rules",
            "",
            "| sample | rule_id | rule_label | fp_removed | tp_removed | f1_after | delta_f1_vs_final |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in best_rules:
        lines.append(
            f"| {row['sample']} | {row['rule_id']} | {row['rule_label']} | {row['fp_removed']} | {row['tp_removed']} | {row['f1_after']} | {row['delta_f1_vs_final']} |"
        )
    (output_dir / "analysis_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    samples = args.sample or DEFAULT_SAMPLES
    for sample in samples:
        if sample not in SAMPLE_CONFIGS:
            raise SystemExit(f"Unsupported sample: {sample}")

    output_dir = ensure_dir(Path(args.output_dir).resolve())
    sample_manifest_rows = []
    sample_summary_rows = []
    class_summary_rows = []
    paired_summary_rows = []
    crosstab_rows = []
    feature_group_rows = []
    rule_sweep_rows = []

    per_sample_results = {}
    for sample in samples:
        config = load_configs(sample)
        config["verification_rule_source"] = args.verification_rule_source
        sample_manifest_rows.append(
            {
                "sample": sample,
                "platform": config["platform"],
                "verification_rule_source": args.verification_rule_source,
                "truth_total": config["truth_total"],
                "to_round_dir": str(config["to_round_dir"]),
                "paired_dir": str(config["paired_dir"]),
                "to_raw_vcf": str(config["to_raw_vcf"]),
                "paired_raw_vcf": str(config["paired_raw_vcf"]),
                "truth_vcf": str(config["truth_vcf"]),
                "truth_bed": str(config["truth_bed"]),
            }
        )
        result = analyze_sample(config, output_dir)
        per_sample_results[sample] = result
        sample_summary_rows.append(result["summary_row"])
        class_summary_rows.extend(result["class_summary_rows"])
        paired_summary_rows.extend(result["paired_summary_rows"])
        crosstab_rows.extend(result["crosstab_rows"])
        feature_group_rows.extend(result["feature_group_rows"])
        rule_sweep_rows.extend(result["sweep_rows"])

    # discovery -> validation: take top 10 HCC1395 rules and apply to DORADO current final set
    hcc1395_top = [row for row in per_sample_results["HCC1395"]["sweep_rows"] if float(row["delta_f1_vs_final"]) > -1.0][:10]
    validation_rules = []
    lookup = {
        row["rule_id"]: row
        for row in build_candidate_rules("HCC1395", args.verification_rule_source)
    }
    dorado_result = per_sample_results.get("HCC1395_DORADO")
    if dorado_result:
        selected_rules = [lookup[row["rule_id"]] for row in hcc1395_top if row["rule_id"] in lookup]
        validation_rules = run_rule_sweep(
            "HCC1395_DORADO",
            dorado_result["final_tp_rows"],
            dorado_result["final_fp_rows"],
            load_configs("HCC1395_DORADO")["metrics_to"],
            selected_rules,
            phase="validation",
            verification_rule_source=args.verification_rule_source,
        )
        for rank, row in enumerate(validation_rules, start=1):
            row["rule_rank"] = rank
        rule_sweep_rows.extend(validation_rules)

    write_tsv_rows(output_dir / "sample_manifest.tsv", [
        "sample", "platform", "verification_rule_source", "truth_total", "to_round_dir", "paired_dir", "to_raw_vcf", "paired_raw_vcf", "truth_vcf", "truth_bed"
    ], sample_manifest_rows)
    write_tsv_rows(output_dir / "sample_level_summary.tsv", SAMPLE_SUMMARY_FIELDS, sample_summary_rows)
    write_tsv_rows(output_dir / "fp_primary_class_summary.tsv", CLASS_SUMMARY_FIELDS, class_summary_rows)
    write_tsv_rows(output_dir / "paired_oracle_resolution_summary.tsv", PAIRED_REASON_FIELDS, paired_summary_rows)
    write_tsv_rows(output_dir / "class_by_paired_reason.tsv", CLASS_PAIRED_CROSSTAB_FIELDS, crosstab_rows)
    write_tsv_rows(output_dir / "tp_fp_feature_comparison.tsv", FEATURE_GROUP_FIELDS, feature_group_rows)
    write_tsv_rows(output_dir / "to_rule_sweep_summary.tsv", RULE_SWEEP_FIELDS, rule_sweep_rows)

    build_workspace_report(output_dir, sample_summary_rows, [row for row in rule_sweep_rows if row["phase"] == "discovery"][:10])
    print(f"[build_to_fp_provenance_analysis] Wrote workspace to {output_dir}")


if __name__ == "__main__":
    main()
