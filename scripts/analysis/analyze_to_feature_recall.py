#!/usr/bin/env python3
"""Audit tumor-only label upgrades/conflicts against TP/FP and caller+methylation features."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import Counter
from pathlib import Path
from typing import Callable, Dict, Iterable, List

from research_common import compute_metrics, read_json, to_float, to_int, write_tsv_rows


PER_REGION_FIELDS = [
    "sample",
    "platform",
    "source_scope",
    "region_key",
    "agreement_type",
    "class_shift",
    "cluster_class",
    "label_class",
    "verification_class",
    "reads",
    "num_cpgs",
    "passed_gating",
    "vaf",
    "qual",
    "gq",
    "dp",
    "ad_ref",
    "ad_alt",
    "allele_delta",
    "cramers_v",
    "pairwise_median_dist",
    "pairwise_mean_dist",
    "global_p",
    "quality_score",
    "low_vaf_high_ad",
    "low_vaf_high_ad_cv",
    "legacy_combo",
    "focus_tags",
]

FOCUS_MATRIX_FIELDS = [
    "sample",
    "platform",
    "focus_group",
    "tp_count",
    "fp_count",
    "total_count",
    "tp_fraction_within_group",
    "fp_fraction_within_group",
    "tp_scope_fraction",
    "fp_scope_fraction",
    "median_vaf",
    "median_allele_delta",
    "median_cramers_v",
    "median_pairwise_median_dist",
    "median_gq",
    "median_qual",
    "top_class_shift",
    "top_agreement_type",
]

RULE_FIELDS = [
    "sample",
    "platform",
    "rule_id",
    "trigger_count",
    "tp_removed",
    "fp_removed",
    "precision",
    "recall",
    "f1",
    "delta_f1_vs_intersubmod",
    "baseline_f1",
    "notes",
]

TOP_FIELDS = [
    "sample",
    "platform",
    "source_scope",
    "region_key",
    "agreement_type",
    "class_shift",
    "cluster_class",
    "label_class",
    "vaf",
    "qual",
    "gq",
    "dp",
    "ad_ref",
    "ad_alt",
    "allele_delta",
    "cramers_v",
    "pairwise_median_dist",
    "quality_score",
    "focus_tags",
]

SHIFT_FIELDS = [
    "sample",
    "platform",
    "agreement_type",
    "class_shift",
    "tp_count",
    "fp_count",
    "total_count",
    "tp_fraction_within_shift",
    "fp_fraction_within_shift",
    "tp_scope_fraction",
    "fp_scope_fraction",
    "median_vaf",
    "median_allele_delta",
    "median_cramers_v",
    "median_pairwise_median_dist",
    "median_gq",
    "median_qual",
]


FocusFunc = Callable[[Dict[str, object]], bool]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", required=True, help="Tumor-only run directory")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--vaf-threshold", type=float, default=0.15, help="Low-VAF threshold")
    parser.add_argument("--allele-delta-threshold", type=float, default=0.15, help="High AlleleDelta threshold")
    parser.add_argument("--cv-threshold", type=float, default=0.05, help="Low CramersV threshold")
    parser.add_argument("--legacy-allele-delta-threshold", type=float, default=0.25, help="Legacy rule AlleleDelta threshold")
    parser.add_argument("--legacy-cv-threshold", type=float, default=0.05, help="Legacy rule CramersV threshold")
    parser.add_argument("--legacy-vaf-threshold", type=float, default=0.24, help="Legacy rule VAF threshold")
    parser.add_argument("--legacy-qual-threshold", type=float, default=0.75, help="Legacy rule QUAL threshold")
    parser.add_argument("--top-n", type=int, default=100, help="Top suspicious TP/FP rows to export")
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def load_summary_rows(path: Path, scope: str) -> Dict[str, Dict[str, object]]:
    rows: Dict[str, Dict[str, object]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            region_key = f"{row['Chr']}:{row['Pos']}:{row['Ref']}:{row['Alt']}"
            rows[region_key] = {
                "source_scope": scope,
                "verification_class": row.get("VerificationClass", ""),
                "reads": to_int(row.get("NumReads")),
                "num_cpgs": to_int(row.get("NumCpGs")),
                "passed_gating": row.get("PassedGating", ""),
                "allele_delta": to_float(row.get("AlleleDelta")),
                "cramers_v": to_float(row.get("CramersV")),
                "pairwise_median_dist": to_float(row.get("PairwiseMedianDist")),
                "pairwise_mean_dist": to_float(row.get("PairwiseMeanDist")),
                "global_p": to_float(row.get("GlobalP")),
                "quality_score": to_float(row.get("Quality_Score")),
            }
    return rows


def load_vcf_features(path: Path, scope: str) -> Dict[str, Dict[str, object]]:
    rows: Dict[str, Dict[str, object]] = {}
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            chrom, pos, _vid, ref, alt, qual, _filt, _info, fmt, sample = fields[:10]
            fmt_map = dict(zip(fmt.split(":"), sample.split(":")))
            ad_ref = math.nan
            ad_alt = math.nan
            if "AD" in fmt_map:
                parts = fmt_map["AD"].split(",")
                if parts:
                    ad_ref = to_float(parts[0])
                if len(parts) > 1:
                    ad_alt = to_float(parts[1])
            vaf = math.nan
            for key in ("AF", "VAF"):
                if key in fmt_map:
                    vaf = to_float(fmt_map[key].split(",")[0])
                    break
            rows[f"{chrom}:{pos}:{ref}:{alt}"] = {
                "source_scope": scope,
                "qual": to_float(None if qual == "." else qual),
                "gq": to_int(fmt_map.get("GQ")),
                "dp": to_int(fmt_map.get("DP")),
                "vaf": vaf,
                "ad_ref": ad_ref,
                "ad_alt": ad_alt,
            }
    return rows


def median(values: Iterable[float]) -> str:
    clean = sorted(value for value in values if not math.isnan(value))
    if not clean:
        return ""
    size = len(clean)
    mid = size // 2
    if size % 2:
        return f"{clean[mid]:.6f}"
    return f"{((clean[mid - 1] + clean[mid]) / 2.0):.6f}"


def build_focus_groups() -> Dict[str, FocusFunc]:
    return {
        "agreement_label_upgrade": lambda row: row["agreement_type"] == "label_upgrade",
        "agreement_conflict": lambda row: row["agreement_type"] == "conflict",
        "shift_Weak_to_Strong": lambda row: row["class_shift"] == "Weak->Strong",
        "shift_Noise_to_Strong": lambda row: row["class_shift"] == "Noise->Strong",
        "shift_Weak_to_Subclone": lambda row: row["class_shift"] == "Weak->Subclone",
        "shift_Noise_to_Subclone": lambda row: row["class_shift"] == "Noise->Subclone",
        "rule_low_vaf_high_ad": lambda row: row["low_vaf_high_ad"],
        "rule_low_vaf_high_ad_cv": lambda row: row["low_vaf_high_ad_cv"],
        "rule_legacy_combo": lambda row: row["legacy_combo"],
        "combo_label_upgrade_low_vaf_high_ad": lambda row: row["agreement_type"] == "label_upgrade"
        and row["low_vaf_high_ad"],
        "combo_conflict_low_vaf_high_ad": lambda row: row["agreement_type"] == "conflict"
        and row["low_vaf_high_ad"],
        "combo_Weak_to_Strong_low_vaf_high_ad": lambda row: row["class_shift"] == "Weak->Strong"
        and row["low_vaf_high_ad"],
        "combo_Noise_to_Strong_low_vaf_high_ad": lambda row: row["class_shift"] == "Noise->Strong"
        and row["low_vaf_high_ad"],
        "combo_Noise_to_Strong_low_vaf_high_ad_cv": lambda row: row["class_shift"] == "Noise->Strong"
        and row["low_vaf_high_ad_cv"],
        "combo_conflict_Noise_to_Strong_low_vaf_high_ad": lambda row: row["agreement_type"] == "conflict"
        and row["class_shift"] == "Noise->Strong"
        and row["low_vaf_high_ad"],
    }


def row_to_output(row: Dict[str, object]) -> Dict[str, object]:
    output: Dict[str, object] = {}
    for field in TOP_FIELDS:
        value = row.get(field, "")
        if isinstance(value, float):
            output[field] = "" if math.isnan(value) else f"{value:.6f}"
        else:
            output[field] = value
    return output


def main() -> None:
    args = parse_args()
    sample_dir = Path(args.sample_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    context = read_json(sample_dir / "round_context.json")
    metrics = read_json(sample_dir / "metrics.json")
    sample = str(context.get("sample") or sample_dir.name)
    platform = str(context.get("platform") or "")
    truth_total = int(metrics.get("truth_total", 0))
    filtered = metrics.get("filtered", {})
    baseline_tp = int(filtered.get("tp", 0))
    baseline_fp = int(filtered.get("fp", 0))
    baseline_f1 = float(filtered.get("f1", 0.0))

    tp_summary = load_summary_rows(sample_dir / "step05_intersubmod" / "intersubmod_tp" / "significance_summary.csv", "tp")
    fp_summary = load_summary_rows(sample_dir / "step05_intersubmod" / "intersubmod_fp" / "significance_summary.csv", "fp")
    summary_index = {**tp_summary, **fp_summary}
    vcf_index = {
        **load_vcf_features(sample_dir / "step04_benchmark_longphase_to" / "filtered_snv_tp.vcf.gz", "tp"),
        **load_vcf_features(sample_dir / "step04_benchmark_longphase_to" / "filtered_snv_fp.vcf.gz", "fp"),
    }

    per_region_rows: List[Dict[str, object]] = []
    enriched_rows: List[Dict[str, object]] = []
    with (sample_dir / "label_cluster_agreement.tsv").open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            region_key = str(row.get("region_key", ""))
            summary = summary_index.get(region_key, {})
            vcf = vcf_index.get(region_key, {})
            if not summary:
                continue

            vaf = to_float(vcf.get("vaf"))
            qual = to_float(vcf.get("qual"))
            allele_delta = to_float(summary.get("allele_delta"))
            cramers_v = to_float(summary.get("cramers_v"))

            low_vaf_high_ad = (
                not math.isnan(vaf)
                and vaf < args.vaf_threshold
                and not math.isnan(allele_delta)
                and allele_delta > args.allele_delta_threshold
            )
            low_vaf_high_ad_cv = low_vaf_high_ad and not math.isnan(cramers_v) and cramers_v < args.cv_threshold
            legacy_combo = (
                (not math.isnan(qual) and qual < args.legacy_qual_threshold)
                or (
                    not math.isnan(allele_delta)
                    and allele_delta > args.legacy_allele_delta_threshold
                    and not math.isnan(cramers_v)
                    and cramers_v < args.legacy_cv_threshold
                    and not math.isnan(vaf)
                    and vaf < args.legacy_vaf_threshold
                )
            )

            focus_tags: List[str] = []
            if row.get("agreement_type") == "label_upgrade":
                focus_tags.append("label_upgrade")
            if row.get("agreement_type") == "conflict":
                focus_tags.append("conflict")
            if row.get("class_shift") == "Weak->Strong":
                focus_tags.append("weak_to_strong")
            if row.get("class_shift") == "Noise->Strong":
                focus_tags.append("noise_to_strong")
            if low_vaf_high_ad:
                focus_tags.append("low_vaf_high_ad")
            if low_vaf_high_ad_cv:
                focus_tags.append("low_vaf_high_ad_cv")
            if legacy_combo:
                focus_tags.append("legacy_combo")

            enriched = {
                "sample": sample,
                "platform": platform,
                "source_scope": summary.get("source_scope", row.get("source_scope", "")),
                "region_key": region_key,
                "agreement_type": row.get("agreement_type", ""),
                "class_shift": row.get("class_shift", ""),
                "cluster_class": row.get("cluster_class", ""),
                "label_class": row.get("label_class", ""),
                "verification_class": summary.get("verification_class", ""),
                "reads": summary.get("reads", 0),
                "num_cpgs": summary.get("num_cpgs", 0),
                "passed_gating": summary.get("passed_gating", ""),
                "vaf": vaf,
                "qual": qual,
                "gq": to_int(vcf.get("gq")),
                "dp": to_int(vcf.get("dp")),
                "ad_ref": to_float(vcf.get("ad_ref")),
                "ad_alt": to_float(vcf.get("ad_alt")),
                "allele_delta": allele_delta,
                "cramers_v": cramers_v,
                "pairwise_median_dist": to_float(summary.get("pairwise_median_dist")),
                "pairwise_mean_dist": to_float(summary.get("pairwise_mean_dist")),
                "global_p": to_float(summary.get("global_p")),
                "quality_score": to_float(summary.get("quality_score")),
                "low_vaf_high_ad": low_vaf_high_ad,
                "low_vaf_high_ad_cv": low_vaf_high_ad_cv,
                "legacy_combo": legacy_combo,
                "focus_tags": ",".join(focus_tags),
            }
            enriched_rows.append(enriched)

            output_row: Dict[str, object] = {}
            for field in PER_REGION_FIELDS:
                value = enriched.get(field, "")
                if isinstance(value, float):
                    output_row[field] = "" if math.isnan(value) else f"{value:.6f}"
                else:
                    output_row[field] = value
            per_region_rows.append(output_row)

    focus_groups = build_focus_groups()
    tp_total = sum(1 for row in enriched_rows if row["source_scope"] == "tp")
    fp_total = sum(1 for row in enriched_rows if row["source_scope"] == "fp")

    matrix_rows: List[Dict[str, object]] = []
    rule_rows: List[Dict[str, object]] = []
    shift_rows: List[Dict[str, object]] = []
    for focus_group, matcher in focus_groups.items():
        subset = [row for row in enriched_rows if matcher(row)]
        tp_subset = [row for row in subset if row["source_scope"] == "tp"]
        fp_subset = [row for row in subset if row["source_scope"] == "fp"]
        class_shifts = Counter(str(row["class_shift"]) for row in subset if row["class_shift"])
        agreement_types = Counter(str(row["agreement_type"]) for row in subset if row["agreement_type"])
        matrix_rows.append(
            {
                "sample": sample,
                "platform": platform,
                "focus_group": focus_group,
                "tp_count": len(tp_subset),
                "fp_count": len(fp_subset),
                "total_count": len(subset),
                "tp_fraction_within_group": f"{(len(tp_subset) / len(subset)):.6f}" if subset else "",
                "fp_fraction_within_group": f"{(len(fp_subset) / len(subset)):.6f}" if subset else "",
                "tp_scope_fraction": f"{(len(tp_subset) / tp_total):.6f}" if tp_total else "",
                "fp_scope_fraction": f"{(len(fp_subset) / fp_total):.6f}" if fp_total else "",
                "median_vaf": median(row["vaf"] for row in subset),
                "median_allele_delta": median(row["allele_delta"] for row in subset),
                "median_cramers_v": median(row["cramers_v"] for row in subset),
                "median_pairwise_median_dist": median(row["pairwise_median_dist"] for row in subset),
                "median_gq": median(float(row["gq"]) for row in subset),
                "median_qual": median(row["qual"] for row in subset),
                "top_class_shift": "; ".join(f"{key}:{value}" for key, value in class_shifts.most_common(3)),
                "top_agreement_type": "; ".join(f"{key}:{value}" for key, value in agreement_types.most_common(3)),
            }
        )

        metrics_after = compute_metrics(baseline_tp - len(tp_subset), baseline_fp - len(fp_subset), truth_total)
        notes: List[str] = []
        if "low_vaf_high_ad" in focus_group:
            notes.append("contains low_vaf_high_ad")
        if "conflict" in focus_group:
            notes.append("conflict subset")
        if "label_upgrade" in focus_group:
            notes.append("label_upgrade subset")
        if "Noise_to_Strong" in focus_group:
            notes.append("targets Noise->Strong")
        if "Weak_to_Strong" in focus_group:
            notes.append("targets Weak->Strong")
        if "legacy" in focus_group:
            notes.append("legacy rule")
        rule_rows.append(
            {
                "sample": sample,
                "platform": platform,
                "rule_id": focus_group,
                "trigger_count": len(subset),
                "tp_removed": len(tp_subset),
                "fp_removed": len(fp_subset),
                "precision": f"{metrics_after['precision']:.6f}",
                "recall": f"{metrics_after['recall']:.6f}",
                "f1": f"{metrics_after['f1']:.6f}",
                "delta_f1_vs_intersubmod": f"{(metrics_after['f1'] - baseline_f1):.6f}",
                "baseline_f1": f"{baseline_f1:.6f}",
                "notes": "; ".join(notes),
            }
        )

    grouped: Dict[tuple[str, str], List[Dict[str, object]]] = {}
    for row in enriched_rows:
        key = (str(row["agreement_type"]), str(row["class_shift"]))
        grouped.setdefault(key, []).append(row)
    for (agreement_type, class_shift), subset in sorted(grouped.items()):
        if agreement_type not in {"label_upgrade", "conflict"} and class_shift in {"", "unchanged"}:
            continue
        tp_subset = [row for row in subset if row["source_scope"] == "tp"]
        fp_subset = [row for row in subset if row["source_scope"] == "fp"]
        total_count = len(subset)
        shift_rows.append(
            {
                "sample": sample,
                "platform": platform,
                "agreement_type": agreement_type,
                "class_shift": class_shift,
                "tp_count": len(tp_subset),
                "fp_count": len(fp_subset),
                "total_count": total_count,
                "tp_fraction_within_shift": f"{(len(tp_subset) / total_count):.6f}" if total_count else "",
                "fp_fraction_within_shift": f"{(len(fp_subset) / total_count):.6f}" if total_count else "",
                "tp_scope_fraction": f"{(len(tp_subset) / tp_total):.6f}" if tp_total else "",
                "fp_scope_fraction": f"{(len(fp_subset) / fp_total):.6f}" if fp_total else "",
                "median_vaf": median(row["vaf"] for row in subset),
                "median_allele_delta": median(row["allele_delta"] for row in subset),
                "median_cramers_v": median(row["cramers_v"] for row in subset),
                "median_pairwise_median_dist": median(row["pairwise_median_dist"] for row in subset),
                "median_gq": median(float(row["gq"]) for row in subset),
                "median_qual": median(row["qual"] for row in subset),
            }
        )

    top_fp = [
        row
        for row in enriched_rows
        if row["source_scope"] == "fp"
        and (
            row["agreement_type"] in {"label_upgrade", "conflict"}
            or row["class_shift"] in {"Weak->Strong", "Noise->Strong"}
            or row["low_vaf_high_ad"]
        )
    ]
    top_tp = [
        row
        for row in enriched_rows
        if row["source_scope"] == "tp"
        and (
            row["agreement_type"] in {"label_upgrade", "conflict"}
            or row["class_shift"] in {"Weak->Strong", "Noise->Strong"}
            or row["low_vaf_high_ad"]
        )
    ]

    def top_sort_key(row: Dict[str, object]):
        shift_priority = {
            "Noise->Strong": 3,
            "Weak->Strong": 2,
            "Noise->Subclone": 1,
            "Weak->Subclone": 1,
        }.get(str(row["class_shift"]), 0)
        agreement_priority = {"conflict": 2, "label_upgrade": 1}.get(str(row["agreement_type"]), 0)
        return (
            int(bool(row["low_vaf_high_ad_cv"])),
            int(bool(row["low_vaf_high_ad"])),
            shift_priority,
            agreement_priority,
            -(row["allele_delta"] if not math.isnan(row["allele_delta"]) else -1.0),
            row["vaf"] if not math.isnan(row["vaf"]) else 1.0,
            -(row["pairwise_median_dist"] if not math.isnan(row["pairwise_median_dist"]) else -1.0),
        )

    top_fp = sorted(top_fp, key=top_sort_key, reverse=True)[: args.top_n]
    top_tp = sorted(top_tp, key=top_sort_key, reverse=True)[: args.top_n]

    write_tsv_rows(output_dir / "to_feature_per_region.tsv", PER_REGION_FIELDS, per_region_rows)
    write_tsv_rows(output_dir / "to_feature_focus_matrix.tsv", FOCUS_MATRIX_FIELDS, matrix_rows)
    write_tsv_rows(output_dir / "to_feature_shift_breakdown.tsv", SHIFT_FIELDS, shift_rows)
    write_tsv_rows(
        output_dir / "to_feature_rule_comparison.tsv",
        RULE_FIELDS,
        sorted(rule_rows, key=lambda row: float(row["delta_f1_vs_intersubmod"]), reverse=True),
    )
    write_tsv_rows(output_dir / "to_feature_top_fp_candidates.tsv", TOP_FIELDS, [row_to_output(row) for row in top_fp])
    write_tsv_rows(output_dir / "to_feature_top_tp_controls.tsv", TOP_FIELDS, [row_to_output(row) for row in top_tp])

    selected = {row["focus_group"]: row for row in matrix_rows}
    md_lines = [
        f"# TO Feature Recall Summary: {sample}",
        "",
        f"- sample dir: `{sample_dir}`",
        f"- baseline InterSubMod F1: `{baseline_f1:.6f}`",
        f"- analyzed TP rows: `{tp_total}`",
        f"- analyzed FP rows: `{fp_total}`",
        f"- thresholds: `VAF<{args.vaf_threshold}`, `AlleleDelta>{args.allele_delta_threshold}`, `CramersV<{args.cv_threshold}`",
        "",
        "## Key focus groups",
        "",
        "| focus_group | tp_count | fp_count | fp_fraction_within_group | tp_scope_fraction | fp_scope_fraction | median_vaf | median_allele_delta | top_class_shift |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for key in (
        "agreement_label_upgrade",
        "agreement_conflict",
        "shift_Weak_to_Strong",
        "shift_Noise_to_Strong",
        "rule_low_vaf_high_ad",
        "combo_label_upgrade_low_vaf_high_ad",
        "combo_conflict_low_vaf_high_ad",
        "combo_Weak_to_Strong_low_vaf_high_ad",
        "combo_Noise_to_Strong_low_vaf_high_ad",
    ):
        row = selected.get(key)
        if not row:
            continue
        md_lines.append(
            f"| {key} | {row['tp_count']} | {row['fp_count']} | {row['fp_fraction_within_group']} | {row['tp_scope_fraction']} | {row['fp_scope_fraction']} | {row['median_vaf']} | {row['median_allele_delta']} | {row['top_class_shift']} |"
        )
    md_lines.extend(
        [
        "",
        "## Best removal rules",
        "",
            "| rule_id | trigger_count | tp_removed | fp_removed | f1 | delta_f1_vs_intersubmod | notes |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    sorted_rules = sorted(rule_rows, key=lambda row: float(row["delta_f1_vs_intersubmod"]), reverse=True)
    for row in sorted_rules[:10]:
        md_lines.append(
            f"| {row['rule_id']} | {row['trigger_count']} | {row['tp_removed']} | {row['fp_removed']} | {row['f1']} | {row['delta_f1_vs_intersubmod']} | {row['notes']} |"
        )
    md_lines.extend(
        [
            "",
            "## Label Upgrade / Conflict class shifts",
            "",
            "| agreement_type | class_shift | tp_count | fp_count | fp_fraction_within_shift | median_vaf | median_allele_delta |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in sorted(shift_rows, key=lambda item: (item["agreement_type"], -item["total_count"], item["class_shift"]))[:12]:
        md_lines.append(
            f"| {row['agreement_type']} | {row['class_shift']} | {row['tp_count']} | {row['fp_count']} | {row['fp_fraction_within_shift']} | {row['median_vaf']} | {row['median_allele_delta']} |"
        )
    (output_dir / "to_feature_summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"[analyze_to_feature_recall] Wrote {output_dir / 'to_feature_per_region.tsv'}")
    print(f"[analyze_to_feature_recall] Wrote {output_dir / 'to_feature_focus_matrix.tsv'}")
    print(f"[analyze_to_feature_recall] Wrote {output_dir / 'to_feature_shift_breakdown.tsv'}")
    print(f"[analyze_to_feature_recall] Wrote {output_dir / 'to_feature_rule_comparison.tsv'}")
    print(f"[analyze_to_feature_recall] Wrote {output_dir / 'to_feature_top_fp_candidates.tsv'}")


if __name__ == "__main__":
    main()
