#!/usr/bin/env python3
"""Phase 2 analysis: paired raw pileup/full benchmarks plus finer feature interval/orthogonality validation."""

from __future__ import annotations

import argparse
import gzip
import itertools
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import pandas as pd

from research_common import compute_metrics, ensure_dir, markdown_table, parse_variant_counts, write_tsv_rows


REPO_ROOT = Path("/big8_disk/liaoyoyo2001/InterSubMod")
OUTPUT_ROOT_DEFAULT = (
    Path("/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds")
    / "20260311_phase2_paired_model_feature_analysis"
)
TRUTH_VCF = Path("/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz")
TRUTH_BED = Path("/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed")
BENCHMARK_SPLIT_SCRIPT = REPO_ROOT / "scripts/pipeline/utils/benchmark_split_snv_vcf.sh"


@dataclass
class PairedRawDataset:
    dataset_id: str
    label: str
    sample: str
    platform: str
    pileup_vcf: Path
    full_vcf: Path
    merged_vcf: Path
    longphase_f1: float
    intersubmod_f1: float
    benchmark_comparison_tsv: Path


@dataclass
class RescueDataset:
    dataset_id: str
    label: str
    sample: str
    platform: str
    mode: str
    caller: str
    joined_tsv: Path
    baseline_tp: int
    baseline_fp: int
    truth_total: int
    source_scope_note: str


PAIRED_RAW_DATASETS: List[PairedRawDataset] = [
    PairedRawDataset(
        dataset_id="hcc1395_5khz_paired",
        label="HCC1395 5kHz paired",
        sample="HCC1395",
        platform="ONT_5kHz",
        pileup_vcf=Path("/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/tmp/vcf_output/pileup_filter.vcf"),
        full_vcf=Path("/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/tmp/vcf_output/full_alignment_filter.vcf"),
        merged_vcf=Path("/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/output.vcf.gz"),
        longphase_f1=0.852208,
        intersubmod_f1=0.853200,
        benchmark_comparison_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
            "20260307_paired_pure_full_columns_official/HCC1395/benchmark_comparison.tsv"
        ),
    ),
    PairedRawDataset(
        dataset_id="hcc1395_dorado_paired",
        label="HCC1395 DORADO paired",
        sample="HCC1395_DORADO",
        platform="ONT_Dorado",
        pileup_vcf=Path("/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/tmp/vcf_output/pileup_filter.vcf"),
        full_vcf=Path("/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/tmp/vcf_output/full_alignment_filter.vcf"),
        merged_vcf=Path("/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/output.vcf.gz"),
        longphase_f1=0.859176,
        intersubmod_f1=0.859000,
        benchmark_comparison_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
            "20260307_paired_pure_full_columns_official/HCC1395_DORADO/benchmark_comparison.tsv"
        ),
    ),
]


RESCUE_DATASETS: List[RescueDataset] = [
    RescueDataset(
        dataset_id="hcc1395_5khz_to",
        label="HCC1395 5kHz TO",
        sample="HCC1395",
        platform="ONT_5kHz",
        mode="to-pure",
        caller="ClairS-TO",
        joined_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
            "20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_joined_features.tsv"
        ),
        baseline_tp=28396,
        baseline_fp=11843,
        truth_total=39447,
        source_scope_note="full tagged BAM + candidate-specific rescue join",
    ),
    RescueDataset(
        dataset_id="hcc1395_5khz_paired",
        label="HCC1395 5kHz paired",
        sample="HCC1395",
        platform="ONT_5kHz",
        mode="paired-pure",
        caller="ClairS",
        joined_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
            "20260308_borderline_rescue_analysis/HCC1395_paired/eval/rescue_joined_features.tsv"
        ),
        baseline_tp=29754,
        baseline_fp=627,
        truth_total=39447,
        source_scope_note="candidate-specific rescue join with lower methylation coverage ceiling",
    ),
    RescueDataset(
        dataset_id="hcc1395_dorado_paired",
        label="HCC1395 DORADO paired",
        sample="HCC1395_DORADO",
        platform="ONT_Dorado",
        mode="paired-pure",
        caller="ClairS",
        joined_tsv=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260309_hcc1395_dorado_paired_candidate_rescue/eval/rescue_joined_features.tsv"
        ),
        baseline_tp=29889,
        baseline_fp=240,
        truth_total=39447,
        source_scope_note="candidate-specific rescue join with lower methylation coverage ceiling",
    ),
    RescueDataset(
        dataset_id="hcc1395_dorado_to",
        label="HCC1395 DORADO TO",
        sample="HCC1395_DORADO",
        platform="ONT_Dorado",
        mode="to",
        caller="ClairS-TO",
        joined_tsv=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260311_hcc1395_dorado_to_candidate_rescue/eval/rescue_joined_features.tsv"
        ),
        baseline_tp=28861,
        baseline_fp=11576,
        truth_total=39447,
        source_scope_note="candidate-window subset tagged BAM + candidate-specific rescue join",
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(OUTPUT_ROOT_DEFAULT), help="Output directory")
    parser.add_argument("--skip-benchmark", action="store_true", help="Skip rerunning paired raw benchmarks")
    return parser.parse_args()


def run_cmd(cmd: Sequence[str], cwd: Path | None = None) -> None:
    subprocess.run(list(cmd), check=True, cwd=str(cwd) if cwd else None)


def format_float(value: float) -> str:
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "inf"
    return f"{value:.6f}"


def load_joined(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
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
        else:
            df[column] = pd.NA

    for column in ["VerificationClass", "agreement_type", "cluster_class", "label_class", "DominantLabel", "class_shift", "source_scope", "downstream_status"]:
        if column in df.columns:
            df[column] = df[column].fillna("NA").replace("", "NA")
        else:
            df[column] = "NA"

    if "candidate_eligible" in df.columns:
        df["candidate_eligible"] = df["candidate_eligible"].fillna(False).map(
            lambda value: str(value).strip().lower() in {"1", "true", "yes", "y"}
        )
    else:
        df["candidate_eligible"] = False

    df["analysis_available"] = df["source_scope"].isin({"tp", "fp"})
    df["strong_subclone"] = df["VerificationClass"].isin({"Strong", "Subclone"})
    df["agreement_positive"] = df["agreement_type"].isin({"label_upgrade", "consistent_strong", "consistent_subclone"})
    return df


def evaluate_mask(df: pd.DataFrame, mask: pd.Series, baseline_tp: int, baseline_fp: int, truth_total: int) -> Dict[str, float | int]:
    selected = df[mask]
    tp_rescued = int((selected["downstream_status"] == "caller_lost_tp").sum())
    fp_reintroduced = int((selected["downstream_status"] == "caller_removed_fp").sum())
    metrics = compute_metrics(baseline_tp + tp_rescued, baseline_fp + fp_reintroduced, truth_total)
    baseline = compute_metrics(baseline_tp, baseline_fp, truth_total)
    return {
        "trigger_count": int(mask.sum()),
        "tp_rescued": tp_rescued,
        "fp_reintroduced": fp_reintroduced,
        "precision": metrics["precision"],
        "recall": metrics["recall"],
        "f1": metrics["f1"],
        "delta_f1_vs_baseline": metrics["f1"] - baseline["f1"],
    }


def candidate_eval_frame(joined_df: pd.DataFrame) -> pd.DataFrame:
    keep_status = {"caller_lost_tp", "caller_removed_fp"}
    return joined_df[joined_df["candidate_eligible"] & joined_df["downstream_status"].isin(keep_status)].copy()


def ensure_bcftools() -> str:
    bcftools = shutil.which("bcftools")
    if not bcftools:
        raise RuntimeError("bcftools not found in PATH")
    return bcftools


def open_text_auto(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def prepare_snv_input(input_vcf: Path, output_vcf: Path) -> Path:
    ensure_dir(output_vcf.parent)
    with open_text_auto(input_vcf) as src, output_vcf.open("w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith("##"):
                dst.write(line)
                continue
            if line.startswith("#CHROM"):
                dst.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            ref = parts[3]
            alt = parts[4]
            if len(ref) != 1:
                continue
            alts = alt.split(",")
            if not alts or any(len(a) != 1 for a in alts):
                continue
            dst.write("\t".join(parts[:8]) + "\n")
    return output_vcf


def benchmark_one_source(dataset: PairedRawDataset, source_name: str, source_vcf: Path, output_dir: Path, bcftools: str) -> Dict[str, object]:
    source_dir = output_dir / dataset.dataset_id / source_name
    ensure_dir(source_dir)
    snv_input = source_dir / "input_snv.vcf"
    prepare_snv_input(source_vcf, snv_input)
    run_cmd(
        [
            "bash",
            str(BENCHMARK_SPLIT_SCRIPT),
            "--input-vcf",
            str(snv_input),
            "--truth-vcf",
            str(TRUTH_VCF),
            "--truth-bed",
            str(TRUTH_BED),
            "--output-dir",
            str(source_dir / "benchmark"),
            "--no-plain",
        ]
    )
    counts = parse_variant_counts(source_dir / "benchmark" / "variant_counts.txt")
    metrics = compute_metrics(counts["TP_COUNT"], counts["FP_COUNT"], counts["TRUTH_TOTAL"])
    return {
        "dataset_id": dataset.dataset_id,
        "label": dataset.label,
        "sample": dataset.sample,
        "platform": dataset.platform,
        "source_name": source_name,
        "source_vcf": str(source_vcf),
        "snv_input_vcf": str(snv_input),
        "benchmark_dir": str(source_dir / "benchmark"),
        "pass_total": counts["PASS_TOTAL"],
        "truth_total": counts["TRUTH_TOTAL"],
        "tp": counts["TP_COUNT"],
        "fp": counts["FP_COUNT"],
        "fn": counts["FN_COUNT"],
        "precision": format_float(metrics["precision"]),
        "recall": format_float(metrics["recall"]),
        "f1": format_float(metrics["f1"]),
    }


def build_paired_raw_benchmark(output_dir: Path, skip_benchmark: bool) -> pd.DataFrame:
    bcftools = ensure_bcftools()
    rows: List[Dict[str, object]] = []
    sources = {
        "pileup_filter": lambda ds: ds.pileup_vcf,
        "full_alignment_filter": lambda ds: ds.full_vcf,
        "merged_output": lambda ds: ds.merged_vcf,
    }
    for dataset in PAIRED_RAW_DATASETS:
        for source_name, source_getter in sources.items():
            source_vcf = source_getter(dataset)
            benchmark_dir = output_dir / dataset.dataset_id / source_name / "benchmark"
            counts_path = benchmark_dir / "variant_counts.txt"
            if skip_benchmark and counts_path.exists():
                counts = parse_variant_counts(counts_path)
                metrics = compute_metrics(counts["TP_COUNT"], counts["FP_COUNT"], counts["TRUTH_TOTAL"])
                row = {
                    "dataset_id": dataset.dataset_id,
                    "label": dataset.label,
                    "sample": dataset.sample,
                    "platform": dataset.platform,
                    "source_name": source_name,
                    "source_vcf": str(source_vcf),
                    "snv_input_vcf": str((output_dir / dataset.dataset_id / source_name / "input_snv.vcf.gz")),
                    "benchmark_dir": str(benchmark_dir),
                    "pass_total": counts["PASS_TOTAL"],
                    "truth_total": counts["TRUTH_TOTAL"],
                    "tp": counts["TP_COUNT"],
                    "fp": counts["FP_COUNT"],
                    "fn": counts["FN_COUNT"],
                    "precision": format_float(metrics["precision"]),
                    "recall": format_float(metrics["recall"]),
                    "f1": format_float(metrics["f1"]),
                }
            else:
                row = benchmark_one_source(dataset, source_name, source_vcf, output_dir, bcftools)
            rows.append(row)

        # Append downstream reference methods for paired context.
        benchmark_df = pd.read_csv(dataset.benchmark_comparison_tsv, sep="\t")
        for method_name in ["ClairS", "LongPhase-S", "InterSubMod"]:
            part = benchmark_df[benchmark_df["method"] == method_name]
            if part.empty:
                continue
            part = part.iloc[0]
            rows.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "label": dataset.label,
                    "sample": dataset.sample,
                    "platform": dataset.platform,
                    "source_name": method_name,
                    "source_vcf": str(dataset.benchmark_comparison_tsv),
                    "snv_input_vcf": "",
                    "benchmark_dir": str(dataset.benchmark_comparison_tsv.parent),
                    "pass_total": int(part["calls_total"]),
                    "truth_total": int(part["truth_total"]),
                    "tp": int(part["tp"]),
                    "fp": int(part["fp"]),
                    "fn": int(part["fn"]),
                    "precision": format_float(float(part["precision"])),
                    "recall": format_float(float(part["recall"])),
                    "f1": format_float(float(part["f1"])),
                }
            )

    df = pd.DataFrame(rows)
    merged_f1 = (
        df[df["source_name"] == "merged_output"][["dataset_id", "f1", "tp", "fp", "fn"]]
        .rename(columns={"f1": "merged_f1", "tp": "merged_tp", "fp": "merged_fp", "fn": "merged_fn"})
    )
    df = df.merge(merged_f1, on="dataset_id", how="left")
    df["delta_f1_vs_merged_output"] = df.apply(
        lambda row: format_float(float(row["f1"]) - float(row["merged_f1"])) if row["source_name"] not in {"merged_output"} else "0.000000",
        axis=1,
    )
    df["delta_tp_vs_merged_output"] = df["tp"] - df["merged_tp"]
    df["delta_fp_vs_merged_output"] = df["fp"] - df["merged_fp"]
    df["delta_fn_vs_merged_output"] = df["fn"] - df["merged_fn"]
    return df


def fine_threshold_defs() -> List[Dict[str, object]]:
    def pct_tag(value: float) -> str:
        return f"{int(round(value * 100)):03d}"

    defs: List[Dict[str, object]] = []
    for threshold in [5, 8, 10, 12, 15, 18, 20, 25]:
        defs.append({"feature": "gq", "rule_id": f"gq_ge_{threshold}", "notes": f"GQ >= {threshold}", "func": lambda df, t=threshold: df["gq"] >= t})
    for threshold in [40, 50, 55, 60, 65, 70, 80]:
        defs.append({"feature": "Quality_Score", "rule_id": f"quality_ge_{threshold}", "notes": f"Quality_Score >= {threshold}", "func": lambda df, t=threshold: df["Quality_Score"] >= t})
    for threshold in [0.05, 0.10, 0.12, 0.15, 0.18, 0.20, 0.25, 0.30]:
        tag = pct_tag(threshold)
        defs.append({"feature": "PairwiseMedianDist", "rule_id": f"pairwise_ge_{tag}", "notes": f"PairwiseMedianDist >= {threshold:.2f}", "func": lambda df, t=threshold: df["PairwiseMedianDist"] >= t})
        defs.append({"feature": "PairwiseMedianDist", "rule_id": f"pairwise_le_{tag}", "notes": f"PairwiseMedianDist <= {threshold:.2f}", "func": lambda df, t=threshold: df["PairwiseMedianDist"] <= t})
    for threshold in [0.50, 0.70, 0.85, 0.90, 0.95, 0.99]:
        tag = pct_tag(threshold)
        defs.append({"feature": "hp_assign_rate", "rule_id": f"hp_assign_ge_{tag}", "notes": f"hp_assign_rate >= {threshold:.2f}", "func": lambda df, t=threshold: df["hp_assign_rate"] >= t})
    defs.extend(
        [
            {"feature": "agreement_positive", "rule_id": "agreement_positive", "notes": "agreement_type in label_upgrade/consistent_strong/consistent_subclone", "func": lambda df: df["agreement_positive"]},
            {"feature": "VerificationClass", "rule_id": "strong_subclone", "notes": "VerificationClass in Strong/Subclone", "func": lambda df: df["strong_subclone"]},
        ]
    )
    return defs


def feature_bins() -> Dict[str, Sequence[float]]:
    return {
        "gq": [0, 5, 8, 10, 12, 15, 18, 20, 25, math.inf],
        "Quality_Score": [0, 40, 50, 55, 60, 65, 70, 80, math.inf],
        "PairwiseMedianDist": [0.0, 0.05, 0.10, 0.12, 0.15, 0.18, 0.20, 0.25, 0.30, math.inf],
        "af": [0.0, 0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.50, 1.01],
        "AlleleDelta": [0.0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 1.01],
        "CramersV": [0.0, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 1.01],
        "hp_assign_rate": [0.0, 0.50, 0.70, 0.85, 0.90, 0.95, 0.99, 1.01],
    }


def bin_labels(edges: Sequence[float]) -> List[str]:
    labels: List[str] = []
    for left, right in zip(edges[:-1], edges[1:]):
        right_text = "inf" if math.isinf(right) else f"{right:g}"
        labels.append(f"[{left:g},{right_text})")
    return labels


def cut_series(series: pd.Series, edges: Sequence[float]) -> pd.Series:
    return pd.cut(series, bins=edges, right=False, labels=bin_labels(edges), include_lowest=True)


def artifact_low_vaf_high_adelta_low_cv(df: pd.DataFrame) -> pd.Series:
    return (df["af"] < 0.24) & (df["AlleleDelta"] > 0.25) & (df["CramersV"] < 0.05)


def orthogonal_masks(df: pd.DataFrame) -> Dict[str, pd.Series]:
    return {
        "gq_ge_10": df["gq"] >= 10,
        "gq_ge_15": df["gq"] >= 15,
        "gq_ge_20": df["gq"] >= 20,
        "quality_ge_60": df["Quality_Score"] >= 60,
        "quality_ge_70": df["Quality_Score"] >= 70,
        "pairwise_ge_015": df["PairwiseMedianDist"] >= 0.15,
        "pairwise_ge_020": df["PairwiseMedianDist"] >= 0.20,
        "pairwise_le_015": df["PairwiseMedianDist"] <= 0.15,
        "pairwise_le_020": df["PairwiseMedianDist"] <= 0.20,
        "hp_assign_ge_095": df["hp_assign_rate"] >= 0.95,
        "hp_assign_ge_099": df["hp_assign_rate"] >= 0.99,
        "agreement_positive": df["agreement_positive"],
        "strong_subclone": df["strong_subclone"],
        "lowvaf_highadelta_lowcv": artifact_low_vaf_high_adelta_low_cv(df),
    }


def is_complementary_threshold_pair(name_a: str, name_b: str) -> bool:
    pair = {name_a, name_b}
    if pair in [
        {"pairwise_ge_015", "pairwise_le_015"},
        {"pairwise_ge_020", "pairwise_le_020"},
    ]:
        return True
    return False


def build_feature_outputs(output_dir: Path) -> Dict[str, pd.DataFrame]:
    threshold_rows: List[Dict[str, object]] = []
    interval_rows: List[Dict[str, object]] = []
    top_bin_rows: List[Dict[str, object]] = []
    orth_rows: List[Dict[str, object]] = []
    top_orth_rows: List[Dict[str, object]] = []

    for dataset in RESCUE_DATASETS:
        joined_df = load_joined(dataset.joined_tsv)
        eval_df = candidate_eval_frame(joined_df)
        analyzed_df = eval_df[eval_df["analysis_available"]].copy()

        if analyzed_df.empty:
            continue

        single_metrics: Dict[str, Dict[str, float | int]] = {}
        for rule_def in fine_threshold_defs():
            mask = rule_def["func"](analyzed_df)
            metrics = evaluate_mask(analyzed_df, mask, dataset.baseline_tp, dataset.baseline_fp, dataset.truth_total)
            threshold_rows.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "label": dataset.label,
                    "sample": dataset.sample,
                    "platform": dataset.platform,
                    "mode": dataset.mode,
                    "caller": dataset.caller,
                    "source_scope_note": dataset.source_scope_note,
                    "feature": rule_def["feature"],
                    "rule_id": rule_def["rule_id"],
                    "rule_notes": rule_def["notes"],
                    "trigger_count": metrics["trigger_count"],
                    "tp_rescued": metrics["tp_rescued"],
                    "fp_reintroduced": metrics["fp_reintroduced"],
                    "precision": format_float(float(metrics["precision"])),
                    "recall": format_float(float(metrics["recall"])),
                    "f1": format_float(float(metrics["f1"])),
                    "delta_f1_vs_baseline": format_float(float(metrics["delta_f1_vs_baseline"])),
                    "fp_per_tp": format_float(metrics["fp_reintroduced"] / metrics["tp_rescued"] if metrics["tp_rescued"] else float("inf")),
                }
            )
            single_metrics[rule_def["rule_id"]] = metrics

        for feature, edges in feature_bins().items():
            if feature not in analyzed_df.columns:
                continue
            base = analyzed_df[pd.notna(analyzed_df[feature])].copy()
            if base.empty:
                continue
            base["bin_label"] = cut_series(base[feature], edges)
            total_tp = int((base["downstream_status"] == "caller_lost_tp").sum())
            total_fp = int((base["downstream_status"] == "caller_removed_fp").sum())
            feature_bin_rows: List[Dict[str, object]] = []
            for bin_label, part in base.groupby("bin_label", observed=False):
                tp_count = int((part["downstream_status"] == "caller_lost_tp").sum())
                fp_count = int((part["downstream_status"] == "caller_removed_fp").sum())
                total_count = tp_count + fp_count
                tp_fraction = tp_count / total_tp if total_tp else 0.0
                fp_fraction = fp_count / total_fp if total_fp else 0.0
                tp_fp_ratio = tp_count / fp_count if fp_count else (float("inf") if tp_count else 0.0)
                enrichment = tp_fraction / fp_fraction if fp_fraction else (float("inf") if tp_fraction else 0.0)
                row = {
                    "dataset_id": dataset.dataset_id,
                    "label": dataset.label,
                    "sample": dataset.sample,
                    "platform": dataset.platform,
                    "mode": dataset.mode,
                    "feature": feature,
                    "bin_label": str(bin_label),
                    "tp_count": tp_count,
                    "fp_count": fp_count,
                    "total_count": total_count,
                    "tp_fraction": format_float(tp_fraction),
                    "fp_fraction": format_float(fp_fraction),
                    "tp_to_fp_ratio": format_float(tp_fp_ratio),
                    "enrichment_ratio": format_float(enrichment),
                    "source_scope_note": dataset.source_scope_note,
                }
                interval_rows.append(row)
                feature_bin_rows.append(row)

            filtered = [row for row in feature_bin_rows if row["total_count"] >= 10]
            if filtered:
                best = max(filtered, key=lambda row: (float(row["enrichment_ratio"]) if row["enrichment_ratio"] != "inf" else 1e9, row["tp_count"], -row["fp_count"]))
                top_bin_rows.append(best)

        masks = orthogonal_masks(analyzed_df)
        for name_a, name_b in itertools.combinations(masks.keys(), 2):
            mask_a = masks[name_a]
            mask_b = masks[name_b]
            inter = mask_a & mask_b
            union = mask_a | mask_b
            a_only = mask_a & (~mask_b)
            b_only = mask_b & (~mask_a)
            union_metrics = evaluate_mask(analyzed_df, union, dataset.baseline_tp, dataset.baseline_fp, dataset.truth_total)
            single_a = single_metrics.get(name_a)
            single_b = single_metrics.get(name_b)
            best_single_f1 = max(
                [float(item["f1"]) for item in [single_a, single_b] if item is not None] or [compute_metrics(dataset.baseline_tp, dataset.baseline_fp, dataset.truth_total)["f1"]]
            )
            union_f1 = float(union_metrics["f1"])
            row = {
                "dataset_id": dataset.dataset_id,
                "label": dataset.label,
                "sample": dataset.sample,
                "platform": dataset.platform,
                "mode": dataset.mode,
                "feature_a": name_a,
                "feature_b": name_b,
                "jaccard": format_float(inter.sum() / union.sum() if int(union.sum()) else 0.0),
                "a_only_tp": int(((analyzed_df["downstream_status"] == "caller_lost_tp") & a_only).sum()),
                "a_only_fp": int(((analyzed_df["downstream_status"] == "caller_removed_fp") & a_only).sum()),
                "b_only_tp": int(((analyzed_df["downstream_status"] == "caller_lost_tp") & b_only).sum()),
                "b_only_fp": int(((analyzed_df["downstream_status"] == "caller_removed_fp") & b_only).sum()),
                "both_tp": int(((analyzed_df["downstream_status"] == "caller_lost_tp") & inter).sum()),
                "both_fp": int(((analyzed_df["downstream_status"] == "caller_removed_fp") & inter).sum()),
                "union_tp": int(union_metrics["tp_rescued"]),
                "union_fp": int(union_metrics["fp_reintroduced"]),
                "union_f1": format_float(union_f1),
                "delta_f1_vs_baseline": format_float(float(union_metrics["delta_f1_vs_baseline"])),
                "delta_f1_vs_best_single": format_float(union_f1 - best_single_f1),
                "source_scope_note": dataset.source_scope_note,
            }
            orth_rows.append(row)

        dataset_orth = [row for row in orth_rows if row["dataset_id"] == dataset.dataset_id]
        if dataset_orth:
            dataset_orth = [
                row
                for row in dataset_orth
                if not is_complementary_threshold_pair(str(row["feature_a"]), str(row["feature_b"]))
                and float(row["delta_f1_vs_best_single"]) > 0
            ]
            top_orth_rows.extend(
                sorted(
                    dataset_orth,
                    key=lambda row: (float(row["delta_f1_vs_best_single"]), float(row["delta_f1_vs_baseline"]), -float(row["jaccard"])),
                    reverse=True,
                )[:8]
            )

    return {
        "threshold": pd.DataFrame(threshold_rows),
        "interval": pd.DataFrame(interval_rows),
        "top_bins": pd.DataFrame(top_bin_rows),
        "orthogonality": pd.DataFrame(orth_rows),
        "top_orthogonality": pd.DataFrame(top_orth_rows),
    }


def write_phase2_markdown(
    output_dir: Path,
    paired_df: pd.DataFrame,
    feature_outputs: Dict[str, pd.DataFrame],
) -> None:
    report_path = output_dir / "phase2_summary.md"
    paired_focus = paired_df[
        [
            "label",
            "source_name",
            "tp",
            "fp",
            "fn",
            "precision",
            "recall",
            "f1",
            "delta_f1_vs_merged_output",
        ]
    ].copy()
    top_bins = feature_outputs["top_bins"][
        [
            "dataset_id",
            "feature",
            "bin_label",
            "tp_count",
            "fp_count",
            "tp_to_fp_ratio",
            "enrichment_ratio",
        ]
    ].copy()
    top_orth = feature_outputs["top_orthogonality"][
        [
            "dataset_id",
            "feature_a",
            "feature_b",
            "jaccard",
            "union_tp",
            "union_fp",
            "delta_f1_vs_baseline",
            "delta_f1_vs_best_single",
        ]
    ].copy()

    lines = [
        "# Phase 2：paired raw model 與 finer feature validation",
        "",
        "## 1. paired raw pileup vs full model 直接 benchmark",
        "",
        markdown_table(list(paired_focus.columns), paired_focus.to_dict(orient="records")),
        "",
        "解讀重點：",
        "- `source_name=pileup_filter/full_alignment_filter/merged_output` 代表同一份 paired raw ClairS 輸出的三個 SNV 路徑。",
        "- `delta_f1_vs_merged_output` 為該路徑相對 merged output 的 F1 差值；正值表示優於 merge 後最終 caller 輸出。",
        "",
        "## 2. feature interval 最佳分箱",
        "",
        markdown_table(list(top_bins.columns), top_bins.to_dict(orient="records")),
        "",
        "解讀重點：",
        "- `bin_label` 為特徵分箱區間。",
        "- `tp_to_fp_ratio` 越高表示該箱中的 rescued TP 相對 FP 更集中。",
        "- `enrichment_ratio` 代表 TP fraction / FP fraction，>1 表示該區間相對富集 TP。",
        "",
        "## 3. orthogonality top pairs",
        "",
        markdown_table(list(top_orth.columns), top_orth.to_dict(orient="records")),
        "",
        "解讀重點：",
        "- `jaccard` 越低表示兩特徵重疊越少。",
        "- `delta_f1_vs_best_single` > 0 代表兩特徵聯集比單獨最佳者更好，較接近真正的正交補強。",
    ]
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    ensure_dir(output_dir)

    paired_df = build_paired_raw_benchmark(output_dir, skip_benchmark=args.skip_benchmark)
    feature_outputs = build_feature_outputs(output_dir)

    write_tsv_rows(
        output_dir / "paired_raw_model_benchmark.tsv",
        list(paired_df.columns),
        paired_df.to_dict(orient="records"),
    )
    write_tsv_rows(
        output_dir / "fine_feature_threshold_sweep.tsv",
        list(feature_outputs["threshold"].columns),
        feature_outputs["threshold"].to_dict(orient="records"),
    )
    write_tsv_rows(
        output_dir / "fine_feature_interval_detail.tsv",
        list(feature_outputs["interval"].columns),
        feature_outputs["interval"].to_dict(orient="records"),
    )
    write_tsv_rows(
        output_dir / "fine_feature_interval_top_bins.tsv",
        list(feature_outputs["top_bins"].columns),
        feature_outputs["top_bins"].to_dict(orient="records"),
    )
    write_tsv_rows(
        output_dir / "feature_orthogonality_detail.tsv",
        list(feature_outputs["orthogonality"].columns),
        feature_outputs["orthogonality"].to_dict(orient="records"),
    )
    write_tsv_rows(
        output_dir / "feature_orthogonality_top_pairs.tsv",
        list(feature_outputs["top_orthogonality"].columns),
        feature_outputs["top_orthogonality"].to_dict(orient="records"),
    )

    write_phase2_markdown(output_dir, paired_df, feature_outputs)


if __name__ == "__main__":
    main()
