#!/usr/bin/env python3
"""Integrate HCC1395 5kHz / DORADO paired+TO rescue analyses into a comparable matrix."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd

from research_common import compute_metrics, ensure_dir, markdown_table, parse_variant_counts, write_tsv_rows


REPO_ROOT = Path("/big8_disk/liaoyoyo2001/InterSubMod")
OUTPUT_ROOT_DEFAULT = (
    REPO_ROOT / "output" / "bip8_disk_output" / "research_rounds" / "20260311_hcc1395_cross_quadrant_matrix"
)
MATRIX_ROOT = Path(
    "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
    "20260311_gq_methylation_rescue_matrix_with_dorado_to"
)
TO_DIAG_ROOT = Path(
    "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
    "20260311_to_support_feature_diagnostics"
)


@dataclass
class DatasetConfig:
    dataset_id: str
    label: str
    sample: str
    platform: str
    mode: str
    caller_name: str
    benchmark_comparison_tsv: Optional[Path]
    clairs_variant_counts: Optional[Path]
    longphase_variant_counts: Optional[Path]
    intersubmod_available: bool
    intersubmod_note: str
    raw_pool_tsv: Path
    joined_tsv: Path
    rescue_rule_tsv: Path
    snapshot_summary_tsv: Optional[Path]
    snapshot_bam_scope: str
    snapshot_bam_path: Optional[Path]
    pairwise_support_direction: str
    candidate_specific_round_root: Path


@dataclass
class RawModelConfig:
    sample: str
    platform: str
    mode: str
    caller_family: str
    source_name: str
    caller_dir: Path
    log_path: Path
    final_vcf_path: Path
    pileup_intermediate_path: Optional[Path]
    full_intermediate_path: Optional[Path]
    model_scope: str
    notes: str


DATASETS: List[DatasetConfig] = [
    DatasetConfig(
        dataset_id="hcc1395_5khz_paired",
        label="HCC1395 5kHz paired",
        sample="HCC1395",
        platform="ONT_5kHz",
        mode="paired-pure",
        caller_name="ClairS",
        benchmark_comparison_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
            "20260307_paired_pure_full_columns_official/HCC1395/benchmark_comparison.tsv"
        ),
        clairs_variant_counts=None,
        longphase_variant_counts=None,
        intersubmod_available=True,
        intersubmod_note="paired pure official rerun 已含完整 InterSubMod benchmark",
        raw_pool_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
            "20260308_borderline_rescue_analysis/HCC1395_paired/extract/borderline_candidate_pool.tsv"
        ),
        joined_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
            "20260308_borderline_rescue_analysis/HCC1395_paired/eval/rescue_joined_features.tsv"
        ),
        rescue_rule_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
            "20260308_borderline_rescue_analysis/HCC1395_paired/eval/rescue_rule_comparison.tsv"
        ),
        snapshot_summary_tsv=None,
        snapshot_bam_scope="none",
        snapshot_bam_path=None,
        pairwise_support_direction="not_supportive",
        candidate_specific_round_root=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
            "20260308_borderline_rescue_analysis/HCC1395_paired"
        ),
    ),
    DatasetConfig(
        dataset_id="hcc1395_5khz_to",
        label="HCC1395 5kHz TO",
        sample="HCC1395",
        platform="ONT_5kHz",
        mode="to-pure",
        caller_name="ClairS-TO",
        benchmark_comparison_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
            "20260307_hcc1395_to_pilot_1/benchmark_comparison.tsv"
        ),
        clairs_variant_counts=None,
        longphase_variant_counts=None,
        intersubmod_available=True,
        intersubmod_note="TO pilot 已含完整 InterSubMod benchmark",
        raw_pool_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
            "20260308_borderline_rescue_analysis/HCC1395_to/extract/borderline_candidate_pool.tsv"
        ),
        joined_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
            "20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_joined_features.tsv"
        ),
        rescue_rule_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
            "20260308_hcc1395_to_candidate_rescue_methylation/eval/rescue_rule_comparison.tsv"
        ),
        snapshot_summary_tsv=TO_DIAG_ROOT / "hcc1395_5khz_to" / "snapshot_summary" / "snapshot_summary.tsv",
        snapshot_bam_scope="full_tagged_bam",
        snapshot_bam_path=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
            "20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam"
        ),
        pairwise_support_direction="higher_pairwise_support",
        candidate_specific_round_root=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/"
            "20260308_hcc1395_to_candidate_rescue_methylation"
        ),
    ),
    DatasetConfig(
        dataset_id="hcc1395_dorado_paired",
        label="HCC1395 DORADO paired",
        sample="HCC1395_DORADO",
        platform="ONT_Dorado",
        mode="paired-pure",
        caller_name="ClairS",
        benchmark_comparison_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
            "20260307_paired_pure_full_columns_official/HCC1395_DORADO/benchmark_comparison.tsv"
        ),
        clairs_variant_counts=None,
        longphase_variant_counts=None,
        intersubmod_available=True,
        intersubmod_note="paired pure official rerun 已含完整 InterSubMod benchmark",
        raw_pool_tsv=Path(
            "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/research_rounds/"
            "20260308_borderline_rescue_analysis/HCC1395_DORADO_paired/extract/borderline_candidate_pool.tsv"
        ),
        joined_tsv=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260309_hcc1395_dorado_paired_candidate_rescue/eval/rescue_joined_features.tsv"
        ),
        rescue_rule_tsv=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260309_hcc1395_dorado_paired_candidate_rescue/eval/rescue_rule_comparison.tsv"
        ),
        snapshot_summary_tsv=None,
        snapshot_bam_scope="none",
        snapshot_bam_path=None,
        pairwise_support_direction="lower_pairwise_support",
        candidate_specific_round_root=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260309_hcc1395_dorado_paired_candidate_rescue"
        ),
    ),
    DatasetConfig(
        dataset_id="hcc1395_dorado_to",
        label="HCC1395 DORADO TO",
        sample="HCC1395_DORADO",
        platform="ONT_Dorado",
        mode="to",
        caller_name="ClairS-TO",
        benchmark_comparison_tsv=None,
        clairs_variant_counts=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260311_hcc1395_dorado_to_candidate_rescue/step02_benchmark_clairs_to/variant_counts.txt"
        ),
        longphase_variant_counts=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260311_hcc1395_dorado_to_candidate_rescue/step04_benchmark_longphase_to/variant_counts.txt"
        ),
        intersubmod_available=False,
        intersubmod_note="僅完成 candidate-specific InterSubMod；缺 full baseline InterSubMod benchmark",
        raw_pool_tsv=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260311_hcc1395_dorado_to_candidate_rescue/extract/borderline_candidate_pool.tsv"
        ),
        joined_tsv=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260311_hcc1395_dorado_to_candidate_rescue/eval/rescue_joined_features.tsv"
        ),
        rescue_rule_tsv=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260311_hcc1395_dorado_to_candidate_rescue/eval/rescue_rule_comparison.tsv"
        ),
        snapshot_summary_tsv=TO_DIAG_ROOT / "hcc1395_dorado_to" / "snapshot_summary" / "snapshot_summary.tsv",
        snapshot_bam_scope="candidate_window_subset_tagged_bam",
        snapshot_bam_path=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260311_hcc1395_dorado_to_candidate_rescue/step03_longphase_to/tumor_candidate_windows_tagged.bam"
        ),
        pairwise_support_direction="lower_pairwise_support",
        candidate_specific_round_root=Path(
            "/home/liaoyoyo2001/InterSubMod_runs/output/research_rounds/"
            "20260311_hcc1395_dorado_to_candidate_rescue"
        ),
    ),
]


RAW_MODELS: List[RawModelConfig] = [
    RawModelConfig(
        sample="HCC1395",
        platform="ONT_5kHz",
        mode="paired",
        caller_family="ClairS",
        source_name="ClairS_v0_4_1",
        caller_dir=Path("/big8_disk/data/HCC1395/ONT/ClairS_v0_4_1"),
        log_path=Path("/big8_disk/data/HCC1395/ONT/ClairS_v0_4_1/run_clairs.log"),
        final_vcf_path=Path("/big8_disk/data/HCC1395/ONT/ClairS_v0_4_1/output.vcf.gz"),
        pileup_intermediate_path=Path("/big8_disk/data/HCC1395/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"),
        full_intermediate_path=Path("/big8_disk/data/HCC1395/ONT/ClairS_v0_4_1/tmp/vcf_output/full_alignment_filter.vcf"),
        model_scope="paired_integrated_with_explicit_pileup_and_full_intermediates",
        notes="paired raw caller 同時有 pileup 與 full-alignment 中間產物，但目前主研究結果多來自 final output 或 downstream benchmark",
    ),
    RawModelConfig(
        sample="HCC1395_DORADO",
        platform="ONT_Dorado",
        mode="paired",
        caller_family="ClairS",
        source_name="ClairS_v0_4_1",
        caller_dir=Path("/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_1"),
        log_path=Path("/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_1/run_clairs.log"),
        final_vcf_path=Path("/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_1/output.vcf.gz"),
        pileup_intermediate_path=Path("/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf"),
        full_intermediate_path=Path("/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_1/tmp/vcf_output/full_alignment_filter.vcf"),
        model_scope="paired_integrated_with_explicit_pileup_and_full_intermediates",
        notes="paired raw caller 同時有 pileup 與 full-alignment 中間產物，但尚未完成四象限 full-vs-pileup benchmark 對照",
    ),
    RawModelConfig(
        sample="HCC1395",
        platform="ONT_5kHz",
        mode="to",
        caller_family="ClairS-TO",
        source_name="ClairS_TO_v0_3_0",
        caller_dir=Path("/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0"),
        log_path=Path("/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/run_clairs_to.log"),
        final_vcf_path=Path("/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz"),
        pileup_intermediate_path=None,
        full_intermediate_path=None,
        model_scope="to_pileup_only_ssrs",
        notes="run_clairs_to log 僅宣告 pileup affirmative/negational models，沒有 full-alignment model",
    ),
    RawModelConfig(
        sample="HCC1395",
        platform="ONT_5kHz",
        mode="to",
        caller_family="ClairS-TO",
        source_name="ClairS_TO_ss_v0_3_0",
        caller_dir=Path("/big8_disk/data/HCC1395/ONT/ClairS_TO_ss_v0_3_0"),
        log_path=Path("/big8_disk/data/HCC1395/ONT/ClairS_TO_ss_v0_3_0/run_clairs_to.log"),
        final_vcf_path=Path("/big8_disk/data/HCC1395/ONT/ClairS_TO_ss_v0_3_0/snv.vcf.gz"),
        pileup_intermediate_path=None,
        full_intermediate_path=None,
        model_scope="to_pileup_only_ss",
        notes="ss variant 是不同 model family，不等於 full model；目前仍屬 pileup-only TO 路線",
    ),
    RawModelConfig(
        sample="HCC1395_DORADO",
        platform="ONT_Dorado",
        mode="to",
        caller_family="ClairS-TO",
        source_name="ClairS_TO_v0_3_0",
        caller_dir=Path("/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_v0_3_0"),
        log_path=Path("/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_v0_3_0/run_clairs_to.log"),
        final_vcf_path=Path("/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_v0_3_0/snv.vcf.gz"),
        pileup_intermediate_path=None,
        full_intermediate_path=None,
        model_scope="to_pileup_only_ssrs",
        notes="DORADO TO 現有 candidate-specific 驗證採此路線；目前沒有 full model 對照",
    ),
    RawModelConfig(
        sample="HCC1395_DORADO",
        platform="ONT_Dorado",
        mode="to",
        caller_family="ClairS-TO",
        source_name="ClairS_TO_ss_v0_3_0",
        caller_dir=Path("/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_ss_v0_3_0"),
        log_path=Path("/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_ss_v0_3_0/run_clairs_to.log"),
        final_vcf_path=Path("/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_ss_v0_3_0/snv.vcf.gz"),
        pileup_intermediate_path=None,
        full_intermediate_path=None,
        model_scope="to_pileup_only_ss",
        notes="ss variant 是另一條 pileup-only model 路線；尚未納入本輪主結果 benchmark",
    ),
]


SELECTED_GQ_RULES = {"gq_ge_10", "gq_ge_15", "gq_ge_18", "gq_ge_20", "qual_ge_10_or_gq_ge_20"}
SELECTED_METHYL_RULES = {
    "quality_ge_60",
    "pairwise_ge_015",
    "pairwise_ge_020",
    "pairwise_le_020",
    "hp_assign_ge_099",
    "agreement_positive",
    "strong_subclone",
}
SELECTED_COMBO_RULES = {
    "gq_ge_10__quality_ge_60",
    "gq_ge_10__pairwise_ge_020",
    "gq_ge_10__pairwise_le_020",
    "gq_ge_10__agreement_positive",
    "gq_ge_10__strong_subclone",
    "gq_ge_10__hp_assign_ge_099",
    "gq_ge_15__quality_ge_60",
    "gq_ge_15__pairwise_ge_020",
    "gq_ge_15__pairwise_le_020",
    "gq_ge_15__agreement_positive",
    "gq_ge_15__strong_subclone",
    "gq_ge_15__hp_assign_ge_099",
    "gq_ge_20__quality_ge_60",
    "gq_ge_20__pairwise_ge_020",
    "gq_ge_20__pairwise_le_020",
    "gq_ge_20__agreement_positive",
    "gq_ge_20__strong_subclone",
    "gq_ge_20__hp_assign_ge_099",
}
SELECTED_ARTIFACT_RULES = {"lowvaf_highadelta", "lowvaf_highadelta_lowcv", "combined_artifact_veto"}
KEY_FEATURES = [
    "gq",
    "qual",
    "af",
    "PairwiseMedianDist",
    "Quality_Score",
    "AlleleDelta",
    "CramersV",
    "hp_assign_rate",
    "allele_assign_rate",
]
KEY_ORTHOGONAL_FEATURES = {
    "pairwise_ge_020",
    "pairwise_le_020",
    "quality_ge_60",
    "agreement_positive",
    "strong_subclone",
    "hp_assign_ge_099",
    "lowvaf_highadelta_lowcv",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(OUTPUT_ROOT_DEFAULT), help="Output directory")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def load_matrix_table(name: str) -> pd.DataFrame:
    return read_tsv(MATRIX_ROOT / name)


def metric_row_from_counts(
    dataset: DatasetConfig,
    layer: str,
    method: str,
    counts_path: Path,
    notes: str,
    available: bool = True,
) -> Dict[str, object]:
    counts = parse_variant_counts(counts_path)
    tp = int(counts.get("TP_COUNT", 0))
    fp = int(counts.get("FP_COUNT", 0))
    truth_total = int(counts.get("TRUTH_TOTAL", 0))
    metrics = compute_metrics(tp, fp, truth_total)
    return {
        "dataset_id": dataset.dataset_id,
        "label": dataset.label,
        "sample": dataset.sample,
        "platform": dataset.platform,
        "mode": dataset.mode,
        "layer": layer,
        "method": method,
        "caller_family": dataset.caller_name,
        "available": available,
        "truth_total": metrics["truth_total"],
        "calls_total": metrics["calls_total"],
        "tp": metrics["tp"],
        "fp": metrics["fp"],
        "fn": metrics["fn"],
        "precision": round(metrics["precision"], 6),
        "recall": round(metrics["recall"], 6),
        "f1": round(metrics["f1"], 6),
        "source_path": str(counts_path),
        "notes": notes,
    }


def build_model_availability_matrix() -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for model in RAW_MODELS:
        log_text = model.log_path.read_text(encoding="utf-8", errors="ignore") if model.log_path.exists() else ""
        rows.append(
            {
                "sample": model.sample,
                "platform": model.platform,
                "mode": model.mode,
                "caller_family": model.caller_family,
                "source_name": model.source_name,
                "caller_dir": str(model.caller_dir),
                "log_path": str(model.log_path),
                "final_vcf_path": str(model.final_vcf_path),
                "final_vcf_exists": model.final_vcf_path.exists(),
                "pileup_intermediate_path": str(model.pileup_intermediate_path) if model.pileup_intermediate_path else "",
                "pileup_intermediate_exists": bool(model.pileup_intermediate_path and model.pileup_intermediate_path.exists()),
                "full_intermediate_path": str(model.full_intermediate_path) if model.full_intermediate_path else "",
                "full_intermediate_exists": bool(model.full_intermediate_path and model.full_intermediate_path.exists()),
                "declares_pileup_model": "PILEUP MODEL PATH" in log_text or "PILEUP AFFIRMATIVE MODEL PATH" in log_text,
                "declares_full_model": "FULL-ALIGNMENT MODEL PATH" in log_text or "full_alignment.pkl" in log_text,
                "model_scope": model.model_scope,
                "notes": model.notes,
            }
        )
    return rows


def build_availability_matrix(dataset_overview: pd.DataFrame) -> List[Dict[str, object]]:
    overview_lookup = {row["dataset_id"]: row for _, row in dataset_overview.iterrows()}
    rows: List[Dict[str, object]] = []
    for dataset in DATASETS:
        overview = overview_lookup.get(dataset.dataset_id)
        rows.append(
            {
                "dataset_id": dataset.dataset_id,
                "label": dataset.label,
                "sample": dataset.sample,
                "platform": dataset.platform,
                "mode": dataset.mode,
                "caller_name": dataset.caller_name,
                "benchmark_comparison_path": str(dataset.benchmark_comparison_tsv) if dataset.benchmark_comparison_tsv else "",
                "clairs_variant_counts_path": str(dataset.clairs_variant_counts) if dataset.clairs_variant_counts else "",
                "longphase_variant_counts_path": str(dataset.longphase_variant_counts) if dataset.longphase_variant_counts else "",
                "candidate_round_root": str(dataset.candidate_specific_round_root),
                "raw_pool_tsv": str(dataset.raw_pool_tsv),
                "joined_tsv": str(dataset.joined_tsv),
                "rescue_rule_tsv": str(dataset.rescue_rule_tsv),
                "snapshot_summary_tsv": str(dataset.snapshot_summary_tsv) if dataset.snapshot_summary_tsv else "",
                "snapshot_bam_scope": dataset.snapshot_bam_scope,
                "snapshot_bam_path": str(dataset.snapshot_bam_path) if dataset.snapshot_bam_path else "",
                "candidate_specific_available": dataset.joined_tsv.exists(),
                "snapshot_available": bool(dataset.snapshot_summary_tsv and dataset.snapshot_summary_tsv.exists()),
                "lost_tp_coverage": round(float(overview["lost_tp_coverage"]), 6) if overview is not None else math.nan,
                "removed_fp_coverage": round(float(overview["removed_fp_coverage"]), 6) if overview is not None else math.nan,
                "analyzed_rows": int(overview["analyzed_rows"]) if overview is not None else 0,
                "pairwise_support_direction": dataset.pairwise_support_direction,
                "notes": dataset.intersubmod_note,
            }
        )
    return rows


def build_comparability_matrix(dataset_overview: pd.DataFrame) -> List[Dict[str, object]]:
    overview_lookup = {row["dataset_id"]: row for _, row in dataset_overview.iterrows()}
    rows: List[Dict[str, object]] = []
    for dataset in DATASETS:
        overview = overview_lookup.get(dataset.dataset_id)
        if dataset.dataset_id == "hcc1395_5khz_to":
            warning = "TO read-level diagnostics 使用 full tagged BAM；可做 within-dataset ranking，但不可與 DORADO TO subset snapshot 做絕對值硬比較"
            read_level_comparable = False
        elif dataset.dataset_id == "hcc1395_dorado_to":
            warning = "TO read-level diagnostics 使用 candidate-window subset tagged BAM；只適合觀察局部結構與 presence/absence，不適合與 5kHz TO 做絕對 coverage / hp_assign 比較"
            read_level_comparable = False
        elif dataset.dataset_id == "hcc1395_5khz_paired":
            warning = "paired 有完整 benchmark 與 rescue 矩陣，但目前沒有同規格 TO-support snapshot 診斷"
            read_level_comparable = False
        else:
            warning = "paired 有完整 benchmark 與 rescue 矩陣，但目前沒有同規格 TO-support snapshot 診斷"
            read_level_comparable = False

        rows.append(
            {
                "dataset_id": dataset.dataset_id,
                "label": dataset.label,
                "sample": dataset.sample,
                "platform": dataset.platform,
                "mode": dataset.mode,
                "caller_baseline_available": bool(dataset.benchmark_comparison_tsv or dataset.clairs_variant_counts),
                "longphase_baseline_available": bool(dataset.benchmark_comparison_tsv or dataset.longphase_variant_counts),
                "intersubmod_baseline_available": dataset.intersubmod_available,
                "candidate_specific_available": dataset.joined_tsv.exists(),
                "read_level_diagnostics_available": bool(dataset.snapshot_summary_tsv and dataset.snapshot_summary_tsv.exists()),
                "lost_tp_coverage": round(float(overview["lost_tp_coverage"]), 6) if overview is not None else math.nan,
                "removed_fp_coverage": round(float(overview["removed_fp_coverage"]), 6) if overview is not None else math.nan,
                "benchmark_cross_dataset_comparable": True,
                "candidate_specific_cross_dataset_comparable": bool(
                    overview is not None and float(overview["lost_tp_coverage"]) > 0 and float(overview["removed_fp_coverage"]) > 0
                ),
                "read_level_absolute_cross_dataset_comparable": read_level_comparable,
                "pairwise_direction_stable": False,
                "comparability_warning": warning,
            }
        )
    return rows


def build_layered_performance_matrix() -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for dataset in DATASETS:
        if dataset.benchmark_comparison_tsv and dataset.benchmark_comparison_tsv.exists():
            df = read_tsv(dataset.benchmark_comparison_tsv)
            for _, row in df.iterrows():
                rows.append(
                    {
                        "dataset_id": dataset.dataset_id,
                        "label": dataset.label,
                        "sample": dataset.sample,
                        "platform": dataset.platform,
                        "mode": dataset.mode,
                        "layer": {
                            "ClairS": "caller_baseline",
                            "ClairS-TO": "caller_baseline",
                            "LongPhase-S": "longphase_baseline",
                            "LongPhase-TO": "longphase_baseline",
                            "InterSubMod": "intersubmod_baseline",
                        }.get(str(row["method"]), "unknown"),
                        "method": row["method"],
                        "caller_family": row["caller"],
                        "available": True,
                        "truth_total": int(row["truth_total"]),
                        "calls_total": int(row["calls_total"]),
                        "tp": int(row["tp"]),
                        "fp": int(row["fp"]),
                        "fn": int(row["fn"]),
                        "precision": float(row["precision"]),
                        "recall": float(row["recall"]),
                        "f1": float(row["f1"]),
                        "source_path": str(dataset.benchmark_comparison_tsv),
                        "notes": str(row["notes"]),
                    }
                )
            continue

        if dataset.clairs_variant_counts:
            rows.append(
                metric_row_from_counts(
                    dataset,
                    layer="caller_baseline",
                    method="ClairS-TO",
                    counts_path=dataset.clairs_variant_counts,
                    notes="由 variant_counts.txt 重建 DORADO TO baseline",
                )
            )
        if dataset.longphase_variant_counts:
            rows.append(
                metric_row_from_counts(
                    dataset,
                    layer="longphase_baseline",
                    method="LongPhase-TO",
                    counts_path=dataset.longphase_variant_counts,
                    notes="由 variant_counts.txt 重建 DORADO TO LongPhase-TO baseline",
                )
            )
        rows.append(
            {
                "dataset_id": dataset.dataset_id,
                "label": dataset.label,
                "sample": dataset.sample,
                "platform": dataset.platform,
                "mode": dataset.mode,
                "layer": "intersubmod_baseline",
                "method": "InterSubMod",
                "caller_family": dataset.caller_name,
                "available": dataset.intersubmod_available,
                "truth_total": "",
                "calls_total": "",
                "tp": "",
                "fp": "",
                "fn": "",
                "precision": "",
                "recall": "",
                "f1": "",
                "source_path": "",
                "notes": dataset.intersubmod_note,
            }
        )
    return rows


def build_candidate_pool_layer_matrix(dataset_overview: pd.DataFrame) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for _, row in dataset_overview.iterrows():
        rows.extend(
            [
                {
                    "dataset_id": row["dataset_id"],
                    "label": row["label"],
                    "sample": row["sample"],
                    "platform": row["platform"],
                    "mode": row["mode"],
                    "layer": "candidate_pool",
                    "subset": "caller_lost_tp_raw",
                    "rows_total": int(row["raw_lost_tp_total"]),
                    "rows_with_joined_features": int(row["joined_lost_tp_rows"]),
                    "rows_analyzed": int(row["analyzed_lost_tp_rows"]),
                    "coverage_fraction": round(float(row["lost_tp_coverage"]), 6),
                },
                {
                    "dataset_id": row["dataset_id"],
                    "label": row["label"],
                    "sample": row["sample"],
                    "platform": row["platform"],
                    "mode": row["mode"],
                    "layer": "candidate_pool",
                    "subset": "caller_removed_fp_raw",
                    "rows_total": int(row["raw_removed_fp_total"]),
                    "rows_with_joined_features": int(row["joined_removed_fp_rows"]),
                    "rows_analyzed": int(row["analyzed_removed_fp_rows"]),
                    "coverage_fraction": round(float(row["removed_fp_coverage"]), 6),
                },
            ]
        )
    return rows


def build_feature_distribution_focus() -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for dataset in DATASETS:
        df = read_tsv(dataset.joined_tsv)
        if "candidate_eligible" in df.columns:
            eligible = df["candidate_eligible"].astype(str).str.lower().isin({"true", "1", "yes", "y"})
            df = df[eligible].copy()
        if "downstream_status" in df.columns:
            df = df[df["downstream_status"].isin({"caller_lost_tp", "caller_removed_fp"})].copy()
        for feature in KEY_FEATURES:
            if feature not in df.columns:
                continue
            numeric = pd.to_numeric(df[feature], errors="coerce")
            for truth_status in ("TP", "FP"):
                sub = numeric[df["truth_status"] == truth_status].dropna()
                if sub.empty:
                    continue
                rows.append(
                    {
                        "dataset_id": dataset.dataset_id,
                        "label": dataset.label,
                        "sample": dataset.sample,
                        "platform": dataset.platform,
                        "mode": dataset.mode,
                        "feature": feature,
                        "truth_status": truth_status,
                        "n": int(sub.shape[0]),
                        "mean": round(float(sub.mean()), 6),
                        "median": round(float(sub.median()), 6),
                        "p25": round(float(sub.quantile(0.25)), 6),
                        "p75": round(float(sub.quantile(0.75)), 6),
                        "p10": round(float(sub.quantile(0.10)), 6),
                        "p90": round(float(sub.quantile(0.90)), 6),
                        "source_path": str(dataset.joined_tsv),
                    }
                )
    return rows


def load_selected_rows(df: pd.DataFrame, filters: Dict[str, Iterable[str]]) -> pd.DataFrame:
    work = df.copy()
    for column, values in filters.items():
        work = work[work[column].isin(list(values))]
    return work.copy()


def build_selected_rule_matrix() -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    gq_df = load_matrix_table("gq_threshold_sweep.tsv")
    for _, row in load_selected_rows(gq_df, {"rule_id": SELECTED_GQ_RULES}).iterrows():
        rows.append(
            {
                "dataset_id": row["dataset_id"],
                "label": row["label"],
                "sample": row["sample"],
                "platform": row["platform"],
                "mode": row["mode"],
                "rule_group": "caller_only",
                "rule_id": row["rule_id"],
                "feature_family": "gq",
                "gate_id": "",
                "trigger_count": int(row["trigger_count"]),
                "tp_rescued": int(row["tp_rescued"]),
                "fp_reintroduced": int(row["fp_reintroduced"]),
                "precision": float(row["precision"]),
                "recall": float(row["recall"]),
                "f1": float(row["f1"]),
                "delta_f1_vs_baseline": float(row["delta_f1_vs_baseline"]),
                "delta_f1_vs_gate": math.nan,
                "fp_per_tp": float(row["fp_per_tp"]) if str(row["fp_per_tp"]) != "inf" else math.inf,
                "meets_safety": row["meets_safety"],
                "source_path": str(MATRIX_ROOT / "gq_threshold_sweep.tsv"),
                "notes": f"threshold_type={row['threshold_type']}, threshold_value={row['threshold_value']}",
            }
        )

    methyl_df = load_matrix_table("methylation_only_rule_sweep.tsv")
    for _, row in load_selected_rows(methyl_df, {"rule_id": SELECTED_METHYL_RULES}).iterrows():
        rows.append(
            {
                "dataset_id": row["dataset_id"],
                "label": row["label"],
                "sample": row["sample"],
                "platform": row["platform"],
                "mode": row["mode"],
                "rule_group": "methylation_only",
                "rule_id": row["rule_id"],
                "feature_family": row["feature_family"],
                "gate_id": "",
                "trigger_count": int(row["trigger_count"]),
                "tp_rescued": int(row["tp_rescued"]),
                "fp_reintroduced": int(row["fp_reintroduced"]),
                "precision": float(row["precision"]),
                "recall": float(row["recall"]),
                "f1": float(row["f1"]),
                "delta_f1_vs_baseline": float(row["delta_f1_vs_baseline"]),
                "delta_f1_vs_gate": math.nan,
                "fp_per_tp": float(row["fp_per_tp"]) if str(row["fp_per_tp"]) != "inf" else math.inf,
                "meets_safety": row["meets_safety"],
                "source_path": str(MATRIX_ROOT / "methylation_only_rule_sweep.tsv"),
                "notes": row["rule_notes"],
            }
        )

    combo_df = load_matrix_table("gq_plus_methylation_rule_sweep.tsv")
    for _, row in load_selected_rows(combo_df, {"rule_id": SELECTED_COMBO_RULES}).iterrows():
        rows.append(
            {
                "dataset_id": row["dataset_id"],
                "label": row["label"],
                "sample": row["sample"],
                "platform": row["platform"],
                "mode": row["mode"],
                "rule_group": "caller_plus_methylation",
                "rule_id": row["rule_id"],
                "feature_family": row["feature_family"],
                "gate_id": row["gate_id"],
                "trigger_count": int(row["trigger_count"]),
                "tp_rescued": int(row["tp_rescued"]),
                "fp_reintroduced": int(row["fp_reintroduced"]),
                "precision": float(row["precision"]),
                "recall": float(row["recall"]),
                "f1": float(row["f1"]),
                "delta_f1_vs_baseline": float(row["delta_f1_vs_baseline"]),
                "delta_f1_vs_gate": float(row["delta_f1_vs_gate"]),
                "fp_per_tp": float(row["fp_per_tp"]) if str(row["fp_per_tp"]) != "inf" else math.inf,
                "meets_safety": row["meets_safety"],
                "source_path": str(MATRIX_ROOT / "gq_plus_methylation_rule_sweep.tsv"),
                "notes": row["rule_notes"],
            }
        )

    artifact_df = load_matrix_table("artifact_role_comparison.tsv")
    selected_artifact = artifact_df[
        artifact_df["artifact_rule"].isin(SELECTED_ARTIFACT_RULES)
        & artifact_df["role"].eq("rescue_veto")
        & artifact_df["base_gate"].isin({"gq_ge_10", "gq_ge_15", "gq_ge_20", "qual_ge_10_or_gq_ge_20"})
    ]
    for _, row in selected_artifact.iterrows():
        rows.append(
            {
                "dataset_id": row["dataset_id"],
                "label": row["label"],
                "sample": row["sample"],
                "platform": row["platform"],
                "mode": row["mode"],
                "rule_group": "artifact_veto",
                "rule_id": row["artifact_rule"],
                "feature_family": "artifact",
                "gate_id": row["base_gate"],
                "trigger_count": int(row["tp_selected"]) + int(row["fp_selected"]),
                "tp_rescued": int(row["tp_selected"]),
                "fp_reintroduced": int(row["fp_selected"]),
                "precision": float(row["precision"]),
                "recall": float(row["recall"]),
                "f1": float(row["f1"]),
                "delta_f1_vs_baseline": float(row["delta_f1_vs_baseline"]),
                "delta_f1_vs_gate": float(row["delta_f1_vs_gate"]),
                "fp_per_tp": float(row["fp_selected"]) / float(row["tp_selected"]) if float(row["tp_selected"]) else math.inf,
                "meets_safety": float(row["fp_selected"]) <= float(row["tp_selected"]) if float(row["tp_selected"]) else False,
                "source_path": str(MATRIX_ROOT / "artifact_role_comparison.tsv"),
                "notes": f"tp_removed_by_veto={row['tp_removed_by_veto']}; fp_removed_by_veto={row['fp_removed_by_veto']}",
            }
        )

    return rows


def build_feature_overlap_focus() -> List[Dict[str, object]]:
    orth_df = load_matrix_table("orthogonal_feature_candidates.tsv")
    rows: List[Dict[str, object]] = []
    for _, row in orth_df.iterrows():
        if row["feature_a"] not in KEY_ORTHOGONAL_FEATURES and row["feature_b"] not in KEY_ORTHOGONAL_FEATURES:
            continue
        rows.append({key: row[key] for key in orth_df.columns})
    return rows


def build_snapshot_scope_bias_assessment() -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for dataset in DATASETS:
        if not dataset.snapshot_summary_tsv or not dataset.snapshot_summary_tsv.exists():
            rows.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "label": dataset.label,
                    "mode": dataset.mode,
                    "diagnostics_available": False,
                    "snapshot_bam_scope": dataset.snapshot_bam_scope,
                    "snapshot_bam_path": str(dataset.snapshot_bam_path) if dataset.snapshot_bam_path else "",
                    "regions_profiled": 0,
                    "median_target_depth": "",
                    "median_target_alt_fraction": "",
                    "median_na_hp_fraction": "",
                    "high_na_hp_regions": "",
                    "haplotype_skewed_regions": "",
                    "cross_dataset_absolute_comparable": False,
                    "within_dataset_interpretation": "no_standardized_snapshot",
                    "notes": "本輪沒有同規格 read-level support diagnostics；不可與 TO snapshot 直接類比",
                }
            )
            continue

        df = read_tsv(dataset.snapshot_summary_tsv)
        notes_text = df["notes"].fillna("").astype(str)
        rows.append(
            {
                "dataset_id": dataset.dataset_id,
                "label": dataset.label,
                "mode": dataset.mode,
                "diagnostics_available": True,
                "snapshot_bam_scope": dataset.snapshot_bam_scope,
                "snapshot_bam_path": str(dataset.snapshot_bam_path) if dataset.snapshot_bam_path else "",
                "regions_profiled": int(df.shape[0]),
                "median_target_depth": round(float(pd.to_numeric(df["target_depth"], errors="coerce").median()), 6),
                "median_target_alt_fraction": round(float(pd.to_numeric(df["target_alt_fraction"], errors="coerce").median()), 6),
                "median_na_hp_fraction": round(float(pd.to_numeric(df["na_hp_fraction"], errors="coerce").median()), 6),
                "high_na_hp_regions": int(notes_text.str.contains("high_na_hp").sum()),
                "haplotype_skewed_regions": int(notes_text.str.contains("haplotype_skewed").sum()),
                "cross_dataset_absolute_comparable": False,
                "within_dataset_interpretation": "safe_for_within_dataset_ranking_only",
                "notes": (
                    "5kHz TO 為 full tagged BAM、DORADO TO 為 candidate-window subset tagged BAM；"
                    "可比較 support/artefact 現象方向，但不可把 coverage、hp_assign 或 alt_fraction 當跨 dataset 絕對值比較"
                ),
            }
        )
    return rows


def build_phase2_gap_tracker() -> List[Dict[str, object]]:
    return [
        {
            "gap_id": "snapshot_scope_alignment",
            "priority": "high",
            "affected_scope": "TO read-level diagnostics",
            "blocking_conclusion": "5kHz TO 與 DORADO TO 的 read-level 絕對值差異",
            "required_action": "若要做 cross-dataset 絕對比較，需在同一 dataset 上同時建立 full 與 subset snapshot 對照，或補 DORADO TO full tagged BAM 驗證",
        },
        {
            "gap_id": "paired_support_diagnostics",
            "priority": "high",
            "affected_scope": "5kHz paired / DORADO paired",
            "blocking_conclusion": "paired support 特徵 read-level 現象是否與 TO 一致",
            "required_action": "用與 TO 相同規格，補 paired 的 Quality_Score / Pairwise / hp_assign diagnostics",
        },
        {
            "gap_id": "to_full_model_comparison",
            "priority": "high",
            "affected_scope": "TO pileup vs full model",
            "blocking_conclusion": "TO 下 pileup/full model 是否改變 rescue 特徵方向",
            "required_action": "先盤點是否真有 full model TO 可 benchmark；若沒有，只能維持 availability 結論",
        },
        {
            "gap_id": "paired_direct_pileup_full_benchmark",
            "priority": "medium",
            "affected_scope": "paired pileup vs full",
            "blocking_conclusion": "paired 下 pileup_filter / full_alignment_filter 的直接層級差異",
            "required_action": "以 v0.4.1 的 pileup/full intermediates 建立一致 benchmark split，再重算 layered matrix",
        },
        {
            "gap_id": "full_baseline_intersubmod_dorado_to",
            "priority": "medium",
            "affected_scope": "HCC1395_DORADO TO",
            "blocking_conclusion": "DORADO TO 完整 baseline 的 ClairS-TO -> LongPhase-TO -> InterSubMod 層級比較",
            "required_action": "若空間許可，補 full baseline InterSubMod；若不補，報告中持續標記 candidate-specific only",
        },
    ]


def format_float(value: object, digits: int = 4) -> str:
    try:
        num = float(value)
    except (TypeError, ValueError):
        return ""
    if math.isnan(num):
        return ""
    return f"{num:.{digits}f}"


def build_final_report(
    output_dir: Path,
    model_rows: List[Dict[str, object]],
    availability_rows: List[Dict[str, object]],
    comparability_rows: List[Dict[str, object]],
    layered_rows: List[Dict[str, object]],
    feature_distribution_rows: List[Dict[str, object]],
    selected_rule_rows: List[Dict[str, object]],
    overlap_rows: List[Dict[str, object]],
    snapshot_rows: List[Dict[str, object]],
    phase2_rows: List[Dict[str, object]],
) -> str:
    availability_path = output_dir / "availability_matrix.tsv"
    comparability_path = output_dir / "comparability_matrix.tsv"
    model_path = output_dir / "model_availability_matrix.tsv"
    layered_path = output_dir / "layered_performance_matrix.tsv"
    candidate_pool_path = output_dir / "candidate_pool_layer_matrix.tsv"
    feature_dist_path = output_dir / "feature_distribution_focus.tsv"
    selected_rules_path = output_dir / "selected_rule_matrix.tsv"
    overlap_path = output_dir / "feature_overlap_focus.tsv"
    snapshot_path = output_dir / "snapshot_scope_bias_assessment.tsv"
    phase2_path = output_dir / "phase2_gap_tracker.tsv"

    layered_df = pd.DataFrame(layered_rows)
    key_layer_rows = []
    for dataset_id in [d.dataset_id for d in DATASETS]:
        sub = layered_df[(layered_df["dataset_id"] == dataset_id) & (layered_df["available"] == True)].copy()
        if sub.empty:
            continue
        sub = sub[["label", "method", "tp", "fp", "fn", "f1", "notes", "source_path"]]
        key_layer_rows.extend(sub.to_dict("records"))

    selected_df = pd.DataFrame(selected_rule_rows)
    rule_focus = []
    for dataset_id in [d.dataset_id for d in DATASETS]:
        sub = selected_df[selected_df["dataset_id"] == dataset_id].copy()
        if sub.empty:
            continue
        for rule_group, rule_id in [
            ("caller_only", "gq_ge_10"),
            ("caller_only", "gq_ge_15"),
            ("caller_only", "gq_ge_20"),
            ("methylation_only", "quality_ge_60"),
            ("methylation_only", "pairwise_ge_020"),
            ("methylation_only", "pairwise_le_020"),
            ("methylation_only", "hp_assign_ge_099"),
            ("methylation_only", "agreement_positive"),
            ("methylation_only", "strong_subclone"),
            ("caller_plus_methylation", "gq_ge_10__quality_ge_60"),
            ("caller_plus_methylation", "gq_ge_10__pairwise_ge_020"),
            ("caller_plus_methylation", "gq_ge_10__pairwise_le_020"),
            ("caller_plus_methylation", "gq_ge_10__agreement_positive"),
            ("caller_plus_methylation", "gq_ge_10__strong_subclone"),
            ("caller_plus_methylation", "gq_ge_10__hp_assign_ge_099"),
        ]:
            row = sub[(sub["rule_group"] == rule_group) & (sub["rule_id"] == rule_id)]
            if row.empty:
                continue
            rule_focus.extend(row.to_dict("records"))

    feature_df = pd.DataFrame(feature_distribution_rows)
    dist_focus = feature_df[
        feature_df["feature"].isin(["gq", "qual", "PairwiseMedianDist", "Quality_Score", "AlleleDelta", "hp_assign_rate"])
    ].copy()

    overlap_df = pd.DataFrame(overlap_rows)
    overlap_focus = overlap_df[
        (
            overlap_df["feature_a"].isin(["pairwise_ge_020", "pairwise_le_020", "quality_ge_60", "agreement_positive", "hp_assign_ge_099"])
            | overlap_df["feature_b"].isin(["pairwise_ge_020", "pairwise_le_020", "quality_ge_60", "agreement_positive", "hp_assign_ge_099"])
        )
    ].head(16)

    snapshot_df = pd.DataFrame(snapshot_rows)
    comparability_df = pd.DataFrame(comparability_rows)

    lines: List[str] = []
    lines.append("# HCC1395 四象限甲基 Rescue / Feature Matrix 整合報告")
    lines.append("")
    lines.append("## 1. 破題結論")
    lines.append("")
    lines.append(
        "這份整合報告把目前 `HCC1395 5kHz / HCC1395_DORADO × paired / TO` 的 caller、LongPhase、InterSubMod、candidate-specific rescue、甲基特徵矩陣與 TO read-level diagnostics 放到同一個可比框架內。"
    )
    lines.append("")
    lines.append("目前能下的最穩定結論有 5 個：")
    lines.append("")
    lines.append("1. `caller-first` 仍是四象限中最穩定的主規則，特別是 `GQ`。")
    lines.append(
        "2. 甲基特徵不是全部負效果；`Quality_Score`、`PairwiseMedianDist`、`agreement_positive`、`hp_assign_rate` 在部分 dataset 對 baseline 有正向 rescue 訊號，但多半沒有超過同一 caller gate。"
    )
    lines.append(
        "3. `PairwiseMedianDist` 方向具有明顯 dataset-dependence：`5kHz TO` 偏高 pairwise 是 support，`DORADO paired/TO` 偏低 pairwise 較合理，因此目前不能升級成全域規則。"
    )
    lines.append(
        "4. `low VAF + high AlleleDelta (+ low CramersV)` 仍較適合後段 artifact triage，不適合提前當 TP rescue 主規則。"
    )
    lines.append(
        "5. `5kHz TO` 的 snapshot 來自 full tagged BAM，而 `DORADO TO` 來自 candidate-window subset tagged BAM；因此 TO 的 read-level 圖與 HP 結構只能做方向比較，不能直接拿 coverage、alt fraction、hp_assign 的絕對值硬比。"
    )
    lines.append("")
    lines.append("## 2. 本輪新輸出")
    lines.append("")
    for path in [
        model_path,
        availability_path,
        comparability_path,
        layered_path,
        candidate_pool_path,
        feature_dist_path,
        selected_rules_path,
        overlap_path,
        snapshot_path,
        phase2_path,
    ]:
        lines.append(f"- [{path.name}]({path})")
    lines.append("")
    lines.append("## 3. 四象限資料與模型可比性")
    lines.append("")
    lines.append("### 3.1 Raw caller / model availability")
    lines.append("")
    lines.append(markdown_table(
        ["sample", "platform", "mode", "caller_family", "source_name", "declares_pileup_model", "declares_full_model", "model_scope"],
        model_rows,
    ))
    lines.append("")
    lines.append("重點：")
    lines.append("")
    lines.append("- paired 的 `ClairS_v0_4_1` 明確同時有 `pileup_filter.vcf` 與 `full_alignment_filter.vcf`，表示 raw caller 層確實能做 pileup/full model availability 盤點。")
    lines.append("- TO 的 `ClairS_TO_v0_3_0` 與 `ClairS_TO_ss_v0_3_0` 目前都只看到 pileup affirmative/negational model；本輪不能把 TO 當成已有 full model 直接對照。")
    lines.append("")
    lines.append("### 3.2 Comparability warning")
    lines.append("")
    lines.append(markdown_table(
        [
            "label",
            "mode",
            "caller_baseline_available",
            "longphase_baseline_available",
            "intersubmod_baseline_available",
            "candidate_specific_available",
            "read_level_diagnostics_available",
            "lost_tp_coverage",
            "removed_fp_coverage",
            "read_level_absolute_cross_dataset_comparable",
            "comparability_warning",
        ],
        comparability_rows,
    ))
    lines.append("")
    lines.append("## 4. 分層 benchmark 結果")
    lines.append("")
    lines.append(markdown_table(
        ["label", "method", "tp", "fp", "fn", "f1", "notes"],
        key_layer_rows,
    ))
    lines.append("")
    lines.append("重點：")
    lines.append("")
    lines.append("- `HCC1395 5kHz paired`：`InterSubMod F1=0.8532`，高於 `LongPhase-S 0.8522`。")
    lines.append("- `HCC1395_DORADO paired`：`InterSubMod F1=0.8590`，略低於 `LongPhase-S 0.8592`。")
    lines.append("- `HCC1395 5kHz TO`：`ClairS-TO 0.7127 -> LongPhase-TO 0.7127 -> InterSubMod 0.7130`，TO 主正增益來自 InterSubMod 過濾層。")
    lines.append("- `HCC1395_DORADO TO`：目前有 `ClairS-TO` 與 `LongPhase-TO` baseline，但缺 full baseline InterSubMod benchmark，只能在 candidate-specific 層下結論。")
    lines.append("")
    lines.append("## 5. Candidate-specific coverage")
    lines.append("")
    lines.append(markdown_table(
        ["dataset_id", "subset", "rows_total", "rows_with_joined_features", "rows_analyzed", "coverage_fraction"],
        read_tsv(candidate_pool_path).to_dict("records"),
    ))
    lines.append("")
    lines.append("重點：")
    lines.append("")
    lines.append("- `5kHz TO` 與 `DORADO TO` 的 candidate-specific coverage 足夠高，TO rescue 結論可信。")
    lines.append("- `5kHz paired` 與 `DORADO paired` 的 lost_tp/removed_fp coverage 明顯較低，因此 paired 的甲基 rescue 更容易受 coverage ceiling 限制。")
    lines.append("")
    lines.append("## 6. 關鍵規則比較")
    lines.append("")
    lines.append(markdown_table(
        [
            "dataset_id",
            "rule_group",
            "rule_id",
            "gate_id",
            "tp_rescued",
            "fp_reintroduced",
            "delta_f1_vs_baseline",
            "delta_f1_vs_gate",
            "notes",
        ],
        rule_focus,
    ))
    lines.append("")
    lines.append("重點：")
    lines.append("")
    lines.append("- `GQ` 仍是最穩定的 caller-only 規則來源，但最佳 threshold 依 dataset 不同。")
    lines.append("- `Quality_Score>=60` 在 `5kHz TO` 與 `DORADO TO` 對 baseline 都是正向，但固定 `gq>=10` gate 後未超過 caller-only。")
    lines.append("- `agreement_positive` 與 `Strong/Subclone` 在 `5kHz TO` 屬較乾淨的 support，但在 `DORADO paired/TO` 沒有形成跨樣本主規則。")
    lines.append("- `hp_assign_rate>=0.99` 在 `DORADO paired/TO` 有小幅正訊號，但 read-level diagnostics 顯示它更接近 phase/QC annotation，而不是 truth-discriminative 主特徵。")
    lines.append("")
    lines.append("## 7. 特徵分佈與 dataset-dependent 方向")
    lines.append("")
    lines.append(markdown_table(
        ["dataset_id", "feature", "truth_status", "n", "median", "p25", "p75", "source_path"],
        dist_focus.to_dict("records")[:32],
    ))
    lines.append("")
    lines.append("解讀方式：")
    lines.append("")
    lines.append("- `gq / qual`：觀察 caller 邊緣候選在 TP/FP 間是否存在穩定分離。")
    lines.append("- `PairwiseMedianDist`：觀察結構異質性方向是否隨 dataset 改變。")
    lines.append("- `Quality_Score`：觀察甲基訊號整體品質是否能當 soft support。")
    lines.append("- `AlleleDelta` 與 `CramersV`：主要用於 artifact triage，而非主 rescue。")
    lines.append("- `hp_assign_rate`：需配合 snapshot scope 解讀，不能把 full-tagged 與 subset-tagged 的絕對值直接混比。")
    lines.append("- paired 兩個 dataset 在 candidate-specific 可分析 coverage 較低，許多甲基欄位在 joined table 中會呈現大量 `0/NA`；這表示目前 paired 的甲基分佈比較適合解讀為「coverage ceiling 與方向性不足」，而不是直接當成甲基完全無訊號。")
    lines.append("")
    lines.append("## 8. 交集與正交性")
    lines.append("")
    if not overlap_focus.empty:
        lines.append(markdown_table(
            ["dataset_id", "feature_a", "feature_b", "jaccard", "a_only_tp", "a_only_fp", "b_only_tp", "b_only_fp", "both_tp", "both_fp"],
            overlap_focus.to_dict("records"),
        ))
    else:
        lines.append("_No overlap rows_")
    lines.append("")
    lines.append("解讀方式：")
    lines.append("")
    lines.append("- `jaccard` 越低，代表兩個 feature 較正交。")
    lines.append("- `a_only_tp / b_only_tp` 若明顯大於 `a_only_fp / b_only_fp`，代表可能存在互補 support。")
    lines.append("- `both_fp` 高時，代表交集規則未必更乾淨，可能只是把同一群 noisy candidate 重複描述。")
    lines.append("")
    lines.append("## 9. Snapshot scope 與 read-level diagnostics")
    lines.append("")
    lines.append(markdown_table(
        [
            "dataset_id",
            "snapshot_bam_scope",
            "regions_profiled",
            "median_target_depth",
            "median_target_alt_fraction",
            "median_na_hp_fraction",
            "cross_dataset_absolute_comparable",
            "notes",
        ],
        snapshot_rows,
    ))
    lines.append("")
    lines.append("重點：")
    lines.append("")
    lines.append("- `5kHz TO` snapshot 來自 full tagged BAM，能較完整反映全體 read 結構。")
    lines.append("- `DORADO TO` snapshot 來自 candidate-window subset tagged BAM，較適合觀察局部 support / artifact 現象，但不適合與 5kHz TO 做絕對深度與 HP 指標比較。")
    lines.append("- 因此 TO diagnostics 目前只能回答「這個特徵在該 dataset 的 read-level 現象長什麼樣」，不能回答「兩個 dataset 的 read-level 數值誰更高」。")
    lines.append("")
    lines.append("## 10. 目前最合理的高層結論")
    lines.append("")
    lines.append("1. `caller-first` 仍是主流程核心，四象限都沒有證據支持甲基取代 caller。")
    lines.append("2. `methylation-support` 目前最合理的定位是：")
    lines.append("   - `Quality_Score`: soft support / ranking / annotation")
    lines.append("   - `PairwiseMedianDist`: dataset-aware annotation，不可全域固定方向")
    lines.append("   - `agreement_positive`: 只在部分 dataset 提供較乾淨的 support 子集")
    lines.append("   - `hp_assign_rate`: phase/QC annotation，而非主判決規則")
    lines.append("3. `artifact triage` 仍以 `low VAF + high AlleleDelta (+ low CramersV)` 為核心，但應放在後段 FP triage，不該拿來當主 rescue 規則。")
    lines.append("4. `pileup / full model` 現階段能做的是 availability 與 phase 2 規劃，還不能直接宣稱四象限都已有完整對照。")
    lines.append("")
    lines.append("## 11. Phase 2 缺口")
    lines.append("")
    lines.append(markdown_table(
        ["gap_id", "priority", "affected_scope", "blocking_conclusion", "required_action"],
        phase2_rows,
    ))
    lines.append("")
    lines.append("## 12. 主要可複查路徑")
    lines.append("")
    lines.append(f"- [model_availability_matrix.tsv]({model_path})")
    lines.append(f"- [availability_matrix.tsv]({availability_path})")
    lines.append(f"- [comparability_matrix.tsv]({comparability_path})")
    lines.append(f"- [layered_performance_matrix.tsv]({layered_path})")
    lines.append(f"- [feature_distribution_focus.tsv]({feature_dist_path})")
    lines.append(f"- [selected_rule_matrix.tsv]({selected_rules_path})")
    lines.append(f"- [feature_overlap_focus.tsv]({overlap_path})")
    lines.append(f"- [snapshot_scope_bias_assessment.tsv]({snapshot_path})")
    lines.append(f"- [phase2_gap_tracker.tsv]({phase2_path})")
    lines.append(f"- [既有 GQ / 甲基矩陣總表]({MATRIX_ROOT / 'analysis_summary.md'})")
    lines.append(f"- [TO support read-level diagnostics]({TO_DIAG_ROOT})")
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    output_dir = ensure_dir(Path(args.output_dir).resolve())

    dataset_overview = load_matrix_table("dataset_overview.tsv")
    availability_rows = build_availability_matrix(dataset_overview)
    comparability_rows = build_comparability_matrix(dataset_overview)
    model_rows = build_model_availability_matrix()
    layered_rows = build_layered_performance_matrix()
    candidate_pool_rows = build_candidate_pool_layer_matrix(dataset_overview)
    feature_distribution_rows = build_feature_distribution_focus()
    selected_rule_rows = build_selected_rule_matrix()
    overlap_rows = build_feature_overlap_focus()
    snapshot_rows = build_snapshot_scope_bias_assessment()
    phase2_rows = build_phase2_gap_tracker()

    write_tsv_rows(
        output_dir / "model_availability_matrix.tsv",
        [
            "sample",
            "platform",
            "mode",
            "caller_family",
            "source_name",
            "caller_dir",
            "log_path",
            "final_vcf_path",
            "final_vcf_exists",
            "pileup_intermediate_path",
            "pileup_intermediate_exists",
            "full_intermediate_path",
            "full_intermediate_exists",
            "declares_pileup_model",
            "declares_full_model",
            "model_scope",
            "notes",
        ],
        model_rows,
    )
    write_tsv_rows(
        output_dir / "availability_matrix.tsv",
        [
            "dataset_id",
            "label",
            "sample",
            "platform",
            "mode",
            "caller_name",
            "benchmark_comparison_path",
            "clairs_variant_counts_path",
            "longphase_variant_counts_path",
            "candidate_round_root",
            "raw_pool_tsv",
            "joined_tsv",
            "rescue_rule_tsv",
            "snapshot_summary_tsv",
            "snapshot_bam_scope",
            "snapshot_bam_path",
            "candidate_specific_available",
            "snapshot_available",
            "lost_tp_coverage",
            "removed_fp_coverage",
            "analyzed_rows",
            "pairwise_support_direction",
            "notes",
        ],
        availability_rows,
    )
    write_tsv_rows(
        output_dir / "comparability_matrix.tsv",
        [
            "dataset_id",
            "label",
            "sample",
            "platform",
            "mode",
            "caller_baseline_available",
            "longphase_baseline_available",
            "intersubmod_baseline_available",
            "candidate_specific_available",
            "read_level_diagnostics_available",
            "lost_tp_coverage",
            "removed_fp_coverage",
            "benchmark_cross_dataset_comparable",
            "candidate_specific_cross_dataset_comparable",
            "read_level_absolute_cross_dataset_comparable",
            "pairwise_direction_stable",
            "comparability_warning",
        ],
        comparability_rows,
    )
    write_tsv_rows(
        output_dir / "layered_performance_matrix.tsv",
        [
            "dataset_id",
            "label",
            "sample",
            "platform",
            "mode",
            "layer",
            "method",
            "caller_family",
            "available",
            "truth_total",
            "calls_total",
            "tp",
            "fp",
            "fn",
            "precision",
            "recall",
            "f1",
            "source_path",
            "notes",
        ],
        layered_rows,
    )
    write_tsv_rows(
        output_dir / "candidate_pool_layer_matrix.tsv",
        [
            "dataset_id",
            "label",
            "sample",
            "platform",
            "mode",
            "layer",
            "subset",
            "rows_total",
            "rows_with_joined_features",
            "rows_analyzed",
            "coverage_fraction",
        ],
        candidate_pool_rows,
    )
    write_tsv_rows(
        output_dir / "feature_distribution_focus.tsv",
        [
            "dataset_id",
            "label",
            "sample",
            "platform",
            "mode",
            "feature",
            "truth_status",
            "n",
            "mean",
            "median",
            "p25",
            "p75",
            "p10",
            "p90",
            "source_path",
        ],
        feature_distribution_rows,
    )
    write_tsv_rows(
        output_dir / "selected_rule_matrix.tsv",
        [
            "dataset_id",
            "label",
            "sample",
            "platform",
            "mode",
            "rule_group",
            "rule_id",
            "feature_family",
            "gate_id",
            "trigger_count",
            "tp_rescued",
            "fp_reintroduced",
            "precision",
            "recall",
            "f1",
            "delta_f1_vs_baseline",
            "delta_f1_vs_gate",
            "fp_per_tp",
            "meets_safety",
            "source_path",
            "notes",
        ],
        selected_rule_rows,
    )
    write_tsv_rows(output_dir / "feature_overlap_focus.tsv", list(pd.DataFrame(overlap_rows).columns), overlap_rows)
    write_tsv_rows(
        output_dir / "snapshot_scope_bias_assessment.tsv",
        [
            "dataset_id",
            "label",
            "mode",
            "diagnostics_available",
            "snapshot_bam_scope",
            "snapshot_bam_path",
            "regions_profiled",
            "median_target_depth",
            "median_target_alt_fraction",
            "median_na_hp_fraction",
            "high_na_hp_regions",
            "haplotype_skewed_regions",
            "cross_dataset_absolute_comparable",
            "within_dataset_interpretation",
            "notes",
        ],
        snapshot_rows,
    )
    write_tsv_rows(
        output_dir / "phase2_gap_tracker.tsv",
        ["gap_id", "priority", "affected_scope", "blocking_conclusion", "required_action"],
        phase2_rows,
    )

    report_text = build_final_report(
        output_dir,
        model_rows,
        availability_rows,
        comparability_rows,
        layered_rows,
        feature_distribution_rows,
        selected_rule_rows,
        overlap_rows,
        snapshot_rows,
        phase2_rows,
    )
    report_path = output_dir / "final_integrated_report.md"
    report_path.write_text(report_text, encoding="utf-8")
    print(f"[build_hcc1395_cross_quadrant_matrix] Wrote outputs to {output_dir}")


if __name__ == "__main__":
    main()
