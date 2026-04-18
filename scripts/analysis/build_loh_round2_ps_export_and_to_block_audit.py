#!/usr/bin/env python3
"""
LOH Round 2:
1. paired PS export feasibility audit
2. TO-only FP PS-block clustering audit

This script intentionally separates two questions:
  - 5.2 schema/export: where paired PS is lost, and whether a read-level
    fallback is immediately feasible from tagged BAMs
  - 5.1 research: whether TO-only FP loci cluster within PS blocks and whether
    those blocks look different from TP-heavy blocks
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import random
import shutil
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns


ROUND1_WS = Path(
    os.environ.get(
        "LOH_ROUND1_WS",
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit",
    )
)
OUT_DIR = Path(
    os.environ.get(
        "LOH_ROUND2_PS_OUT_DIR",
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_loh_round2_ps_export_and_to_block_audit",
    )
)
FIG_DIR = OUT_DIR / "figures"
REPORT_ASSET_DIR = Path(
    os.environ.get(
        "LOH_ROUND2_PS_REPORT_ASSET_DIR",
        "/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260330_loh_round2_ps_export_and_to_block_audit",
    )
)

ROUND1_ALL_ROWS = ROUND1_WS / "all_region_rows.tsv.gz"
GENERATED_AT = datetime.now().strftime("%Y-%m-%d %H:%M:%S %Z")
PILOT_RANDOM_SEED = 20260330
PAIRED_PILOT_PER_SAMPLE_PER_TRUTH = 10
FOCUS_SAMPLES = ["HCC1954", "COLO829", "H1437", "HCC1395_DORADO"]
SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
SAMPLE_LABELS = {
    "HCC1395": "HCC1395 5kHz",
    "HCC1395_DORADO": "HCC1395 DORADO",
    "COLO829": "COLO829",
    "H1437": "H1437",
    "H2009": "H2009",
    "HCC1937": "HCC1937",
    "HCC1954": "HCC1954",
}

TRUTH_COLORS = {"TP": "#2563eb", "FP": "#dc2626"}
STATE_ORDER = [
    "HP1-1-like",
    "HP2-1-like",
    "HP3-like",
    "LOH hap1-side",
    "LOH hap2-side",
    "Other/Unknown",
]
STATE_COLORS = {
    "HP1-1-like": "#2563eb",
    "HP2-1-like": "#0f766e",
    "HP3-like": "#f59e0b",
    "LOH hap1-side": "#7c3aed",
    "LOH hap2-side": "#db2777",
    "Other/Unknown": "#6b7280",
}


@dataclass(frozen=True)
class PairedArtifacts:
    sample: str
    base_dir: Path
    tagged_bam: Path
    germline_phased_vcf: Path
    somatic_pass_vcf: Path
    filtered_tp_vcf: Path
    filtered_fp_vcf: Path


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def normalize_missing(value) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "none"}:
        return ""
    return text


def status_with_default(value) -> str:
    text = normalize_missing(value)
    return text if text in {"TP", "FP"} else ""


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


def collapse_to_pseudo_state(gt: str, gt2: str, gt3: str) -> str:
    def canon(value: str) -> str:
        text = normalize_missing(value)
        mapping = {
            "0/0": "0|0",
            "0/1": "0|1",
            "1/0": "1|0",
            "1/1": "1|1",
            "0/.": "0|.",
            "./0": ".|0",
            "1/.": "1|.",
            "./1": ".|1",
        }
        return mapping.get(text, text)

    gt = canon(gt)
    gt2 = canon(gt2)
    gt3 = canon(gt3)

    if gt == "0|0":
        if gt2 == "1|." and gt3 in {"", "./."}:
            return "HP1-1-like"
        if gt2 == ".|1" and gt3 in {"", "./."}:
            return "HP2-1-like"
        if gt2 == "1|1" and gt3 in {"", "./."}:
            return "HP3-like"
        return "Other/Unknown"

    if gt == ".|0":
        if gt2 in {"0|.", "1|."} and gt3 in {"0|.", "1|."}:
            if gt2 == "1|." and gt3 == "1|.":
                return "HP3-like"
            return "LOH hap1-side"
        return "LOH hap1-side"

    if gt == "0|.":
        if gt2 in {".|0", ".|1"} and gt3 in {".|0", ".|1"}:
            if gt2 == ".|1" and gt3 == ".|1":
                return "HP3-like"
            return "LOH hap2-side"
        return "LOH hap2-side"

    return "Other/Unknown"


def format_pct(value: float) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return f"{value * 100:.1f}%"


def load_round1_all_rows() -> pd.DataFrame:
    print("[1/7] Loading round1 all_region_rows ...")
    usecols = [
        "sample",
        "platform",
        "mode",
        "truth_label",
        "variant_key",
        "Chr",
        "Pos",
        "Ref",
        "Alt",
        "source_vcf_file",
        "source_tagged_bam",
        "effective_hp_reads",
        "core_loh_like",
        "verification_class",
        "loh_subtype",
        "hp_ratio_core",
        "hp_assign_rate",
        "hp0_ratio",
        "hp3_ratio",
        "allele_delta",
        "pairwise_median_dist",
        "quality_score",
        "caller_af",
        "caller_filter",
        "to_phase_ps",
        "to_phase_gt",
        "to_phase_gt2",
        "to_phase_ps2",
        "to_phase_gt3",
        "to_phase_ps3",
        "to_loh_bed_hit",
        "phase_block_status",
    ]
    df = pd.read_csv(ROUND1_ALL_ROWS, sep="\t", usecols=usecols, low_memory=False)
    df = df[df["truth_label"].isin(["TP", "FP"])].copy()
    df["Pos"] = pd.to_numeric(df["Pos"], errors="coerce")
    numeric_cols = [
        "effective_hp_reads",
        "hp_ratio_core",
        "hp_assign_rate",
        "hp0_ratio",
        "hp3_ratio",
        "allele_delta",
        "pairwise_median_dist",
        "quality_score",
        "caller_af",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in [
        "to_phase_ps",
        "to_phase_gt",
        "to_phase_gt2",
        "to_phase_ps2",
        "to_phase_gt3",
        "to_phase_ps3",
        "phase_block_status",
        "verification_class",
        "loh_subtype",
        "caller_filter",
    ]:
        df[col] = df[col].map(normalize_missing)
    print(
        f"  rows={len(df):,}, paired={df['mode'].eq('paired').sum():,}, to={df['mode'].eq('to').sum():,}"
    )
    return df


def derive_paired_artifacts(all_df: pd.DataFrame) -> List[PairedArtifacts]:
    artifacts: List[PairedArtifacts] = []
    paired_rows = (
        all_df[all_df["mode"] == "paired"][["sample", "source_vcf_file", "source_tagged_bam"]]
        .drop_duplicates()
        .sort_values(["sample", "source_vcf_file"])
    )
    grouped = paired_rows.groupby("sample", sort=False)
    for sample, group in grouped:
        row = group.iloc[0]
        base_dir = Path(row["source_vcf_file"]).parent
        artifacts.append(
            PairedArtifacts(
                sample=sample,
                base_dir=base_dir,
                tagged_bam=Path(row["source_tagged_bam"]),
                germline_phased_vcf=base_dir / "germline_phased_merged.vcf.gz",
                somatic_pass_vcf=base_dir / "somatic_pass.vcf.gz",
                filtered_tp_vcf=base_dir / "filtered_snv_tp.vcf.gz",
                filtered_fp_vcf=base_dir / "filtered_snv_fp.vcf.gz",
            )
        )
    return artifacts


def audit_vcf_ps(path: Path, sample: str, artifact_type: str) -> Dict:
    total_records = 0
    ps_nonempty_records = 0
    ps_header_present = False
    exists = path.exists()
    if exists:
        with open_text(path) as handle:
            for line in handle:
                if line.startswith("##FORMAT=<ID=PS,"):
                    ps_header_present = True
                    continue
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                total_records += 1
                if len(fields) < 10:
                    continue
                fmt = fields[8].split(":")
                sample_values = fields[9].split(":")
                value_map = dict(zip(fmt, sample_values))
                if normalize_missing(value_map.get("PS")) not in {"", "."}:
                    ps_nonempty_records += 1
    return {
        "sample": sample,
        "artifact_type": artifact_type,
        "path": str(path),
        "exists": exists,
        "ps_header_present": ps_header_present,
        "total_records": total_records,
        "ps_nonempty_records": ps_nonempty_records,
        "ps_nonempty_fraction": (
            ps_nonempty_records / total_records if total_records else math.nan
        ),
    }


def build_paired_variant_ps_audit(all_df: pd.DataFrame) -> pd.DataFrame:
    print("[2/7] Auditing paired variant-level PS propagation ...")
    rows: List[Dict] = []
    for bundle in derive_paired_artifacts(all_df):
        rows.append(audit_vcf_ps(bundle.germline_phased_vcf, bundle.sample, "germline_phased_merged"))
        rows.append(audit_vcf_ps(bundle.somatic_pass_vcf, bundle.sample, "somatic_pass"))
        rows.append(audit_vcf_ps(bundle.filtered_tp_vcf, bundle.sample, "filtered_snv_tp"))
        rows.append(audit_vcf_ps(bundle.filtered_fp_vcf, bundle.sample, "filtered_snv_fp"))
    df = pd.DataFrame(rows)
    df["sample"] = pd.Categorical(df["sample"], categories=SAMPLE_ORDER, ordered=True)
    df = df.sort_values(["sample", "artifact_type"]).reset_index(drop=True)
    df.to_csv(OUT_DIR / "paired_variant_ps_audit.tsv", sep="\t", index=False)

    summary = (
        df.groupby("artifact_type", dropna=False)
        .agg(
            samples=("sample", "count"),
            total_records=("total_records", "sum"),
            ps_nonempty_records=("ps_nonempty_records", "sum"),
            mean_ps_nonempty_fraction=("ps_nonempty_fraction", "mean"),
            samples_with_ps_header=("ps_header_present", "sum"),
        )
        .reset_index()
    )
    summary.to_csv(OUT_DIR / "paired_variant_ps_audit_summary.tsv", sep="\t", index=False)
    print(f"  wrote {len(df)} audit rows")
    return df


def choose_paired_pilot_sites(all_df: pd.DataFrame) -> pd.DataFrame:
    rng = random.Random(PILOT_RANDOM_SEED)
    paired_df = all_df[all_df["mode"] == "paired"].copy()
    paired_df = paired_df[
        [
            "sample",
            "truth_label",
            "variant_key",
            "Chr",
            "Pos",
            "Ref",
            "Alt",
            "source_tagged_bam",
            "effective_hp_reads",
            "quality_score",
        ]
    ].drop_duplicates(subset=["sample", "truth_label", "variant_key"])

    picked_frames: List[pd.DataFrame] = []
    for (sample, truth_label), group in paired_df.groupby(["sample", "truth_label"], sort=False):
        records = group.sort_values(["Chr", "Pos", "variant_key"]).to_dict("records")
        pick_n = min(PAIRED_PILOT_PER_SAMPLE_PER_TRUTH, len(records))
        selected = rng.sample(records, pick_n) if len(records) > pick_n else records
        picked_frames.append(pd.DataFrame(selected))

    pilot_df = pd.concat(picked_frames, ignore_index=True)
    pilot_df["sample"] = pd.Categorical(pilot_df["sample"], categories=SAMPLE_ORDER, ordered=True)
    pilot_df = pilot_df.sort_values(["sample", "truth_label", "Chr", "Pos"]).reset_index(drop=True)
    pilot_df.to_csv(OUT_DIR / "paired_read_ps_pilot_sites.tsv", sep="\t", index=False)
    return pilot_df


def pileup_read_ps_fraction(bam: pysam.AlignmentFile, chrom: str, pos_1based: int) -> Dict[str, float]:
    per_read: Dict[str, Dict] = {}
    for column in bam.pileup(
        chrom,
        pos_1based - 1,
        pos_1based,
        truncate=True,
        stepper="samtools",
        max_depth=100000,
    ):
        if column.reference_pos != pos_1based - 1:
            continue
        for pileup_read in column.pileups:
            aln = pileup_read.alignment
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.is_duplicate:
                continue
            if aln.query_name in per_read:
                continue
            ps_value = normalize_missing(aln.get_tag("PS") if aln.has_tag("PS") else "")
            hp_value = normalize_missing(aln.get_tag("HP") if aln.has_tag("HP") else "")
            per_read[aln.query_name] = {
                "ps": ps_value,
                "hp": hp_value,
            }

    total_reads = len(per_read)
    ps_reads = sum(1 for row in per_read.values() if row["ps"] not in {"", "."})
    hp_reads = sum(1 for row in per_read.values() if row["hp"] not in {"", "."})
    unique_ps = sorted({row["ps"] for row in per_read.values() if row["ps"] not in {"", "."}})
    return {
        "pileup_total_reads": total_reads,
        "pileup_ps_reads": ps_reads,
        "pileup_hp_reads": hp_reads,
        "pileup_ps_fraction": ps_reads / total_reads if total_reads else math.nan,
        "pileup_hp_fraction": hp_reads / total_reads if total_reads else math.nan,
        "pileup_unique_ps": len(unique_ps),
    }


def build_paired_read_ps_pilot(pilot_sites: pd.DataFrame) -> pd.DataFrame:
    print("[3/7] Running paired read-level PS fallback pilot ...")
    rows: List[Dict] = []
    bam_cache: Dict[str, pysam.AlignmentFile] = {}
    try:
        for row in pilot_sites.itertuples(index=False):
            bam_path = str(row.source_tagged_bam)
            if bam_path not in bam_cache:
                bam_cache[bam_path] = pysam.AlignmentFile(bam_path, "rb")
            bam = bam_cache[bam_path]
            stats = pileup_read_ps_fraction(bam, row.Chr, int(row.Pos))
            rows.append(
                {
                    "sample": row.sample,
                    "truth_label": row.truth_label,
                    "variant_key": row.variant_key,
                    "Chr": row.Chr,
                    "Pos": int(row.Pos),
                    "Ref": row.Ref,
                    "Alt": row.Alt,
                    "source_tagged_bam": bam_path,
                    **stats,
                }
            )
    finally:
        for bam in bam_cache.values():
            bam.close()

    pilot_df = pd.DataFrame(rows)
    pilot_df["sample"] = pd.Categorical(pilot_df["sample"], categories=SAMPLE_ORDER, ordered=True)
    pilot_df = pilot_df.sort_values(["sample", "truth_label", "Chr", "Pos"]).reset_index(drop=True)
    pilot_df.to_csv(OUT_DIR / "paired_read_ps_pilot.tsv", sep="\t", index=False)

    summary = (
        pilot_df.groupby(["sample", "truth_label"], dropna=False)
        .agg(
            pilot_sites=("variant_key", "count"),
            median_pileup_total_reads=("pileup_total_reads", "median"),
            median_pileup_ps_fraction=("pileup_ps_fraction", "median"),
            loci_with_any_ps=("pileup_ps_reads", lambda s: int((s > 0).sum())),
            loci_with_multiple_ps=("pileup_unique_ps", lambda s: int((s >= 2).sum())),
        )
        .reset_index()
    )
    summary["loci_with_any_ps_fraction"] = summary["loci_with_any_ps"] / summary["pilot_sites"]
    summary["loci_with_multiple_ps_fraction"] = summary["loci_with_multiple_ps"] / summary["pilot_sites"]
    summary.to_csv(OUT_DIR / "paired_read_ps_pilot_summary.tsv", sep="\t", index=False)
    return pilot_df


def build_status_wide(all_df: pd.DataFrame) -> pd.DataFrame:
    status_df = (
        all_df[["sample", "variant_key", "mode", "truth_label"]]
        .drop_duplicates(subset=["sample", "variant_key", "mode"])
        .copy()
    )
    wide = (
        status_df.pivot(index=["sample", "variant_key"], columns="mode", values="truth_label")
        .reset_index()
        .rename(columns={"paired": "paired_status", "to": "to_status"})
    )
    for col in ["paired_status", "to_status"]:
        if col not in wide.columns:
            wide[col] = ""
        wide[col] = wide[col].map(status_with_default)
    wide["paired_present"] = wide["paired_status"] != ""
    wide["to_present"] = wide["to_status"] != ""
    wide["locus_overlap_class"] = np.select(
        [wide["paired_present"] & wide["to_present"], wide["paired_present"], wide["to_present"]],
        ["both", "paired_only", "to_only"],
        default="none",
    )
    wide["paired_to_concordance"] = [
        compare_status(paired_status, to_status, overlap)
        for paired_status, to_status, overlap in zip(
            wide["paired_status"].tolist(),
            wide["to_status"].tolist(),
            wide["locus_overlap_class"].tolist(),
        )
    ]
    return wide


def build_to_block_audit(all_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    print("[4/7] Building TO-only FP PS-block audit ...")
    status_wide = build_status_wide(all_df)
    to_df = all_df[all_df["mode"] == "to"].copy()
    to_df = to_df.merge(status_wide, on=["sample", "variant_key"], how="left")
    to_df["phase_block_id"] = to_df["to_phase_ps"].map(normalize_missing)
    to_df["has_ps"] = to_df["phase_block_id"] != ""
    to_df["is_to_only_fp"] = to_df["paired_to_concordance"] == "to_only_fp"
    to_df["is_to_tp"] = to_df["truth_label"] == "TP"
    to_df["is_to_fp"] = to_df["truth_label"] == "FP"
    to_df["pseudo_phase_state"] = [
        collapse_to_pseudo_state(gt, gt2, gt3)
        for gt, gt2, gt3 in zip(
            to_df["to_phase_gt"].tolist(),
            to_df["to_phase_gt2"].tolist(),
            to_df["to_phase_gt3"].tolist(),
        )
    ]
    to_df.to_csv(OUT_DIR / "to_locus_status_table.tsv", sep="\t", index=False)

    ps_df = to_df[to_df["has_ps"]].copy()
    block_rows: List[Dict] = []
    for (sample, phase_block_id), group in ps_df.groupby(["sample", "phase_block_id"], sort=False):
        fp_only = group[group["is_to_only_fp"]]
        pseudo_counts = Counter(fp_only["pseudo_phase_state"].tolist())
        dominant_state = pseudo_counts.most_common(1)[0][0] if pseudo_counts else ""
        block_rows.append(
            {
                "sample": sample,
                "phase_block_id": phase_block_id,
                "block_start": int(group["Pos"].min()),
                "block_end": int(group["Pos"].max()),
                "block_total_loci": len(group),
                "tp_n": int(group["is_to_tp"].sum()),
                "fp_n": int(group["is_to_fp"].sum()),
                "to_only_fp_n": int(group["is_to_only_fp"].sum()),
                "to_only_fp_fraction_within_block": group["is_to_only_fp"].mean(),
                "core_loh_like_to_only_fp_n": int(fp_only["core_loh_like"].sum()),
                "core_loh_like_to_only_fp_fraction": (
                    fp_only["core_loh_like"].mean() if len(fp_only) else math.nan
                ),
                "mean_effective_hp_reads": group["effective_hp_reads"].mean(),
                "mean_hp_ratio_core": group["hp_ratio_core"].mean(),
                "mean_quality_score": group["quality_score"].mean(),
                "dominant_to_only_fp_pseudo_state": dominant_state,
            }
        )

    block_df = pd.DataFrame(block_rows)
    block_df["sample"] = pd.Categorical(block_df["sample"], categories=SAMPLE_ORDER, ordered=True)
    block_df = block_df.sort_values(
        ["sample", "to_only_fp_n", "block_total_loci", "phase_block_id"],
        ascending=[True, False, False, True],
    ).reset_index(drop=True)
    block_df.to_csv(OUT_DIR / "to_ps_block_summary.tsv", sep="\t", index=False)

    state_summary = (
        ps_df[ps_df["is_to_only_fp"]]
        .groupby(["sample", "pseudo_phase_state"], dropna=False)
        .size()
        .reset_index(name="count")
    )
    state_summary["sample"] = pd.Categorical(state_summary["sample"], categories=SAMPLE_ORDER, ordered=True)
    state_summary = state_summary.sort_values(["sample", "pseudo_phase_state"]).reset_index(drop=True)
    state_summary.to_csv(OUT_DIR / "to_only_fp_pseudo_state_summary.tsv", sep="\t", index=False)

    concentration_rows: List[Dict] = []
    for sample, sample_df in ps_df.groupby("sample", sort=False):
        for label_name, flag_col in [("to_only_fp", "is_to_only_fp"), ("to_tp", "is_to_tp")]:
            counts = (
                sample_df.groupby("phase_block_id")[flag_col].sum().sort_values(ascending=False)
            )
            occupied = counts[counts > 0]
            total = int(occupied.sum())
            top1 = int(occupied.iloc[:1].sum()) if len(occupied) else 0
            top3 = int(occupied.iloc[:3].sum()) if len(occupied) else 0
            top5 = int(occupied.iloc[:5].sum()) if len(occupied) else 0
            ge2 = int(occupied[occupied >= 2].sum()) if len(occupied) else 0
            concentration_rows.append(
                {
                    "sample": sample,
                    "label_group": label_name,
                    "total_loci": total,
                    "occupied_blocks": int(len(occupied)),
                    "top1_share": top1 / total if total else math.nan,
                    "top3_share": top3 / total if total else math.nan,
                    "top5_share": top5 / total if total else math.nan,
                    "ge2_block_share": ge2 / total if total else math.nan,
                    "mean_loci_per_occupied_block": occupied.mean() if len(occupied) else math.nan,
                    "median_loci_per_occupied_block": occupied.median() if len(occupied) else math.nan,
                }
            )
    concentration_df = pd.DataFrame(concentration_rows)
    concentration_df["sample"] = pd.Categorical(
        concentration_df["sample"], categories=SAMPLE_ORDER, ordered=True
    )
    concentration_df = concentration_df.sort_values(["sample", "label_group"]).reset_index(drop=True)
    concentration_df.to_csv(OUT_DIR / "to_ps_block_concentration_summary.tsv", sep="\t", index=False)

    return to_df, block_df, concentration_df


def make_paired_variant_ps_heatmap(df: pd.DataFrame) -> Path:
    artifact_order = [
        "germline_phased_merged",
        "somatic_pass",
        "filtered_snv_tp",
        "filtered_snv_fp",
    ]
    plot_df = df.copy()
    plot_df["sample"] = plot_df["sample"].astype(str)
    pivot = (
        plot_df.pivot(index="sample", columns="artifact_type", values="ps_nonempty_fraction")
        .reindex(SAMPLE_ORDER)
        .reindex(columns=artifact_order)
    )
    annot = []
    for sample in pivot.index:
        row = []
        for artifact in pivot.columns:
            sub = plot_df[(plot_df["sample"] == sample) & (plot_df["artifact_type"] == artifact)]
            if sub.empty:
                row.append("NA")
                continue
            rec = sub.iloc[0]
            row.append(f"{rec['ps_nonempty_records']:,}/{rec['total_records']:,}\n{format_pct(rec['ps_nonempty_fraction'])}")
        annot.append(row)
    annot = np.array(annot)

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(
        pivot,
        annot=annot,
        fmt="",
        cmap="YlGnBu",
        vmin=0,
        vmax=1,
        cbar_kws={"label": "records with non-empty PS"},
        ax=ax,
    )
    ax.set_title("Paired variant-level PS audit across upstream and split VCF artifacts")
    ax.set_xlabel("artifact")
    ax.set_ylabel("sample")
    ax.set_yticklabels([SAMPLE_LABELS.get(x, x) for x in pivot.index], rotation=0)
    fig.tight_layout()
    out_path = FIG_DIR / "fig01_paired_variant_ps_audit_heatmap.png"
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out_path


def make_paired_read_ps_pilot_plot(df: pd.DataFrame) -> Path:
    plot_df = df.copy()
    plot_df["sample"] = plot_df["sample"].astype(str)
    plot_df["sample_label"] = plot_df["sample"].map(lambda x: SAMPLE_LABELS.get(x, x))

    fig, ax = plt.subplots(figsize=(14, 6))
    sns.boxplot(
        data=plot_df,
        x="sample_label",
        y="pileup_ps_fraction",
        hue="truth_label",
        palette=TRUTH_COLORS,
        showfliers=False,
        ax=ax,
    )
    sns.stripplot(
        data=plot_df,
        x="sample_label",
        y="pileup_ps_fraction",
        hue="truth_label",
        dodge=True,
        palette=TRUTH_COLORS,
        alpha=0.45,
        size=4,
        ax=ax,
    )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], title="truth")
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("sample")
    ax.set_ylabel("fraction of pileup reads carrying PS tag")
    ax.set_title("Paired tagged BAM read-level PS fallback pilot")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    out_path = FIG_DIR / "fig02_paired_read_ps_fallback_pilot.png"
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out_path


def make_to_block_concentration_plot(ps_df: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    axes = axes.flatten()
    for ax, sample in zip(axes, FOCUS_SAMPLES):
        sample_df = ps_df[ps_df["sample"] == sample].copy()
        for flag_col, label, color in [
            ("is_to_only_fp", "TO-only FP", "#dc2626"),
            ("is_to_tp", "TO TP", "#2563eb"),
        ]:
            counts = sample_df.groupby("phase_block_id")[flag_col].sum().sort_values(ascending=False)
            counts = counts[counts > 0]
            if len(counts) == 0:
                continue
            cumulative = counts.cumsum() / counts.sum()
            ax.plot(
                np.arange(1, len(cumulative) + 1),
                cumulative.values,
                marker="o",
                linewidth=2,
                markersize=3,
                label=label,
                color=color,
            )
        ax.set_title(SAMPLE_LABELS.get(sample, sample))
        ax.set_xlabel("top-k PS blocks")
        ax.set_ylabel("cumulative share of loci")
        ax.grid(True, alpha=0.2)
        ax.legend()
    fig.suptitle("TO-only FP vs TP PS-block concentration (focus samples)", y=1.02, fontsize=13)
    fig.tight_layout()
    out_path = FIG_DIR / "fig03_to_ps_block_concentration_curves.png"
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out_path


def make_to_block_scatter(block_df: pd.DataFrame) -> Path:
    focus_df = block_df[block_df["sample"].isin(FOCUS_SAMPLES)].copy()
    focus_df["sample"] = pd.Categorical(focus_df["sample"], categories=FOCUS_SAMPLES, ordered=True)
    focus_df = focus_df[focus_df["to_only_fp_n"] > 0].copy()

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    axes = axes.flatten()
    vmax = max(1.0, focus_df["core_loh_like_to_only_fp_fraction"].dropna().max()) if not focus_df.empty else 1.0
    sc = None

    for ax, sample in zip(axes, FOCUS_SAMPLES):
        sub = focus_df[focus_df["sample"] == sample].copy()
        if sub.empty:
            ax.set_title(SAMPLE_LABELS.get(sample, sample))
            ax.text(0.5, 0.5, "No TO-only FP blocks with PS", ha="center", va="center")
            continue
        sizes = 30 + 20 * np.sqrt(sub["block_total_loci"].clip(lower=1))
        sc = ax.scatter(
            sub["tp_n"],
            sub["to_only_fp_n"],
            s=sizes,
            c=sub["core_loh_like_to_only_fp_fraction"],
            cmap="YlOrRd",
            vmin=0,
            vmax=vmax,
            alpha=0.8,
            edgecolors="black",
            linewidths=0.4,
        )
        lim = max(sub["tp_n"].max(), sub["to_only_fp_n"].max(), 1)
        ax.plot([0, lim], [0, lim], linestyle="--", color="#6b7280", linewidth=1)
        ax.set_title(SAMPLE_LABELS.get(sample, sample))
        ax.set_xlabel("TO TP count within PS block")
        ax.set_ylabel("TO-only FP count within PS block")
        ax.grid(True, alpha=0.2)
    if sc is not None:
        cbar = fig.colorbar(sc, ax=axes.tolist(), fraction=0.025, pad=0.02)
        cbar.set_label("LOH-like fraction among TO-only FP in block")
    fig.suptitle("PS blocks carrying TO-only FP: TP mix vs FP load", y=1.02, fontsize=13)
    fig.tight_layout()
    out_path = FIG_DIR / "fig04_to_only_fp_block_scatter.png"
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out_path


def make_to_pseudo_state_plot(to_df: pd.DataFrame) -> Path:
    state_df = to_df[(to_df["is_to_only_fp"]) & (to_df["has_ps"])].copy()
    summary = (
        state_df.groupby(["sample", "pseudo_phase_state"]).size().reset_index(name="count")
    )
    pivot = (
        summary.pivot(index="sample", columns="pseudo_phase_state", values="count")
        .fillna(0)
        .reindex(SAMPLE_ORDER)
    )
    for state in STATE_ORDER:
        if state not in pivot.columns:
            pivot[state] = 0
    pivot = pivot[STATE_ORDER]
    fractions = pivot.div(pivot.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)

    fig, ax = plt.subplots(figsize=(12, 6))
    bottom = np.zeros(len(fractions))
    for state in STATE_ORDER:
        values = fractions[state].values
        ax.bar(
            range(len(fractions)),
            values,
            bottom=bottom,
            color=STATE_COLORS[state],
            label=state,
        )
        bottom += values
    ax.set_xticks(range(len(fractions)))
    ax.set_xticklabels([SAMPLE_LABELS.get(x, x) for x in fractions.index], rotation=25, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction within TO-only FP loci carrying PS")
    ax.set_xlabel("sample")
    ax.set_title("TO-only FP pseudo-state composition from GT/GT2/GT3")
    ax.legend(ncol=3, fontsize=9)
    fig.tight_layout()
    out_path = FIG_DIR / "fig05_to_only_fp_pseudo_state_composition.png"
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out_path


def write_decision_ledger(
    paired_variant_df: pd.DataFrame,
    paired_pilot_df: pd.DataFrame,
    to_df: pd.DataFrame,
    concentration_df: pd.DataFrame,
) -> pd.DataFrame:
    paired_split_nonempty = paired_variant_df[
        paired_variant_df["artifact_type"].isin(["somatic_pass", "filtered_snv_tp", "filtered_snv_fp"])
    ]["ps_nonempty_records"].sum()
    paired_germline_nonempty = paired_variant_df[
        paired_variant_df["artifact_type"] == "germline_phased_merged"
    ]["ps_nonempty_records"].sum()
    paired_pilot_summary = (
        paired_pilot_df.groupby("truth_label")["pileup_ps_fraction"].median().to_dict()
        if not paired_pilot_df.empty
        else {}
    )
    to_only_fp_total = int(to_df["is_to_only_fp"].sum())
    to_only_fp_with_ps = int((to_df["is_to_only_fp"] & to_df["has_ps"]).sum())
    focus_conc = concentration_df[
        concentration_df["sample"].isin(FOCUS_SAMPLES) & concentration_df["label_group"].eq("to_only_fp")
    ].copy()
    rows = [
        {
            "decision_id": "D1",
            "topic": "paired_variant_level_ps",
            "decision": "Current paired Round1 artifacts do not provide usable variant-level PS for somatic TP/FP VCFs.",
            "evidence": f"germline_ps_nonempty={paired_germline_nonempty:,}; paired_split_ps_nonempty={paired_split_nonempty:,}",
            "recommended_action": "Treat paired variant-level PS export as a tool/pipeline change request, not a missing merge bug.",
        },
        {
            "decision_id": "D2",
            "topic": "paired_read_level_ps_fallback",
            "decision": "Read-level PS fallback is feasible now from tagged BAM and can be added as an interim summary layer.",
            "evidence": (
                f"pilot_median_ps_fraction_TP={paired_pilot_summary.get('TP', math.nan):.3f}; "
                f"pilot_median_ps_fraction_FP={paired_pilot_summary.get('FP', math.nan):.3f}"
            ),
            "recommended_action": "Add variant_key -> read-level PS summary before attempting full paired same-locus block linkage.",
        },
        {
            "decision_id": "D3",
            "topic": "to_only_fp_block_unit",
            "decision": "TO-only FP should be tracked at the PS-block level, not only as isolated loci.",
            "evidence": f"to_only_fp_with_ps={to_only_fp_with_ps:,}/{to_only_fp_total:,}",
            "recommended_action": "Promote PS-block linkage audit into the next TO mainline diagnostic.",
        },
        {
            "decision_id": "D4",
            "topic": "to_focus_samples",
            "decision": "Focus samples remain HCC1954, COLO829, H1437, and HCC1395_DORADO because they show interpretable PS-block structure.",
            "evidence": "; ".join(
                f"{row.sample}: top5_share={row.top5_share:.3f}"
                for row in focus_conc.itertuples()
            ),
            "recommended_action": "Use these four samples for block-level follow-up with regional snapshots and pseudo-state linkage.",
        },
    ]
    ledger_df = pd.DataFrame(rows)
    ledger_df.to_csv(OUT_DIR / "decision_ledger.tsv", sep="\t", index=False)
    return ledger_df


def copy_figures_to_report_assets(fig_paths: Iterable[Path]) -> List[str]:
    ensure_dir(REPORT_ASSET_DIR)
    rel_paths: List[str] = []
    for fig_path in fig_paths:
        target = REPORT_ASSET_DIR / fig_path.name
        shutil.copy2(fig_path, target)
        rel_paths.append(f"assets/{REPORT_ASSET_DIR.name}/{fig_path.name}")
    return rel_paths


def write_round_context(fig_paths: List[Path]) -> None:
    payload = {
        "generated_at": GENERATED_AT,
        "round1_workspace": str(ROUND1_WS),
        "round1_all_region_rows": str(ROUND1_ALL_ROWS),
        "output_dir": str(OUT_DIR),
        "report_asset_dir": str(REPORT_ASSET_DIR),
        "paired_pilot_random_seed": PILOT_RANDOM_SEED,
        "paired_pilot_per_sample_per_truth": PAIRED_PILOT_PER_SAMPLE_PER_TRUTH,
        "focus_samples": FOCUS_SAMPLES,
        "figures": [str(path) for path in fig_paths],
    }
    with open(OUT_DIR / "round_context.json", "w") as handle:
        json.dump(payload, handle, indent=2)


def main() -> None:
    ensure_dir(OUT_DIR)
    ensure_dir(FIG_DIR)

    all_df = load_round1_all_rows()
    paired_variant_df = build_paired_variant_ps_audit(all_df)
    pilot_sites = choose_paired_pilot_sites(all_df)
    paired_pilot_df = build_paired_read_ps_pilot(pilot_sites)
    to_df, block_df, concentration_df = build_to_block_audit(all_df)

    print("[5/7] Generating figures ...")
    fig_paths = [
        make_paired_variant_ps_heatmap(paired_variant_df),
        make_paired_read_ps_pilot_plot(paired_pilot_df),
        make_to_block_concentration_plot(to_df[to_df["has_ps"]].copy()),
        make_to_block_scatter(block_df),
        make_to_pseudo_state_plot(to_df),
    ]
    copy_figures_to_report_assets(fig_paths)

    print("[6/7] Writing decision ledger and run context ...")
    write_decision_ledger(paired_variant_df, paired_pilot_df, to_df, concentration_df)
    write_round_context(fig_paths)

    print("[7/7] Done")
    print(f"  output_dir={OUT_DIR}")
    for fig_path in fig_paths:
        print(f"  figure={fig_path}")


if __name__ == "__main__":
    main()
