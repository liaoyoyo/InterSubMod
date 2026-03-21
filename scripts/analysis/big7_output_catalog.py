#!/usr/bin/env python3
"""Canonical big7 output metadata for InterSubMod."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable


BIG7_OUTPUT_ROOT = Path("/big7_disk/liaoyoyo2001/big7_disk_output")
CANONICAL_ROOT = BIG7_OUTPUT_ROOT / "canonical"
SYNTHESIS_ROOT = BIG7_OUTPUT_ROOT / "synthesis"
RUNBOOK_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod_big7_runbook")
RUNBOOK_MANIFEST_ROOT = RUNBOOK_ROOT / "manifests"

ARCHIVE_ROOTS: Dict[str, Path] = {
    "big8_output_archive": BIG7_OUTPUT_ROOT / "big8_output_archive",
    "bip8_output_archive": BIG7_OUTPUT_ROOT / "bip8_output_archive",
}


@dataclass(frozen=True)
class ModeSpec:
    canonical_mode: str
    paired_or_to: str
    caller_family: str
    caller_model: str
    baseline_method: str


MODE_SPECS: Dict[str, ModeSpec] = {
    "paired_full": ModeSpec("paired_full", "paired", "ClairS", "full", "LongPhase-S"),
    "paired_pileup": ModeSpec("paired_pileup", "paired", "ClairS", "pileup", "LongPhase-S"),
    "to_full": ModeSpec("to_full", "tumor_only", "ClairS-TO", "full", "LongPhase-TO"),
    "to_pileup": ModeSpec("to_pileup", "tumor_only", "ClairS-TO", "pileup", "LongPhase-TO"),
}


@dataclass(frozen=True)
class SampleSpec:
    sample: str
    platform: str
    methylation_status: str
    truth_set: str
    truth_total: int
    tumor_bam: str
    normal_bam: str
    truth_vcf: str
    truth_bed: str
    somatic_vcf_full: str
    somatic_vcf_pileup: str
    to_somatic_vcf_pileup: str = ""
    to_somatic_vcf_full: str = ""

    @property
    def expected_modes(self) -> tuple[str, ...]:
        return self.available_modes()

    def available_modes(self) -> tuple[str, ...]:
        modes = []
        if self.somatic_vcf_full:
            modes.append("paired_full")
        if self.somatic_vcf_pileup:
            modes.append("paired_pileup")
        if self.to_somatic_vcf_full:
            modes.append("to_full")
        if self.to_somatic_vcf_pileup:
            modes.append("to_pileup")
        return tuple(modes)

    def mode_unavailable_reason(self, canonical_mode: str) -> str:
        if canonical_mode == "paired_full":
            return "paired_full_source_unavailable"
        if canonical_mode == "paired_pileup":
            return "paired_pileup_source_unavailable"
        if canonical_mode == "to_full":
            return "tumor_only_full_model_unavailable" if self.to_somatic_vcf_pileup else "tumor_only_source_unavailable"
        if canonical_mode == "to_pileup":
            return "tumor_only_source_unavailable"
        return "mode_unavailable"


SAMPLE_SPECS: Dict[str, SampleSpec] = {
    "HCC1395": SampleSpec(
        sample="HCC1395",
        platform="ONT_5kHz",
        methylation_status="5mCG+5hmCG",
        truth_set="SEQC2_v1.2.1_HC_SNV",
        truth_total=39447,
        tumor_bam="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam",
        normal_bam="/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam",
        truth_vcf="/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz",
        truth_bed="/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed",
        somatic_vcf_full="/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/output.vcf.gz",
        somatic_vcf_pileup="/big8_disk/data/HCC1395/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf",
        to_somatic_vcf_pileup="/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz",
    ),
    "HCC1395_DORADO": SampleSpec(
        sample="HCC1395_DORADO",
        platform="ONT_Dorado",
        methylation_status="5mCG",
        truth_set="SEQC2_v1.2.1_HC_SNV",
        truth_total=39447,
        tumor_bam="/big8_disk/data/HCC1395/ONT_Dorado/HCC1395.bam",
        normal_bam="/big8_disk/data/HCC1395/ONT_Dorado/HCC1395BL.bam",
        truth_vcf="/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz",
        truth_bed="/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed",
        somatic_vcf_full="/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/output.vcf.gz",
        somatic_vcf_pileup="/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf",
        to_somatic_vcf_pileup="/big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_v0_3_0/snv.vcf.gz",
    ),
    "COLO829": SampleSpec(
        sample="COLO829",
        platform="ONT_PAO",
        methylation_status="5mCG+5hmCG",
        truth_set="NYGC",
        truth_total=41427,
        tumor_bam="/big8_disk/data/COLO829/ONT_PAO/PAO29420.bam",
        normal_bam="/big8_disk/data/COLO829/ONT_PAO/PAO33946.bam",
        truth_vcf="/big8_disk/data/COLO829/NYGC/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf.gz",
        truth_bed="",
        somatic_vcf_full="/big8_disk/data/COLO829/ONT_PAO/ClairS_v0_4_1/output.vcf.gz",
        somatic_vcf_pileup="/big8_disk/data/COLO829/ONT_PAO/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf",
        to_somatic_vcf_pileup="/big8_disk/data/COLO829/ONT_PAO/ClairS_TO_v0_3_0/snv.vcf.gz",
    ),
    "H1437": SampleSpec(
        sample="H1437",
        platform="ONT",
        methylation_status="5mCG",
        truth_set="orthogonal-tools",
        truth_total=90129,
        tumor_bam="/big8_disk/data/H1437/ONT/H1437.bam",
        normal_bam="/big8_disk/data/H1437/ONT/H1437BL.bam",
        truth_vcf="/big8_disk/data/H1437/orthogonal-tools-benchmark/H1437_orthogonal-tools-benchmark_somatic-only.vcf.gz",
        truth_bed="/big8_disk/data/H1437/orthogonal-tools-benchmark/H1437_orthogonal-tools-benchmark.bed",
        somatic_vcf_full="/big8_disk/data/H1437/ONT/ClairS_v0_4_1/output.vcf.gz",
        somatic_vcf_pileup="/big8_disk/data/H1437/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf",
        to_somatic_vcf_pileup="/big8_disk/data/H1437/ONT/ClairS_TO_v0_3_0/snv.vcf.gz",
    ),
    "H2009": SampleSpec(
        sample="H2009",
        platform="ONT",
        methylation_status="5mCG",
        truth_set="orthogonal-tools",
        truth_total=168529,
        tumor_bam="/big8_disk/data/H2009/ONT/H2009.bam",
        normal_bam="/big8_disk/data/H2009/ONT/H2009BL.bam",
        truth_vcf="/big8_disk/data/H2009/orthogonal-tools-benchmark/H2009_orthogonal-tools-benchmark_somatic-only.vcf.gz",
        truth_bed="/big8_disk/data/H2009/orthogonal-tools-benchmark/H2009_orthogonal-tools-benchmark.bed",
        somatic_vcf_full="/big8_disk/data/H2009/ONT/ClairS_v0_4_1/output.vcf.gz",
        somatic_vcf_pileup="/big8_disk/data/H2009/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf",
        to_somatic_vcf_pileup="/big8_disk/data/H2009/ONT/ClairS_TO_v0_3_0/snv.vcf.gz",
    ),
    "HCC1937": SampleSpec(
        sample="HCC1937",
        platform="ONT",
        methylation_status="5mCG",
        truth_set="orthogonal-tools",
        truth_total=60691,
        tumor_bam="/big8_disk/data/HCC1937/ONT/HCC1937.bam",
        normal_bam="/big8_disk/data/HCC1937/ONT/HCC1937BL.bam",
        truth_vcf="/big8_disk/data/HCC1937/orthogonal-tools-benchmark/HCC1937_orthogonal-tools-benchmark_somatic-only.vcf.gz",
        truth_bed="/big8_disk/data/HCC1937/orthogonal-tools-benchmark/HCC1937_orthogonal-tools-benchmark.bed",
        somatic_vcf_full="/big8_disk/data/HCC1937/ONT/ClairS_v0_4_1/output.vcf.gz",
        somatic_vcf_pileup="/big8_disk/data/HCC1937/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf",
        to_somatic_vcf_pileup="/big8_disk/data/HCC1937/ONT/ClairS_TO_v0_3_0/snv.vcf.gz",
    ),
    "HCC1954": SampleSpec(
        sample="HCC1954",
        platform="ONT",
        methylation_status="5mCG",
        truth_set="orthogonal-tools",
        truth_total=26567,
        tumor_bam="/big8_disk/data/HCC1954/ONT/HCC1954.bam",
        normal_bam="/big8_disk/data/HCC1954/ONT/HCC1954BL.bam",
        truth_vcf="/big8_disk/data/HCC1954/orthogonal-tools-benchmark/HCC1954_orthogonal-tools-benchmark_somatic-only.vcf.gz",
        truth_bed="/big8_disk/data/HCC1954/orthogonal-tools-benchmark/HCC1954_orthogonal-tools-benchmark.bed",
        somatic_vcf_full="/big8_disk/data/HCC1954/ONT/ClairS_v0_4_1/output.vcf.gz",
        somatic_vcf_pileup="/big8_disk/data/HCC1954/ONT/ClairS_v0_4_1/tmp/vcf_output/pileup_filter.vcf",
        to_somatic_vcf_pileup="/big8_disk/data/HCC1954/ONT/ClairS_TO_v0_3_0/snv.vcf.gz",
    ),
}


def legacy_mode_to_canonical(legacy_mode: str) -> str:
    if legacy_mode == "s-pure":
        return "paired_full"
    if legacy_mode == "s-pure-pileup":
        return "paired_pileup"
    if legacy_mode in {"to-pure", "to_pure"}:
        return "to_pileup"
    return legacy_mode


def sanitize_component(text: str) -> str:
    allowed = []
    for char in text:
        if char.isalnum():
            allowed.append(char)
        elif char in {"-", "_"}:
            allowed.append(char)
        else:
            allowed.append("_")
    return "".join(allowed).strip("_")


def canonical_run_id(sample: str, canonical_mode: str, caller_model: str, source_name: str) -> str:
    source = sanitize_component(source_name)
    sample_name = sanitize_component(sample)
    return f"{source}_{sample_name}_{canonical_mode}_{caller_model}"


def sample_specs() -> Iterable[SampleSpec]:
    return SAMPLE_SPECS.values()
