#!/usr/bin/env python3
"""Prepare latest LongPhase-S PASS truth-labeled FP/TP inputs without mutating canonical data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pysam


AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
DATASETS = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(command: list[str], command_log: list[list[str]]) -> None:
    command_log.append(command)
    subprocess.run(command, check=True)


def variant_keys(path: Path) -> set[tuple[str, int, str, str]]:
    keys: set[tuple[str, int, str, str]] = set()
    with pysam.VariantFile(str(path)) as variants:
        for record in variants:
            if len(record.ref) != 1 or not record.alts or len(record.alts) != 1 or len(record.alts[0]) != 1:
                continue
            keys.add((record.contig, int(record.pos), record.ref.upper(), record.alts[0].upper()))
    return keys


def count_records(path: Path) -> int:
    with pysam.VariantFile(str(path)) as variants:
        return sum(1 for _ in variants)


def count_matrices(root: Path) -> int:
    return sum(1 for _ in root.glob("**/distance/BERNOULLI/matrix.csv"))


def merge_windows(keys: Iterable[tuple[str, int, str, str]], flank: int = 5000) -> list[tuple[str, int, int]]:
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for chrom, pos, _ref, _alt in keys:
        start0 = max(0, pos - 1 - flank)
        end0 = pos + flank
        by_chrom[chrom].append((start0, end0))
    merged: list[tuple[str, int, int]] = []
    rank = {chrom: index for index, chrom in enumerate(AUTOSOMES)}
    for chrom in sorted(by_chrom, key=lambda value: rank.get(value, 10**9)):
        for start0, end0 in sorted(by_chrom[chrom]):
            if not merged or merged[-1][0] != chrom or start0 > merged[-1][2]:
                merged.append((chrom, start0, end0))
            else:
                prev_chrom, prev_start, prev_end = merged[-1]
                merged[-1] = (prev_chrom, prev_start, max(prev_end, end0))
    return merged


def write_windows(path: Path, windows: list[tuple[str, int, int]]) -> None:
    with path.open("w", encoding="ascii") as handle:
        for chrom, start0, end0 in windows:
            handle.write(f"{chrom}\t{start0}\t{end0}\n")


def artifact(path: Path) -> dict[str, Any]:
    return {"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def find_complete_matrix(repo: Path, sample: str) -> Path:
    matches = sorted((repo / "output" / "canonical" / sample / "paired_full").glob("*_complete_matrix"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected exactly one complete_matrix run for {sample}, found {len(matches)}")
    return matches[0].resolve()


def ensure_index(vcf: Path, command_log: list[list[str]]) -> None:
    csi = Path(str(vcf) + ".csi")
    tbi = Path(str(vcf) + ".tbi")
    if not csi.exists() and not tbi.exists():
        run(["bcftools", "index", "-f", str(vcf)], command_log)


def prepare_sample(
    repo: Path,
    workspace: Path,
    sample_entry: dict[str, Any],
    command_log: list[list[str]],
) -> dict[str, Any]:
    sample = sample_entry["sample"]
    sample_dir = workspace / "truth_split" / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    canonical = find_complete_matrix(repo, sample)
    context_path = canonical / "run_context.json"
    context = json.loads(context_path.read_text(encoding="utf-8"))

    latest_pass = Path(sample_entry["somatic"]["longphase_recalibrated_pass_vcf"]["path"])
    raw_bam = Path(sample_entry["alignment_payload"]["path"])
    raw_bai = Path(sample_entry["alignment_payload"]["index_path"])
    sidecar = Path(sample_entry["read_tags"]["sidecar"]["path"])
    sidecar_index = Path(sample_entry["read_tags"]["index"]["path"])
    truth_vcf = Path(context["truth_vcf"])
    truth_bed = Path(context["truth_bed"]) if context.get("truth_bed") else None
    canonical_fp = canonical / "longphase_s" / "filtered_snv_fp.vcf.gz"
    canonical_tp = canonical / "longphase_s" / "filtered_snv_tp.vcf.gz"
    canonical_matrix_root = canonical / "intersubmod_fp"

    required = [latest_pass, raw_bam, raw_bai, sidecar, sidecar_index, truth_vcf, canonical_fp, canonical_tp]
    if truth_bed is not None:
        required.append(truth_bed)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"{sample} missing required inputs: {missing}")

    latest_snv = sample_dir / f"{sample}.longphase_s.recalibrated.pass.autosomal_biallelic_snv.vcf.gz"
    latest_fp = sample_dir / f"{sample}.latest_lps_pass.truth_fp.vcf.gz"
    latest_tp = sample_dir / f"{sample}.latest_lps_pass.truth_tp.vcf.gz"
    isec_dir = sample_dir / "isec_exact"

    if not latest_snv.exists():
        run(
            [
                "bcftools",
                "view",
                "-f",
                "PASS",
                "-v",
                "snps",
                "-m2",
                "-M2",
                "-r",
                ",".join(AUTOSOMES),
                "-Oz",
                "-o",
                str(latest_snv),
                str(latest_pass),
            ],
            command_log,
        )
    ensure_index(latest_snv, command_log)

    if not latest_fp.exists() or not latest_tp.exists():
        isec_dir.mkdir(parents=True, exist_ok=True)
        isec_outputs = [isec_dir / f"{index:04d}.vcf" for index in range(4)]
        if not all(path.exists() and path.stat().st_size > 0 for path in isec_outputs):
            command = ["bcftools", "isec"]
            if truth_bed is not None:
                command.extend(["-R", str(truth_bed)])
            command.extend(["-p", str(isec_dir), str(latest_snv), str(truth_vcf)])
            run(command, command_log)
        if not latest_fp.exists():
            run(["bcftools", "view", "-Oz", "-o", str(latest_fp), str(isec_dir / "0000.vcf")], command_log)
        if not latest_tp.exists():
            run(["bcftools", "view", "-Oz", "-o", str(latest_tp), str(isec_dir / "0002.vcf")], command_log)
    ensure_index(latest_fp, command_log)
    ensure_index(latest_tp, command_log)

    latest_fp_keys = variant_keys(latest_fp)
    latest_tp_keys = variant_keys(latest_tp)
    canonical_fp_keys = variant_keys(canonical_fp)
    canonical_tp_keys = variant_keys(canonical_tp)
    overlap_fp = latest_fp_keys & canonical_fp_keys
    promoted_fp = latest_fp_keys - canonical_fp_keys
    demoted_canonical_fp = canonical_fp_keys - latest_fp_keys
    fp_windows = merge_windows(latest_fp_keys)
    windows_path = sample_dir / f"{sample}.latest_fp.w5000.merged.bed"
    write_windows(windows_path, fp_windows)

    site_table = sample_dir / f"{sample}.latest_fp_sites.tsv"
    with site_table.open("w", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "chrom", "pos", "ref", "alt", "canonical_fp_matrix_available"])
        for chrom, pos, ref, alt in sorted(
            latest_fp_keys, key=lambda key: (AUTOSOMES.index(key[0]), key[1], key[2], key[3])
        ):
            writer.writerow([sample, chrom, pos, ref, alt, int((chrom, pos, ref, alt) in canonical_fp_keys)])

    return {
        "sample": sample,
        "biological_id": sample_entry["biological_id"],
        "latest_lps_pass_source": artifact(latest_pass),
        "latest_autosomal_biallelic_snv": artifact(latest_snv),
        "latest_truth_fp": artifact(latest_fp),
        "latest_truth_tp": artifact(latest_tp),
        "latest_fp_count": len(latest_fp_keys),
        "latest_tp_count": len(latest_tp_keys),
        "canonical_fp_count": len(canonical_fp_keys),
        "canonical_tp_count": len(canonical_tp_keys),
        "canonical_fp_matrix_count": count_matrices(canonical_matrix_root),
        "latest_fp_overlap_canonical_fp": len(overlap_fp),
        "latest_fp_promoted_or_previously_unmaterialized": len(promoted_fp),
        "canonical_fp_not_in_latest_pass": len(demoted_canonical_fp),
        "latest_fp_overlap_fraction": len(overlap_fp) / len(latest_fp_keys) if latest_fp_keys else None,
        "merged_window_count": len(fp_windows),
        "merged_window_bp": sum(end - start for _chrom, start, end in fp_windows),
        "window_bed": artifact(windows_path),
        "site_table": artifact(site_table),
        "raw_alignment": {"path": str(raw_bam), "index_path": str(raw_bai)},
        "latest_read_tag_sidecar": {"path": str(sidecar), "index_path": str(sidecar_index)},
        "truth_vcf": str(truth_vcf),
        "truth_bed": str(truth_bed) if truth_bed is not None else None,
        "canonical_complete_matrix_run": str(canonical),
        "canonical_fp_vcf": str(canonical_fp),
        "canonical_fp_matrix_root": str(canonical_matrix_root),
        "run_context": str(context_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument(
        "--latest-manifest",
        type=Path,
        default=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
            "20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json"
        ),
    )
    parser.add_argument(
        "--workspace",
        type=Path,
        default=Path(
            "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
            "20260715_single_fp_alt_multicluster_subclone_validation"
        ),
    )
    parser.add_argument("--summary", type=Path, default=None)
    args = parser.parse_args()
    repo = args.repo.resolve()
    workspace = args.workspace.resolve()
    workspace.mkdir(parents=True, exist_ok=True)
    summary_path = args.summary or (
        repo
        / "research"
        / "20260715_single_fp_alt_multicluster_subclone_validation"
        / "results"
        / "latest_input_preflight.json"
    )
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    manifest = json.loads(args.latest_manifest.read_text(encoding="utf-8"))
    by_sample = {entry["sample"]: entry for entry in manifest["samples"]}
    if set(by_sample) != set(DATASETS):
        raise RuntimeError(f"Dataset mismatch: manifest={sorted(by_sample)} expected={sorted(DATASETS)}")

    commands: list[list[str]] = []
    samples = [prepare_sample(repo, workspace, by_sample[sample], commands) for sample in DATASETS]
    summary = {
        "schema_name": "intersubmod.latest_fp_input_preflight",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "task_type": "B_comprehensive_validation",
        "dataset_count": len(samples),
        "biological_sample_count": len({entry["biological_id"] for entry in samples}),
        "latest_manifest": artifact(args.latest_manifest),
        "workspace": str(workspace),
        "scope": "chr1-22 biallelic sSNV; LongPhase-S recalibrated FILTER=PASS; post hoc truth split",
        "truth_isolation": "truth VCF/BED used only after LongPhase-S output; no truth option was supplied to producer",
        "samples": samples,
        "totals": {
            "latest_fp": sum(entry["latest_fp_count"] for entry in samples),
            "latest_tp": sum(entry["latest_tp_count"] for entry in samples),
            "canonical_fp": sum(entry["canonical_fp_count"] for entry in samples),
            "canonical_fp_matrices": sum(entry["canonical_fp_matrix_count"] for entry in samples),
            "latest_fp_overlap_canonical": sum(entry["latest_fp_overlap_canonical_fp"] for entry in samples),
            "latest_fp_promoted_or_previously_unmaterialized": sum(
                entry["latest_fp_promoted_or_previously_unmaterialized"] for entry in samples
            ),
            "canonical_fp_not_in_latest_pass": sum(entry["canonical_fp_not_in_latest_pass"] for entry in samples),
            "merged_window_bp": sum(entry["merged_window_bp"] for entry in samples),
        },
        "commands": commands,
    }
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (workspace / "latest_input_preflight.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"summary": str(summary_path), "totals": summary["totals"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
