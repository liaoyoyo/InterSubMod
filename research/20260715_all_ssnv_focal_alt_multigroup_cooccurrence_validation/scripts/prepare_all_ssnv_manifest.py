#!/usr/bin/env python3
"""Freeze and verify the all-sSNV inputs for focal-ALT methylation analysis."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pysam


DATASETS = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]
LAYERED_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260713_layered_reconstruction_v3_raw_all_lps_pass_v5"
)
SOURCE_PREFLIGHT = Path(
    "/big7_disk/liaoyoyo2001/InterSubMod/research/"
    "20260715_single_fp_alt_multicluster_subclone_validation/results/latest_input_preflight.json"
)
WORKSPACE = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path, hash_content: bool = True) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "mtime_ns": path.stat().st_mtime_ns,
    }
    if hash_content:
        result["sha256"] = sha256(path)
    return result


def variant_keys(path: Path) -> set[tuple[str, int, str, str]]:
    keys: set[tuple[str, int, str, str]] = set()
    with pysam.VariantFile(str(path)) as source:
        for record in source:
            if len(record.ref) != 1 or not record.alts or len(record.alts) != 1 or len(record.alts[0]) != 1:
                raise RuntimeError(f"Non-biallelic SNV in frozen all-sSNV input: {path} {record}")
            if record.contig not in {f"chr{index}" for index in range(1, 23)}:
                raise RuntimeError(f"Non-autosomal record in frozen all-sSNV input: {path} {record}")
            key = (record.contig, int(record.pos), record.ref.upper(), record.alts[0].upper())
            if key in keys:
                raise RuntimeError(f"Duplicate VCF key: {path} {key}")
            keys.add(key)
    return keys


def ledger_pass_keys(path: Path) -> tuple[set[tuple[str, int, str, str]], Counter[str]]:
    keys: set[tuple[str, int, str, str]] = set()
    branches: Counter[str] = Counter()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["longphase_recalibrated_filter"] != "PASS":
                continue
            chrom = row["chrom"]
            if chrom not in {f"chr{index}" for index in range(1, 23)}:
                continue
            ref, alt = row["ref"].upper(), row["alt"].upper()
            if len(ref) != 1 or len(alt) != 1:
                continue
            key = (chrom, int(row["pos"]), ref, alt)
            if key in keys:
                raise RuntimeError(f"Duplicate ledger key: {path} {key}")
            keys.add(key)
            branches[row["ssnv_branch"]] += 1
    return keys, branches


def validate_index(path: Path, explicit_index: Path | None = None) -> Path:
    candidates = [explicit_index] if explicit_index else []
    candidates.extend([Path(str(path) + ".bai"), Path(str(path) + ".csi"), Path(str(path) + ".tbi")])
    for candidate in candidates:
        if candidate is not None and candidate.exists() and candidate.stat().st_size > 0:
            return candidate
    raise FileNotFoundError(f"No non-empty index found for {path}")


def prepare_sample(entry: dict[str, Any], workspace: Path) -> dict[str, Any]:
    sample = entry["sample"]
    all_vcf = Path(entry["latest_autosomal_biallelic_snv"]["path"])
    fp_vcf = Path(entry["latest_truth_fp"]["path"])
    tp_vcf = Path(entry["latest_truth_tp"]["path"])
    raw_bam = Path(entry["raw_alignment"]["path"])
    raw_bai = Path(entry["raw_alignment"]["index_path"])
    sidecar = Path(entry["latest_read_tag_sidecar"]["path"])
    sidecar_index = Path(entry["latest_read_tag_sidecar"]["index_path"])
    ledger = LAYERED_ROOT / "samples" / sample / f"ssnv_site_ledger_{sample}.tsv.gz"
    ledger_index = Path(str(ledger) + ".tbi")
    topology = LAYERED_ROOT / "samples" / sample / f"layered_reconstruction_{sample}.json"
    region_view = LAYERED_ROOT / "samples" / sample / f"layered_region_view_{sample}.json"
    required = [
        all_vcf,
        fp_vcf,
        tp_vcf,
        raw_bam,
        raw_bai,
        sidecar,
        sidecar_index,
        ledger,
        ledger_index,
        topology,
        region_view,
    ]
    missing = [str(path) for path in required if not path.exists() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError(f"{sample} missing/empty inputs: {missing}")
    all_index = validate_index(all_vcf)
    fp_index = validate_index(fp_vcf)
    tp_index = validate_index(tp_vcf)

    all_keys = variant_keys(all_vcf)
    fp_keys = variant_keys(fp_vcf)
    tp_keys = variant_keys(tp_vcf)
    ledger_keys, branch_counts = ledger_pass_keys(ledger)
    if fp_keys & tp_keys:
        raise RuntimeError(f"{sample} truth FP/TP sets overlap")
    if not fp_keys <= all_keys or not tp_keys <= all_keys:
        raise RuntimeError(f"{sample} truth labels are not subsets of all-sSNV input")
    if ledger_keys != all_keys:
        raise RuntimeError(
            f"{sample} VCF/ledger mismatch all={len(all_keys)} ledger={len(ledger_keys)} "
            f"vcf_only={len(all_keys-ledger_keys)} ledger_only={len(ledger_keys-all_keys)}"
        )
    expected_branches = {"retained", "max_snv_excluded", "positional_singleton"}
    if set(branch_counts) - expected_branches:
        raise RuntimeError(f"{sample} unexpected PASS ledger branches: {dict(branch_counts)}")

    sample_dir = workspace / "manifest" / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    site_table = sample_dir / f"{sample}.all_ssnv_site_manifest.tsv.gz"
    with gzip.open(site_table, "wt", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "chrom", "pos", "ref", "alt", "truth_label"])
        rank = {f"chr{index}": index for index in range(1, 23)}
        for key in sorted(all_keys, key=lambda value: (rank[value[0]], value[1], value[2], value[3])):
            truth = "TP" if key in tp_keys else "FP" if key in fp_keys else "UNASSESSED"
            writer.writerow([sample, *key, truth])

    return {
        "sample": sample,
        "biological_id": entry["biological_id"],
        "all_ssnv_vcf": artifact(all_vcf),
        "all_ssnv_vcf_index": artifact(all_index),
        "truth_fp_vcf": artifact(fp_vcf),
        "truth_fp_vcf_index": artifact(fp_index),
        "truth_tp_vcf": artifact(tp_vcf),
        "truth_tp_vcf_index": artifact(tp_index),
        "raw_alignment": artifact(raw_bam, hash_content=False),
        "raw_alignment_index": artifact(raw_bai, hash_content=False),
        "latest_read_tag_sidecar": artifact(sidecar, hash_content=False),
        "latest_read_tag_sidecar_index": artifact(sidecar_index),
        "site_ledger": artifact(ledger),
        "site_ledger_index": artifact(ledger_index),
        "layered_reconstruction": artifact(topology),
        "layered_region_view": artifact(region_view),
        "site_manifest": artifact(site_table),
        "counts": {
            "all_ssnv": len(all_keys),
            "truth_tp": len(tp_keys),
            "truth_fp": len(fp_keys),
            "truth_unassessed": len(all_keys - tp_keys - fp_keys),
            "ledger_branches": dict(sorted(branch_counts.items())),
        },
        "checks": {
            "all_vcf_unique_autosomal_biallelic_snv": True,
            "truth_fp_tp_disjoint": True,
            "truth_labels_subset_of_all": True,
            "all_vcf_equals_latest_pass_ledger": True,
            "raw_bam_and_index_nonempty": True,
            "sidecar_and_index_nonempty": True,
        },
    }


def main() -> None:
    topic_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-preflight", type=Path, default=SOURCE_PREFLIGHT)
    parser.add_argument("--workspace", type=Path, default=WORKSPACE)
    parser.add_argument("--summary", type=Path, default=topic_root / "results" / "all_ssnv_input_manifest.json")
    args = parser.parse_args()

    source = json.loads(args.source_preflight.read_text(encoding="utf-8"))
    if [entry["sample"] for entry in source["samples"]] != DATASETS:
        raise RuntimeError("Source preflight dataset order/content mismatch")
    args.workspace.mkdir(parents=True, exist_ok=True)
    samples = [prepare_sample(entry, args.workspace) for entry in source["samples"]]
    totals = {
        key: sum(sample["counts"][key] for sample in samples)
        for key in ("all_ssnv", "truth_tp", "truth_fp", "truth_unassessed")
    }
    aggregate_branches: Counter[str] = Counter()
    for sample in samples:
        aggregate_branches.update(sample["counts"]["ledger_branches"])
    totals["ledger_branches"] = dict(sorted(aggregate_branches.items()))
    if totals["all_ssnv"] != 469_849:
        raise RuntimeError(f"Frozen all-sSNV total drifted: {totals['all_ssnv']} != 469849")
    if totals["truth_tp"] != 335_296 or totals["truth_fp"] != 7_745:
        raise RuntimeError(f"Frozen truth totals drifted: {totals}")
    if totals["truth_unassessed"] != 126_808:
        raise RuntimeError(f"Frozen unassessed total drifted: {totals}")

    disk = shutil.disk_usage(args.workspace)
    summary = {
        "schema_name": "intersubmod.all_ssnv_focal_alt_input_manifest",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "task_type": "B_comprehensive_validation",
        "scope": "7 datasets; chr1-22; LongPhase-S recalibrated FILTER=PASS; biallelic sSNV",
        "grouping_contract": "truth-blind focal-ALT methyl grouping; truth and co-occurrence labels post hoc",
        "tag_contract": "raw BAM supplies allele+MM/ML; frozen LongPhase-S HP/PS sidecar is exact-joined post hoc",
        "source_preflight": artifact(args.source_preflight),
        "layered_root": artifact(LAYERED_ROOT / "input_manifest.snapshot.json"),
        "workspace": str(args.workspace.resolve()),
        "resources": {
            "logical_cpus": os.cpu_count(),
            "disk_total_bytes": disk.total,
            "disk_used_bytes": disk.used,
            "disk_free_bytes": disk.free,
        },
        "samples": samples,
        "totals": totals,
        "pass": True,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    workspace_copy = args.workspace / "all_ssnv_input_manifest.json"
    workspace_copy.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(args.summary), "totals": totals, "pass": True}, indent=2))


if __name__ == "__main__":
    main()

