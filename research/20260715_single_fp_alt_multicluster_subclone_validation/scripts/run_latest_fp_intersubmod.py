#!/usr/bin/env python3
"""Run InterSubMod on every latest LongPhase-S PASS truth-FP site."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
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


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def vcf_count(path: Path) -> int:
    with pysam.VariantFile(str(path)) as variants:
        return sum(1 for _ in variants)


def csv_data_rows(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as handle:
        return max(0, sum(1 for _ in csv.reader(handle)) - 1)


def validate_output(output_dir: Path, expected_sites: int) -> dict[str, Any]:
    summary = output_dir / "significance_summary.csv"
    statistics = output_dir / "significance_statistics.txt"
    summary_rows = csv_data_rows(summary) if summary.exists() else 0
    reads_files = sum(1 for _ in output_dir.rglob("reads.tsv")) if output_dir.exists() else 0
    matrix_files = sum(
        1 for path in output_dir.rglob("matrix.csv") if path.parent.name == "BERNOULLI"
    ) if output_dir.exists() else 0
    methylation_files = sum(1 for _ in output_dir.rglob("methylation.csv")) if output_dir.exists() else 0
    phylo_summaries = sum(1 for _ in output_dir.rglob("phylo_groups_summary.json")) if output_dir.exists() else 0
    passed = (
        summary.exists()
        and statistics.exists()
        and 0 < summary_rows <= expected_sites
        and reads_files == expected_sites
        and matrix_files == expected_sites
        and methylation_files == expected_sites
    )
    return {
        "expected_vcf_sites": expected_sites,
        "summary_rows": summary_rows,
        "summary_rows_not_emitted": expected_sites - summary_rows,
        "reads_files": reads_files,
        "bernoulli_matrix_files": matrix_files,
        "methylation_files": methylation_files,
        "phylo_groups_summary_files": phylo_summaries,
        "pass": passed,
    }


def load_materialization_receipts(
    results_dir: Path, materialization_label: str
) -> dict[str, dict[str, Any]]:
    candidates: dict[str, list[dict[str, Any]]] = {sample: [] for sample in DATASETS}
    for path in sorted(results_dir.glob("*materialization*.json")):
        batch = json.loads(path.read_text(encoding="utf-8"))
        for receipt in batch.get("receipts", []):
            receipt_label = receipt.get("materialization_label", "latest_lps_fp_w5000")
            if (
                receipt.get("pass") is True
                and receipt.get("sample") in candidates
                and receipt_label == materialization_label
            ):
                receipt = dict(receipt)
                receipt["batch_summary"] = str(path)
                candidates[receipt["sample"]].append(receipt)
    selected: dict[str, dict[str, Any]] = {}
    for sample, values in candidates.items():
        if values:
            selected[sample] = max(values, key=lambda value: value.get("created_at_utc", ""))
    return selected


def run_sample(
    sample_entry: dict[str, Any],
    materialization: dict[str, Any],
    binary: str,
    reference: str,
    output_root: str,
    threads: int,
    vcf_field: str,
    analysis_role: str,
) -> dict[str, Any]:
    sample = sample_entry["sample"]
    bam = Path(materialization["output_bam"])
    bai = Path(materialization["output_bai"])
    vcf = Path(sample_entry[vcf_field]["path"])
    output_dir = Path(output_root) / sample
    receipt_path = output_dir / "run_receipt.json"
    log_path = output_dir / "inter_sub_mod.log"
    expected_sites = vcf_count(vcf)
    command = [
        binary,
        "-t", str(bam),
        "-r", reference,
        "-v", str(vcf),
        "-o", str(output_dir),
        "-w", "5000",
        "-j", str(threads),
        "--distance-metric", "BERNOULLI",
        "--min-common-coverage", "3",
        "--nan-distance-strategy", "SKIP",
        "--methyl-low", "0.2",
        "--methyl-high", "0.8",
        "--log-level", "info",
    ]
    for path in (bam, bai, vcf, Path(reference), Path(reference + ".fai")):
        if not path.exists():
            raise FileNotFoundError(path)
    if not any(Path(str(vcf) + suffix).exists() for suffix in (".csi", ".tbi")):
        raise FileNotFoundError(f"No .csi or .tbi index for {vcf}")

    reused = False
    if output_dir.exists():
        validation = validate_output(output_dir, expected_sites)
        if validation["pass"]:
            reused = True
            exit_code = 0
            started_at = None
            finished_at = now_utc()
            log_tail: list[str] = []
        else:
            raise RuntimeError(f"{sample} output directory exists but is incomplete; refusing overwrite: {output_dir}")
    else:
        output_dir.mkdir(parents=True, exist_ok=False)
        started_at = now_utc()
        with log_path.open("w", encoding="utf-8") as log:
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True, check=False)
        exit_code = completed.returncode
        finished_at = now_utc()
        validation = validate_output(output_dir, expected_sites)
        log_tail = log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-30:]

    passed = exit_code == 0 and validation["pass"]
    receipt = {
        "schema_name": "intersubmod.latest_lps_pass_site_run",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "started_at_utc": started_at,
        "finished_at_utc": finished_at,
        "sample": sample,
        "analysis_role": analysis_role,
        "input_bam": str(bam),
        "input_bai": str(bai),
        "input_vcf": str(vcf),
        "input_vcf_sha256": sha256(vcf),
        "materialization_receipt_batch": materialization["batch_summary"],
        "binary": binary,
        "binary_sha256": sha256(Path(binary)),
        "reference": reference,
        "command": command,
        "output_dir": str(output_dir),
        "log_path": str(log_path) if log_path.exists() else None,
        "log_tail": log_tail,
        "exit_code": exit_code,
        "reused_existing_validated_output": reused,
        "validation": validation,
        "pass": passed,
    }
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if not passed:
        raise RuntimeError(f"{sample} InterSubMod validation failed: {validation}; exit={exit_code}")
    return receipt


def main() -> None:
    script_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", type=Path, default=script_root / "results" / "latest_input_preflight.json")
    parser.add_argument("--results-dir", type=Path, default=script_root / "results")
    parser.add_argument("--binary", type=Path, default=Path("build/bin/inter_sub_mod"))
    parser.add_argument("--reference", type=Path, default=Path("/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"))
    parser.add_argument("--sample", action="append", choices=DATASETS)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--threads-per-sample", type=int, default=7)
    parser.add_argument("--materialization-label", default="latest_lps_fp_w5000")
    parser.add_argument("--vcf-field", default="latest_truth_fp")
    parser.add_argument("--analysis-role", default="truth_fp_full_scope")
    parser.add_argument("--output-root", type=Path, default=None)
    parser.add_argument("--summary", type=Path, default=None)
    args = parser.parse_args()

    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    selected = set(args.sample or DATASETS)
    entries = [entry for entry in preflight["samples"] if entry["sample"] in selected]
    materializations = load_materialization_receipts(args.results_dir, args.materialization_label)
    missing = sorted(selected - set(materializations), key=DATASETS.index)
    if missing:
        raise RuntimeError(f"Missing passing latest-tag materialization receipts: {missing}")
    output_root = args.output_root or (Path(preflight["workspace"]) / "intersubmod_latest_fp_v1")
    output_root.mkdir(parents=True, exist_ok=True)

    receipts: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = {
            executor.submit(
                run_sample,
                entry,
                materializations[entry["sample"]],
                str(args.binary.resolve()),
                str(args.reference),
                str(output_root),
                args.threads_per_sample,
                args.vcf_field,
                args.analysis_role,
            ): entry["sample"]
            for entry in entries
        }
        for future in as_completed(futures):
            sample = futures[future]
            try:
                receipt = future.result()
                receipts.append(receipt)
                validation = receipt["validation"]
                print(
                    f"[{sample}] PASS sites={validation['summary_rows']} matrices={validation['bernoulli_matrix_files']} "
                    f"reused={receipt['reused_existing_validated_output']}",
                    flush=True,
                )
            except Exception as error:
                failures.append({"sample": sample, "error": repr(error)})
                print(f"[{sample}] FAIL {error!r}", flush=True)

    receipts.sort(key=lambda value: DATASETS.index(value["sample"]))
    summary = {
        "schema_name": "intersubmod.latest_lps_pass_site_batch",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "preflight": str(args.preflight),
        "output_root": str(output_root),
        "requested_samples": sorted(selected, key=DATASETS.index),
        "workers": args.workers,
        "threads_per_sample": args.threads_per_sample,
        "materialization_label": args.materialization_label,
        "vcf_field": args.vcf_field,
        "analysis_role": args.analysis_role,
        "receipts": receipts,
        "failures": failures,
        "totals": {
            "expected_vcf_sites": sum(r["validation"]["expected_vcf_sites"] for r in receipts),
            "summary_rows": sum(r["validation"]["summary_rows"] for r in receipts),
            "reads_files": sum(r["validation"]["reads_files"] for r in receipts),
            "bernoulli_matrix_files": sum(r["validation"]["bernoulli_matrix_files"] for r in receipts),
            "methylation_files": sum(r["validation"]["methylation_files"] for r in receipts),
        },
        "pass": len(receipts) == len(entries) and not failures and all(r["pass"] for r in receipts),
    }
    summary_path = args.summary or (args.results_dir / "latest_fp_intersubmod_batch.json")
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(summary_path), "totals": summary["totals"], "pass": summary["pass"]}, indent=2))
    raise SystemExit(0 if summary["pass"] else 1)


if __name__ == "__main__":
    main()
