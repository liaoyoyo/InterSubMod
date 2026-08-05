#!/usr/bin/env python3
"""Run InterSubMod for every frozen LongPhase-S PASS autosomal biallelic sSNV."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


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


def csv_data_rows(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as handle:
        return max(0, sum(1 for _ in csv.reader(handle)) - 1)


def validate_output(output_dir: Path, expected_sites: int) -> dict[str, Any]:
    summary = output_dir / "significance_summary.csv"
    statistics = output_dir / "significance_statistics.txt"
    status = output_dir / "region_stratification_status.tsv"
    log = output_dir / "inter_sub_mod.log"
    summary_rows = csv_data_rows(summary) if summary.exists() else 0
    reads_files = sum(1 for _ in output_dir.rglob("reads.tsv")) if output_dir.exists() else 0
    matrix_files = (
        sum(1 for path in output_dir.rglob("matrix.csv") if path.parent.name == "BERNOULLI")
        if output_dir.exists()
        else 0
    )
    methylation_files = sum(1 for _ in output_dir.rglob("methylation.csv")) if output_dir.exists() else 0
    status_rows: list[dict[str, str]] = []
    if status.exists():
        with status.open(newline="", encoding="utf-8") as handle:
            status_rows = list(csv.DictReader(handle, delimiter="\t"))
    status_value = status_rows[0].get("status") if len(status_rows) == 1 else None
    accepted_status = status_value in {"SUCCESS", "NOT_APPLICABLE_TUMOR_ONLY"}
    log_text = log.read_text(encoding="utf-8", errors="replace") if log.exists() else ""
    region_failures = sum(
        1 for line in log_text.splitlines() if "[ERROR] Region " in line and " failed:" in line
    )
    passed = (
        summary.exists()
        and statistics.exists()
        and status.exists()
        and accepted_status
        and region_failures == 0
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
        "region_failures_in_log": region_failures,
        "region_stratification_status": status_value,
        "run_level_files": {
            "significance_summary": summary.exists(),
            "significance_statistics": statistics.exists(),
            "region_stratification_status_tsv": status.exists(),
        },
        "pass": passed,
    }


def same_path(recorded: Any, expected: Path) -> bool:
    return bool(
        isinstance(recorded, str)
        and recorded.strip()
        and Path(recorded).resolve() == expected.resolve()
    )


def validate_reused_receipt_provenance(
    receipt_path: Path,
    entry: dict[str, Any],
    binary: Path,
    reference: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Validate the original run receipt before treating an existing output as reusable."""
    if not receipt_path.is_file() or receipt_path.stat().st_size == 0:
        raise RuntimeError(f"Existing output lacks its original nonempty receipt: {receipt_path}")
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RuntimeError(f"Cannot read existing output receipt {receipt_path}: {error}") from error
    if not isinstance(receipt, dict):
        raise RuntimeError(f"Existing output receipt is not a JSON object: {receipt_path}")

    sample = entry["sample"]
    bam = Path(entry["raw_alignment"]["path"])
    bai = Path(entry["raw_alignment_index"]["path"])
    vcf = Path(entry["all_ssnv_vcf"]["path"])
    gates = {
        "schema": receipt.get("schema_name") == "intersubmod.all_ssnv_site_run",
        "sample": receipt.get("sample") == sample,
        "pass": receipt.get("pass") is True,
        "exit_code": receipt.get("exit_code") == 0,
        "original_receipt": receipt.get("reused_existing_validated_output") is False,
        "input_bam": same_path(receipt.get("input_bam"), bam),
        "input_bai": same_path(receipt.get("input_bai"), bai),
        "input_vcf": same_path(receipt.get("input_vcf"), vcf),
        "input_vcf_sha256": receipt.get("input_vcf_sha256") == entry["all_ssnv_vcf"]["sha256"],
        "binary": same_path(receipt.get("binary"), binary),
        "binary_sha256": receipt.get("binary_sha256") == sha256(binary),
        "reference": same_path(receipt.get("reference"), reference),
        "output_dir": same_path(receipt.get("output_dir"), output_dir),
    }
    if "reference_sha256" in receipt:
        gates["reference_sha256"] = receipt.get("reference_sha256") == sha256(reference)
    failed = [name for name, passed in gates.items() if not passed]
    if failed:
        raise RuntimeError(f"{sample} existing receipt provenance failed gates: {failed}")
    return receipt


def run_sample(
    entry: dict[str, Any],
    binary: str,
    reference: str,
    output_root: str,
    threads: int,
) -> dict[str, Any]:
    sample = entry["sample"]
    bam = Path(entry["raw_alignment"]["path"])
    bai = Path(entry["raw_alignment_index"]["path"])
    vcf = Path(entry["all_ssnv_vcf"]["path"])
    expected_sites = int(entry["counts"]["all_ssnv"])
    output_dir = Path(output_root) / sample
    receipt_path = output_dir / "run_receipt.json"
    log_path = output_dir / "inter_sub_mod.log"
    command = [
        binary,
        "-t",
        str(bam),
        "-r",
        reference,
        "-v",
        str(vcf),
        "-o",
        str(output_dir),
        "-w",
        "5000",
        "-j",
        str(threads),
        "--distance-metric",
        "BERNOULLI",
        "--min-common-coverage",
        "3",
        "--nan-distance-strategy",
        "SKIP",
        "--methyl-low",
        "0.2",
        "--methyl-high",
        "0.8",
        "--log-level",
        "info",
    ]
    for path in (bam, bai, vcf, Path(reference), Path(reference + ".fai")):
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    if sha256(vcf) != entry["all_ssnv_vcf"]["sha256"]:
        raise RuntimeError(f"{sample} frozen VCF hash drift")

    if output_dir.exists():
        receipt = validate_reused_receipt_provenance(
            receipt_path,
            entry,
            Path(binary),
            Path(reference),
            output_dir,
        )
        validation = validate_output(output_dir, expected_sites)
        if not validation["pass"]:
            raise RuntimeError(f"{sample} output exists but is incomplete; refusing overwrite: {output_dir}")
        if receipt.get("validation") != validation:
            raise RuntimeError(
                f"{sample} existing output validation no longer matches its original receipt"
            )
        return receipt

    output_dir.mkdir(parents=True, exist_ok=False)
    started_at = now_utc()
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True, check=False)
    exit_code = completed.returncode
    finished_at = now_utc()
    validation = validate_output(output_dir, expected_sites)
    log_tail = log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-40:]

    receipt = {
        "schema_name": "intersubmod.all_ssnv_site_run",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "started_at_utc": started_at,
        "finished_at_utc": finished_at,
        "sample": sample,
        "analysis_role": "all_lps_pass_focal_methyl_primary; latest_HP_PS_posthoc_sidecar",
        "input_bam": str(bam),
        "input_bai": str(bai),
        "input_vcf": str(vcf),
        "input_vcf_sha256": entry["all_ssnv_vcf"]["sha256"],
        "latest_read_tag_sidecar": entry["latest_read_tag_sidecar"],
        "binary": binary,
        "binary_sha256": sha256(Path(binary)),
        "reference": reference,
        "command": command,
        "output_dir": str(output_dir),
        "log_path": str(log_path) if log_path.exists() else None,
        "log_tail": log_tail,
        "exit_code": exit_code,
        "reused_existing_validated_output": False,
        "validation": validation,
        "pass": exit_code == 0 and validation["pass"],
        "interpretation": (
            "Raw BAM HP is not consumed as latest tag evidence. Allele/MM/ML/distance are tag-independent; "
            "frozen LongPhase-S HP/PS is exact-joined during downstream focal analysis."
        ),
    }
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if not receipt["pass"]:
        raise RuntimeError(f"{sample} InterSubMod failed: exit={exit_code} validation={validation}")
    return receipt


def main() -> None:
    topic_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=topic_root / "results" / "all_ssnv_input_manifest.json")
    parser.add_argument("--binary", type=Path, default=Path("build/bin/inter_sub_mod"))
    parser.add_argument("--reference", type=Path, default=Path("/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta"))
    parser.add_argument("--output-root", type=Path, default=None)
    parser.add_argument("--sample", action="append", choices=DATASETS)
    parser.add_argument("--workers", type=int, default=7)
    parser.add_argument("--threads-per-sample", type=int, default=6)
    parser.add_argument("--minimum-free-gib", type=int, default=300)
    parser.add_argument("--summary", type=Path, default=topic_root / "results" / "all_ssnv_intersubmod_batch.json")
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    if manifest.get("pass") is not True or manifest["totals"]["all_ssnv"] != 469_849:
        raise RuntimeError("Input manifest is not a passing full-scope manifest")
    selected = set(args.sample or DATASETS)
    entries = [entry for entry in manifest["samples"] if entry["sample"] in selected]
    output_root = args.output_root or (Path(manifest["workspace"]) / "intersubmod_all_ssnv_v1")
    output_root.mkdir(parents=True, exist_ok=True)
    free_gib = shutil.disk_usage(output_root).free / (1024**3)
    if free_gib < args.minimum_free_gib:
        raise RuntimeError(f"Insufficient free disk: {free_gib:.1f} GiB < {args.minimum_free_gib} GiB")
    requested_threads = args.workers * args.threads_per_sample
    if requested_threads > (os.cpu_count() or requested_threads) + 2:
        raise RuntimeError(f"Requested {requested_threads} threads exceeds CPU guardrail")

    receipts: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = {
            executor.submit(
                run_sample,
                entry,
                str(args.binary.resolve()),
                str(args.reference.resolve()),
                str(output_root),
                args.threads_per_sample,
            ): entry["sample"]
            for entry in entries
        }
        for future in as_completed(futures):
            sample = futures[future]
            try:
                receipt = future.result()
                receipts.append(receipt)
                print(
                    f"[{sample}] PASS sites={receipt['validation']['expected_vcf_sites']} "
                    f"matrices={receipt['validation']['bernoulli_matrix_files']} ",
                    flush=True,
                )
            except Exception as error:
                failures.append({"sample": sample, "error": repr(error)})
                print(f"[{sample}] FAIL {error!r}", flush=True)

    receipts.sort(key=lambda value: DATASETS.index(value["sample"]))
    totals = {
        key: sum(receipt["validation"][key] for receipt in receipts)
        for key in ("expected_vcf_sites", "summary_rows", "reads_files", "bernoulli_matrix_files", "methylation_files")
    }
    batch = {
        "schema_name": "intersubmod.all_ssnv_site_batch",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "manifest": str(args.manifest.resolve()),
        "output_root": str(output_root.resolve()),
        "requested_samples": sorted(selected, key=DATASETS.index),
        "workers": args.workers,
        "threads_per_sample": args.threads_per_sample,
        "requested_threads": requested_threads,
        "free_gib_at_launch": free_gib,
        "receipts": receipts,
        "failures": failures,
        "totals": totals,
        "pass": (
            len(receipts) == len(entries)
            and not failures
            and all(receipt["pass"] for receipt in receipts)
            and totals["expected_vcf_sites"] == sum(entry["counts"]["all_ssnv"] for entry in entries)
        ),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(batch, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(args.summary), "totals": totals, "pass": batch["pass"]}, indent=2))
    raise SystemExit(0 if batch["pass"] else 1)


if __name__ == "__main__":
    main()
