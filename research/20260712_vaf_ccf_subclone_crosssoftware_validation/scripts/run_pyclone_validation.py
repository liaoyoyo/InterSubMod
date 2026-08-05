#!/usr/bin/env python3
"""Run reproducible PyClone-VI fits and preserve command/log/QA receipts."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time
from typing import Dict, List, Mapping, Sequence


SEED = 20260712
NUM_CLUSTERS = 40
DEFAULT_NUM_RESTARTS = 20
DENSITY = "beta-binomial"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def read_input_shape(path: Path) -> Mapping[str, object]:
    rows = 0
    mutations = set()
    samples = set()
    conservation_fail = 0
    pairs = set()
    duplicate_pairs = 0
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "mutation_id", "sample_id", "ref_counts", "alt_counts", "normal_cn",
            "major_cn", "minor_cn", "tumour_content",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} misses PyClone columns {sorted(missing)}")
        for row in reader:
            rows += 1
            mutation = row["mutation_id"]
            sample = row["sample_id"]
            pair = (mutation, sample)
            if pair in pairs:
                duplicate_pairs += 1
            pairs.add(pair)
            mutations.add(mutation)
            samples.add(sample)
            if int(row["ref_counts"]) + int(row["alt_counts"]) <= 0:
                conservation_fail += 1
    expected = len(mutations) * len(samples)
    return {
        "rows": rows,
        "mutations": len(mutations),
        "samples": sorted(samples),
        "expected_complete_rows": expected,
        "complete_matrix": rows == expected,
        "duplicate_mutation_sample_pairs": duplicate_pairs,
        "nonpositive_depth_rows": conservation_fail,
    }


def read_result_shape(path: Path) -> Mapping[str, object]:
    rows = 0
    mutations = set()
    samples = set()
    clusters = set()
    probability_out_of_range = 0
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "mutation_id", "sample_id", "cluster_id", "cellular_prevalence",
            "cellular_prevalence_std", "cluster_assignment_prob",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} misses result columns {sorted(missing)}")
        for row in reader:
            rows += 1
            mutations.add(row["mutation_id"])
            samples.add(row["sample_id"])
            clusters.add(int(row["cluster_id"]))
            probability = float(row["cluster_assignment_prob"])
            if not 0 <= probability <= 1:
                probability_out_of_range += 1
    return {
        "rows": rows,
        "mutations": len(mutations),
        "samples": sorted(samples),
        "clusters": len(clusters),
        "cluster_ids": sorted(clusters),
        "probability_out_of_range": probability_out_of_range,
    }


def parse_args() -> argparse.Namespace:
    topic = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=topic / "data" / "pyclone_inputs")
    parser.add_argument("--output-dir", type=Path, default=topic / "runs" / "pyclone_vi")
    parser.add_argument(
        "--pyclone",
        type=Path,
        default=topic / "external_tools" / "pyclone_vi_env" / "bin" / "pyclone-vi",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--num-restarts", type=int, default=DEFAULT_NUM_RESTARTS)
    parser.add_argument(
        "--batch-id",
        default="default",
        help="Receipt filename suffix so concurrent independent batches do not overwrite each other",
    )
    parser.add_argument("--include", action="append", default=[], help="fnmatch-style bundle glob; repeatable")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def select_inputs(input_dir: Path, patterns: Sequence[str]) -> List[Path]:
    candidates = sorted(input_dir.glob("*.pyclone_input.tsv.gz"))
    if not patterns:
        return candidates
    import fnmatch
    return [path for path in candidates if any(fnmatch.fnmatch(path.name, pattern) for pattern in patterns)]


def run_logged(command: Sequence[str], stdout_path: Path, stderr_path: Path) -> int:
    with stdout_path.open("w") as stdout_handle, stderr_path.open("w") as stderr_handle:
        process = subprocess.run(command, stdout=stdout_handle, stderr=stderr_handle, check=False)
    return process.returncode


def main() -> int:
    args = parse_args()
    if not args.pyclone.is_file():
        raise FileNotFoundError(f"PyClone-VI executable is unavailable: {args.pyclone}")
    if args.threads < 1:
        raise ValueError("--threads must be >=1")
    if args.num_restarts < 1:
        raise ValueError("--num-restarts must be >=1")
    if not args.batch_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("--batch-id may contain only letters, digits, underscore, and hyphen")
    inputs = select_inputs(args.input_dir, args.include)
    if not inputs:
        raise ValueError("No PyClone input bundles matched")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    version = subprocess.run([str(args.pyclone), "--version"], capture_output=True, text=True, check=True).stdout.strip()
    overall_start = time.monotonic()
    jobs: List[Mapping[str, object]] = []
    any_failure = False

    for input_path in inputs:
        bundle_id = input_path.name.removesuffix(".pyclone_input.tsv.gz")
        qa_path = args.input_dir / f"{bundle_id}.qa.json"
        input_qa = json.loads(qa_path.read_text())
        if input_qa["status"] != "PASS":
            raise ValueError(f"Refusing non-PASS input {bundle_id}: {input_qa['status']}")
        job_dir = args.output_dir / bundle_id
        job_dir.mkdir(parents=True, exist_ok=True)
        model_path = job_dir / "model.h5"
        results_path = job_dir / "results.tsv.gz"
        status_path = job_dir / "status.json"
        fit_stdout = job_dir / "fit.stdout.log"
        fit_stderr = job_dir / "fit.stderr.log"
        write_stdout = job_dir / "write_results.stdout.log"
        write_stderr = job_dir / "write_results.stderr.log"
        input_shape = read_input_shape(input_path)
        if not input_shape["complete_matrix"] or input_shape["duplicate_mutation_sample_pairs"]:
            raise ValueError(f"Input matrix QA failed for {bundle_id}: {input_shape}")
        if status_path.exists() and not args.force:
            old_status = json.loads(status_path.read_text())
            if old_status.get("status") == "PASS" and results_path.is_file():
                jobs.append(old_status)
                print(f"SKIP complete {bundle_id}")
                continue
        fit_command = [
            str(args.pyclone), "fit", "-i", str(input_path.resolve()), "-o", str(model_path.resolve()),
            "-c", str(NUM_CLUSTERS), "-d", DENSITY, "-r", str(args.num_restarts),
            "-t", str(args.threads), "--seed", str(SEED),
        ]
        write_command = [
            str(args.pyclone), "write-results-file", "-i", str(model_path.resolve()),
            "-o", str(results_path.resolve()), "-c",
        ]
        receipt: Dict[str, object] = {
            "schema_name": "intersubmod.pyclone_vi_fit_receipt",
            "schema_version": "1.0.0",
            "bundle_id": bundle_id,
            "status": "RUNNING",
            "started_at_utc": utc_now(),
            "runner": str(Path(__file__).resolve()),
            "runner_sha256_at_job_start": sha256_file(Path(__file__).resolve()),
            "pyclone_version": version,
            "parameters": {
                "density": DENSITY,
                "num_clusters": NUM_CLUSTERS,
                "num_restarts": args.num_restarts,
                "num_threads": args.threads,
                "seed": SEED,
            },
            "fit_command_argv": fit_command,
            "fit_command_shell_escaped": shlex.join(fit_command),
            "write_results_command_argv": write_command,
            "write_results_command_shell_escaped": shlex.join(write_command),
            "input": {
                "path": str(input_path.resolve()),
                "sha256": sha256_file(input_path),
                "shape": input_shape,
                "qa_path": str(qa_path.resolve()),
                "qa_sha256": sha256_file(qa_path),
            },
            "outputs": {
                "model": str(model_path.resolve()),
                "results": str(results_path.resolve()),
                "fit_stdout": str(fit_stdout.resolve()),
                "fit_stderr": str(fit_stderr.resolve()),
                "write_results_stdout": str(write_stdout.resolve()),
                "write_results_stderr": str(write_stderr.resolve()),
            },
        }
        atomic_json(status_path, receipt)
        print(f"RUN {bundle_id}: mutations={input_shape['mutations']} samples={','.join(input_shape['samples'])}", flush=True)
        start = time.monotonic()
        fit_return = run_logged(fit_command, fit_stdout, fit_stderr)
        receipt["fit_exit_code"] = fit_return
        if fit_return == 0:
            write_return = run_logged(write_command, write_stdout, write_stderr)
        else:
            write_return = None
        receipt["write_results_exit_code"] = write_return
        receipt["elapsed_seconds"] = round(time.monotonic() - start, 3)
        receipt["finished_at_utc"] = utc_now()
        if fit_return == 0 and write_return == 0 and model_path.is_file() and results_path.is_file():
            result_shape = read_result_shape(results_path)
            expected_rows = input_shape["mutations"] * len(input_shape["samples"])
            result_ok = (
                result_shape["rows"] == expected_rows
                and result_shape["mutations"] == input_shape["mutations"]
                and result_shape["samples"] == input_shape["samples"]
                and result_shape["probability_out_of_range"] == 0
            )
            receipt["result_shape"] = result_shape
            receipt["outputs"]["model_size_bytes"] = model_path.stat().st_size
            receipt["outputs"]["model_sha256"] = sha256_file(model_path)
            receipt["outputs"]["results_size_bytes"] = results_path.stat().st_size
            receipt["outputs"]["results_sha256"] = sha256_file(results_path)
            receipt["status"] = "PASS" if result_ok else "FAIL"
            if not result_ok:
                receipt["failure"] = f"result matrix mismatch; expected_rows={expected_rows} result={result_shape}"
        else:
            receipt["status"] = "FAIL"
            receipt["failure"] = "fit or write-results failed; inspect preserved stderr logs"
        atomic_json(status_path, receipt)
        jobs.append(receipt)
        if receipt["status"] != "PASS":
            any_failure = True
            print(f"FAIL {bundle_id}; see {fit_stderr}", file=sys.stderr, flush=True)
        else:
            print(
                f"PASS {bundle_id}: clusters={receipt['result_shape']['clusters']} elapsed={receipt['elapsed_seconds']}s",
                flush=True,
            )

    manifest = {
        "schema_name": "intersubmod.pyclone_vi_run_manifest",
        "schema_version": "1.0.0",
        "created_at_utc": utc_now(),
        "runner": str(Path(__file__).resolve()),
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "pyclone_executable": str(args.pyclone.resolve()),
        "pyclone_version": version,
        "parameters": {
            "density": DENSITY,
            "num_clusters": NUM_CLUSTERS,
            "num_restarts": args.num_restarts,
            "num_threads": args.threads,
            "seed": SEED,
        },
        "total_elapsed_seconds": round(time.monotonic() - overall_start, 3),
        "jobs": jobs,
    }
    manifest_path = args.output_dir / f"command_manifest.{args.batch_id}.json"
    atomic_json(manifest_path, manifest)
    summary_path = args.output_dir / f"run_summary.{args.batch_id}.tsv"
    with summary_path.open("w") as handle:
        handle.write("bundle_id\tstatus\tmutations\tsamples\tclusters\telapsed_seconds\tresults_path\n")
        for job in jobs:
            result_shape = job.get("result_shape", {})
            handle.write(
                f"{job['bundle_id']}\t{job['status']}\t{job['input']['shape']['mutations']}\t"
                f"{','.join(job['input']['shape']['samples'])}\t{result_shape.get('clusters', '')}\t"
                f"{job.get('elapsed_seconds', '')}\t{job['outputs']['results']}\n"
            )
    print(f"Command manifest: {manifest_path}")
    print(f"Run summary: {summary_path}")
    return 2 if any_failure else 0


if __name__ == "__main__":
    raise SystemExit(main())
