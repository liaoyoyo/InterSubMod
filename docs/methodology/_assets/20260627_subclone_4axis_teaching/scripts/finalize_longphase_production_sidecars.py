#!/usr/bin/env python3
"""Fail-closed closeout for the active LongPhase-S production sidecar run."""

import argparse
import csv
import hashlib
import json
import os
import shutil
import stat
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pysam


EXPECTED_HP = {".", "1", "2", "3", "4", "1-1", "2-1", "1-2", "2-2"}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def read_hash_manifest(path):
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, raw_path = line.split(None, 1)
        artifact = Path(raw_path.lstrip("* "))
        if artifact in values:
            raise RuntimeError(f"duplicate hash entry: {artifact}")
        if not artifact.is_file():
            raise RuntimeError(f"hashed artifact missing: {artifact}")
        observed = sha256(artifact)
        if observed != digest:
            raise RuntimeError(f"hash mismatch: {artifact}")
        values[artifact.resolve()] = digest
    if not values:
        raise RuntimeError(f"empty hash manifest: {path}")
    return values


def merge_hashes(target, additions):
    for artifact, digest in additions.items():
        if artifact in target and target[artifact] != digest:
            raise RuntimeError(f"conflicting recorded hashes: {artifact}")
        target[artifact] = digest


def variant_keys(path):
    keys = Counter()
    pass_keys = Counter()
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            key = (record.contig, int(record.pos), record.ref, tuple(record.alts or ()))
            keys[key] += 1
            if list(record.filter.keys()) == ["PASS"]:
                pass_keys[key] += 1
    return keys, pass_keys


def command_version(command):
    proc = subprocess.run(command, text=True, capture_output=True, check=False)
    text = (proc.stdout or proc.stderr).strip().splitlines()
    return {"command": command, "exit_code": proc.returncode, "first_line": text[0] if text else ""}


def stat_mtime(path):
    proc = subprocess.run(["stat", "-c", "%y", str(path)], text=True, capture_output=True, check=True)
    return proc.stdout.strip()


def validate_inventory(path, sample):
    expected = {
        "germline": Path(sample["germline_phased_vcf"]),
        "normal_bam": Path(sample["normal_bam"]),
        "tumor_bam": Path(sample["tumor_bam"]),
        "reference": Path(sample["reference"]),
        "caller_pass": Path(sample["caller_pass_vcf"]),
    }
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    observed = {row["role"]: row for row in rows}
    if set(observed) != set(expected):
        raise RuntimeError(f"{sample['sample']}: input inventory roles differ")
    for role, artifact in expected.items():
        row = observed[role]
        if Path(row["path"]).resolve() != artifact.resolve():
            raise RuntimeError(f"{sample['sample']}: {role} path drift")
        if int(row["size_bytes"]) != artifact.lstat().st_size:
            raise RuntimeError(f"{sample['sample']}: {role} size drift")
        if row["mtime"] != stat_mtime(artifact):
            raise RuntimeError(f"{sample['sample']}: {role} mtime drift")


def validate_indexed_sidecar(path):
    index = Path(f"{path}.tbi")
    if not path.is_file() or not index.is_file():
        raise RuntimeError(f"sidecar/index missing: {path}")
    with pysam.TabixFile(str(path)) as tabix:
        if not tabix.contigs:
            raise RuntimeError(f"sidecar has no indexed contigs: {path}")
        first = next(tabix.fetch(tabix.contigs[0]), None)
        if first is None or len(first.split("\t")) != 9:
            raise RuntimeError(f"sidecar first indexed row is malformed: {path}")


def validate_status(path, expected_samples):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    failures = [row for row in rows if row["status"] == "FAIL"]
    passed = [row["sample"] for row in rows
              if row["stage"] == "production_tagging" and row["status"] == "PASS"]
    all_pass = [row for row in rows
                if row["sample"] == "ALL" and row["stage"] == "verify" and row["status"] == "PASS"]
    if failures or Counter(passed) != Counter(expected_samples) or len(all_pass) != 1:
        raise RuntimeError("run_status.tsv does not prove exactly one PASS per dataset plus ALL verification")
    return len(rows)


def copy_source_bundle(staging, code_paths, binary, closeout_script):
    bundle = staging / "source_bundle"
    bundle.mkdir()
    source_map = []
    inputs = [*code_paths, binary.resolve(), closeout_script.resolve()]
    for index, source in enumerate(inputs, start=1):
        target = bundle / f"{index:02d}_{source.name}"
        shutil.copy2(source, target)
        source_map.append({"source": str(source), "bundle": target.name, "sha256": sha256(target)})
    write_json(bundle / "source_map.json", source_map)
    return bundle, source_map


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--expected-manifest", required=True, type=Path)
    parser.add_argument("--longphase-binary", required=True, type=Path)
    parser.add_argument("--expected-binary-sha256", required=True)
    args = parser.parse_args()
    root = args.run_root.resolve()
    if not root.is_dir():
        raise SystemExit(f"run root missing: {root}")
    if (root / "_SUCCESS").exists() or (root / "closeout").exists():
        raise SystemExit("immutable closeout already exists")
    required_root = ["input_manifest.json", "verification_summary.json", "run_status.tsv", "params.json",
                     "code.sha256", "runtime_executable_receipt.json"]
    missing = [name for name in required_root if not (root / name).is_file()]
    if missing:
        raise SystemExit(f"production run is incomplete: {missing}")

    expected = load_json(args.expected_manifest)
    frozen = load_json(root / "input_manifest.json")
    verification = load_json(root / "verification_summary.json")
    runtime_binary = load_json(root / "runtime_executable_receipt.json")
    samples = [item["sample"] for item in expected.get("samples", [])]
    if expected != frozen:
        raise SystemExit("frozen production manifest differs from the expected manifest")
    if len(samples) != 7 or len(set(samples)) != 7:
        raise SystemExit("production manifest must contain seven unique datasets")
    if (verification.get("dataset_count") != 7 or verification.get("n_pass") != 7
            or verification.get("all_pass") is not True or len(verification.get("samples", [])) != 7):
        raise SystemExit("production verification is not 7/7 all-pass")
    if sha256(args.longphase_binary) != args.expected_binary_sha256:
        raise SystemExit("current LongPhase-S binary hash differs from the expected hash")
    if (runtime_binary.get("binary", {}).get("sha256") != args.expected_binary_sha256
            or runtime_binary.get("all_active_proc_executables_match_path") is not True):
        raise SystemExit("active-process LongPhase-S executable receipt does not match")

    code_hashes = read_hash_manifest(root / "code.sha256")
    status_rows = validate_status(root / "run_status.tsv", samples)
    verified_by_sample = {Path(item.get("sidecar", "")).parent.name: item
                          for item in verification["samples"]}
    if set(verified_by_sample) != set(samples):
        raise SystemExit("verification sample set differs from the manifest")

    all_hashes = dict(code_hashes)
    sample_receipts = []
    for sample in expected["samples"]:
        name = sample["sample"]
        wd = root / "samples" / name
        if not wd.is_dir():
            raise SystemExit(f"sample output directory missing: {name}")
        expected_names = {
            "validation": wd / "sidecar_validation.json",
            "capture": wd / "stream_capture_summary.json",
            "sidecar_plain": wd / f"{name}.read_tags.tsv",
            "sidecar_bgzf": wd / f"{name}.read_tags.tsv.gz",
            "vcf_all": wd / f"{name}.longphase_s.recalibrated.all.vcf.gz",
            "vcf_pass": wd / f"{name}.longphase_s.recalibrated.pass.vcf.gz",
            "vcf_raw": wd / f"{name}_production_sc.vcf",
            "purity": wd / f"{name}_production_purity.out",
            "fifo": wd / "consumed_tagged_bam.fifo",
        }
        absent = [str(path) for key, path in expected_names.items() if key != "fifo" and not path.is_file()]
        if absent:
            raise SystemExit(f"{name}: required production outputs missing: {absent}")
        if not expected_names["fifo"].exists() or not stat.S_ISFIFO(expected_names["fifo"].stat().st_mode):
            raise SystemExit(f"{name}: consumed BAM stream is not preserved as a FIFO receipt")
        regular_bams = [path for path in wd.glob("*.bam") if path.is_file()]
        if regular_bams:
            raise SystemExit(f"{name}: unexpected persisted tagged BAM payload: {regular_bams}")

        validate_inventory(wd / "input_files.tsv", sample)
        merge_hashes(all_hashes, read_hash_manifest(wd / "input.sha256"))
        merge_hashes(all_hashes, read_hash_manifest(wd / "output.sha256"))
        validation = load_json(expected_names["validation"])
        capture = load_json(expected_names["capture"])
        checks = validation.get("checks", {})
        if (validation.get("pass") is not True or validation.get("region") != "all"
                or not checks or not all(checks.values())
                or validation.get("duplicate_exact_alignment_rows") != 0
                or validation.get("duplicate_exact_alignment_conflicts") != 0
                or validation.get("record_key_missing") != 0 or validation.get("record_key_extra") != 0
                or set(validation.get("HP_counts", {})) - EXPECTED_HP
                or capture.get("pass") is not True):
            raise SystemExit(f"{name}: strict sidecar validation contract failed")
        if Path(validation["sidecar"]).resolve() != expected_names["sidecar_plain"].resolve():
            raise SystemExit(f"{name}: verification points to an unexpected sidecar")
        validate_indexed_sidecar(expected_names["sidecar_bgzf"])

        input_keys, _ = variant_keys(Path(sample["caller_pass_vcf"]))
        all_keys, all_pass_keys = variant_keys(expected_names["vcf_all"])
        pass_keys, pass_only_keys = variant_keys(expected_names["vcf_pass"])
        if (input_keys != all_keys or all_pass_keys != pass_keys or pass_keys != pass_only_keys
                or sum(input_keys.values()) != int(sample["caller_pass_records"])):
            raise SystemExit(f"{name}: ClairS-input/LongPhase-all/PASS VCF key contract failed")
        for vcf_path in (expected_names["vcf_all"], expected_names["vcf_pass"]):
            index = Path(f"{vcf_path}.csi")
            if not index.is_file():
                raise SystemExit(f"{name}: VCF CSI missing: {index}")
            with pysam.VariantFile(str(vcf_path)) as vcf:
                if not vcf.header.contigs:
                    raise SystemExit(f"{name}: indexed VCF header has no contigs")
                first_contig = next(iter(vcf.header.contigs))
                next(vcf.fetch(first_contig), None)

        sample_receipts.append({
            "sample": name,
            "input_records": sum(input_keys.values()),
            "longphase_all_records": sum(all_keys.values()),
            "longphase_pass_records": sum(pass_keys.values()),
            "sidecar_rows": capture["rows_mapped"],
            "HP_counts": validation["HP_counts"],
            "HP_with_PS_counts": validation["HP_with_PS_counts"],
            "record_key_missing": 0,
            "record_key_extra": 0,
            "duplicate_exact_alignment_rows": 0,
            "duplicate_exact_alignment_conflicts": 0,
        })

    staging = root / f".closeout.pending.{os.getpid()}"
    if staging.exists():
        raise SystemExit(f"closeout staging path exists: {staging}")
    staging.mkdir()
    code_paths = list(code_hashes)
    _, source_map = copy_source_bundle(staging, code_paths, args.longphase_binary, Path(__file__))
    environment = {
        "python": sys.version.splitlines()[0],
        "pysam": pysam.__version__,
        "samtools": command_version(["samtools", "--version"]),
        "bcftools": command_version(["bcftools", "--version"]),
        "bgzip": command_version(["bgzip", "--version"]),
        "tabix": command_version(["tabix", "--version"]),
    }
    receipt = {
        "schema_version": "1.0",
        "status": "PASS",
        "task_type": "B_COMPREHENSIVE_VALIDATION_PRODUCER",
        "dataset_count": 7,
        "n_pass": 7,
        "truth_flags": False,
        "tree_backbone": "LongPhase-S _sc.vcf FILTER=PASS",
        "binary_sha256": args.expected_binary_sha256,
        "runtime_binary_receipt": str(root / "runtime_executable_receipt.json"),
        "status_rows": status_rows,
        "source_bundle": source_map,
        "environment": environment,
        "samples": sample_receipts,
    }
    write_json(staging / "production_closeout.json", receipt)

    critical = [root / name for name in required_root]
    critical.extend(path for wd in (root / "samples").iterdir() if wd.is_dir()
                    for path in wd.iterdir() if path.is_file())
    critical.extend(path for path in (staging / "source_bundle").iterdir() if path.is_file())
    critical.append(staging / "production_closeout.json")
    final_hashes = dict(all_hashes)
    closeout = root / "closeout"
    for artifact in critical:
        if artifact.is_relative_to(staging):
            recorded_path = (closeout / artifact.relative_to(staging)).resolve()
        else:
            recorded_path = artifact.resolve()
        final_hashes.setdefault(recorded_path, sha256(artifact))
    with (staging / "artifacts.final.sha256").open("w", encoding="utf-8") as handle:
        for artifact, digest in sorted(final_hashes.items(), key=lambda item: str(item[0])):
            handle.write(f"{digest}  {artifact}\n")

    os.replace(staging, closeout)
    closeout_receipt = closeout / "production_closeout.json"
    success = {
        "schema_version": "1.0",
        "status": "SUCCESS",
        "closeout_receipt": str(closeout_receipt),
        "closeout_receipt_sha256": sha256(closeout_receipt),
        "artifacts_manifest": str(closeout / "artifacts.final.sha256"),
        "artifacts_manifest_sha256": sha256(closeout / "artifacts.final.sha256"),
    }
    pending = root / f"_SUCCESS.pending.{os.getpid()}"
    write_json(pending, success)
    os.replace(pending, root / "_SUCCESS")
    print(f"LONGPHASE-S PRODUCTION CLOSEOUT: 7/7 PASS -> {root / '_SUCCESS'}")


if __name__ == "__main__":
    main()
