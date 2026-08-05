#!/usr/bin/env python3
"""Build the strict normalized-raw-all LongPhase-S production manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import shlex
from collections import Counter
from pathlib import Path

import pysam


EXPECTED_DATASETS = (
    "HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"
)
EXPECTED_BINARY_SHA256 = "5ceba723d31c52f01202478de19952219371f0b7f9c136b8882cdd93026789b8"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def index_path(path: Path) -> Path:
    for suffix in (".csi", ".tbi"):
        candidate = Path(f"{path}{suffix}")
        if candidate.is_file():
            return candidate
    raise RuntimeError(f"indexed VCF has no CSI/TBI: {path}")


def option(args: list[str], name: str) -> str | None:
    if name not in args:
        return None
    index = args.index(name)
    return args[index + 1] if index + 1 < len(args) else None


def longphase_args(bam_path: Path) -> list[str]:
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        programs = bam.header.to_dict().get("PG", [])
    commands = [item.get("CL", "") for item in programs if "somatic_haplotag" in item.get("CL", "")]
    if not commands:
        raise RuntimeError(f"LongPhase-S @PG command missing: {bam_path}")
    return shlex.split(commands[-1])


def scan_raw(path: Path) -> dict[str, object]:
    filters: Counter[str] = Counter()
    genotypes: Counter[str] = Counter()
    keys: Counter[tuple[str, int, str, tuple[str, ...]]] = Counter()
    records = biallelic = 0
    with pysam.VariantFile(str(path)) as vcf:
        if len(vcf.header.samples) != 1:
            raise RuntimeError(f"raw ClairS VCF must have exactly one sample: {path}")
        sample = next(iter(vcf.header.samples))
        for record in vcf:
            records += 1
            keys[(record.contig, int(record.pos), record.ref, tuple(record.alts or ()))] += 1
            label = ";".join(record.filter.keys()) or "."
            filters[label] += 1
            gt = record.samples[sample].get("GT")
            genotypes["/".join("." if value is None else str(value) for value in (gt or ()))] += 1
            if len(record.ref) == 1 and len(record.alts or ()) == 1 and len(record.alts[0]) == 1:
                biallelic += 1
    duplicate_excess = sum(count - 1 for count in keys.values() if count > 1)
    result = {
        "records": records,
        "unique_record_keys": len(keys),
        "duplicate_record_key_excess": duplicate_excess,
        "biallelic_snvs": biallelic,
        "filter_counts": dict(sorted(filters.items())),
        "genotype_counts": dict(sorted(genotypes.items())),
    }
    if (records == 0 or biallelic != records or duplicate_excess != 0
            or set(filters) - {"PASS", "LowQual"} or set(genotypes) - {"0/1", "1/1"}):
        raise RuntimeError(f"unsupported raw ClairS contract for {path}: {result}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layered-manifest", required=True, type=Path)
    parser.add_argument("--patched-binary", required=True, type=Path)
    parser.add_argument("--patch-receipt", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"immutable manifest exists: {args.output}")
    if sha256(args.patched_binary) != EXPECTED_BINARY_SHA256:
        raise SystemExit("patched LongPhase-S binary hash differs from reviewed probe binary")
    patch_receipt = json.loads(args.patch_receipt.read_text(encoding="utf-8"))
    if (patch_receipt.get("patched_binary_sha256") != EXPECTED_BINARY_SHA256
            or patch_receipt.get("status") != "APPROVED_FOR_FAIL_CLOSED_7_DATASET_VALIDATION"
            or patch_receipt.get("approval", {}).get("scope") != "FAIL_CLOSED_7_DATASET_VALIDATION_ONLY"):
        raise SystemExit("patch build receipt does not bind the reviewed probe binary")
    base = json.loads(args.layered_manifest.read_text(encoding="utf-8"))
    if tuple(item.get("sample") for item in base.get("samples", [])) != EXPECTED_DATASETS:
        raise SystemExit("layered source manifest is not the exact ordered seven-dataset set")
    samples = []
    for meta in base["samples"]:
        tagged_bam = Path(meta["tumor_bam"])
        longphase_dir = tagged_bam.parent
        context = json.loads((longphase_dir.parent / "run_context.json").read_text(encoding="utf-8"))
        command = longphase_args(tagged_bam)
        raw = Path(context["somatic_vcf"])
        caller_pass = Path(meta["somatic_vcf"])
        raw_scan = scan_raw(raw)
        required = {
            "caller_raw_vcf": raw,
            "caller_pass_vcf": caller_pass,
            "germline_phased_vcf": longphase_dir / "germline_phased_merged.vcf.gz",
            "normal_bam": Path(option(command, "-b") or ""),
            "tumor_bam": Path(option(command, "--tumor-bam-file") or ""),
            "reference": Path(option(command, "-r") or ""),
        }
        missing = [f"{role}:{path}" for role, path in required.items() if not path.is_file()]
        if missing:
            raise RuntimeError(f"{meta['sample']} missing production inputs: {missing}")
        samples.append({
            "sample": meta["sample"],
            "biological_id": meta["biological_id"],
            "caller": "ClairS paired",
            "longphase_input_contract": "normalized_ClairS_raw_all",
            "caller_raw_vcf": str(raw),
            "caller_raw_vcf_index": str(index_path(raw)),
            "caller_raw_scan": raw_scan,
            "caller_pass_vcf": str(caller_pass),
            "caller_pass_vcf_index": str(index_path(caller_pass)),
            "germline_phased_vcf": str(required["germline_phased_vcf"]),
            "normal_bam": str(required["normal_bam"]),
            "tumor_bam": str(required["tumor_bam"]),
            "reference": str(required["reference"]),
        })
    output = {
        "schema_version": "2.0",
        "manifest_id": "20260711_longphase_s_raw_all_production_tag_sidecar",
        "dataset_count": 7,
        "biological_sample_count": 6,
        "contract": {
            "tagging_scope": "all input contigs; no truth VCF/BED flags",
            "somatic_input_role": "normalized raw ClairS all records; bidirectional FILTER recalibration",
            "normalization": "replace FORMAT/GQ header Type=Integer with Type=Float; preserve every record/payload",
            "tree_backbone": "LongPhase-S recalibrated FILTER=PASS",
            "output_mode": "read-tag sidecar plus recalibrated VCF; no persisted BAM payload",
            "mapq": 20,
            "tag_supplementary": True,
            "read_tags": [".", "1", "2", "3", "4", "1-1", "2-1", "1-2", "2-2"],
            "zero_read_patch_semantics": "low-confidence site with no eligible tumor read remains untagged/LowQual",
        },
        "longphase_binary": {
            "path": str(args.patched_binary.resolve()),
            "sha256": EXPECTED_BINARY_SHA256,
            "patch_receipt": str(args.patch_receipt.resolve()),
            "patch_receipt_sha256": sha256(args.patch_receipt),
        },
        "samples": samples,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"RAW-ALL PRODUCTION MANIFEST: 7/7 inputs valid -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
