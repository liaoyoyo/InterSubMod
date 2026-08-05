#!/usr/bin/env python3
"""Build a production LongPhase-S manifest from locked historical inputs."""

import argparse
import json
import shlex
from pathlib import Path

import pysam


def option(args, name):
    if name not in args:
        return None
    index = args.index(name)
    return args[index + 1] if index + 1 < len(args) else None


def longphase_args(bam_path):
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        programs = bam.header.to_dict().get("PG", [])
    commands = [item.get("CL", "") for item in programs if "somatic_haplotag" in item.get("CL", "")]
    if not commands:
        raise RuntimeError(f"LongPhase-S @PG command missing: {bam_path}")
    return shlex.split(commands[-1])


def count_records(path):
    total = snvs = 0
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            total += 1
            if len(record.ref) == 1 and len(record.alts or ()) == 1 and len(record.alts[0]) == 1:
                snvs += 1
    return total, snvs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--layered-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    base = json.loads(args.layered_manifest.read_text(encoding="utf-8"))
    samples = []
    for meta in base["samples"]:
        tagged_bam = Path(meta["tumor_bam"])
        longphase_dir = tagged_bam.parent
        run_root = longphase_dir.parent
        context = json.loads((run_root / "run_context.json").read_text(encoding="utf-8"))
        command = longphase_args(tagged_bam)
        caller_pass = Path(meta["somatic_vcf"])
        total, snvs = count_records(caller_pass)
        sample = {
            "sample": meta["sample"],
            "biological_id": meta["biological_id"],
            "caller": "ClairS paired",
            "caller_raw_vcf": context["somatic_vcf"],
            "caller_pass_vcf": str(caller_pass),
            "caller_pass_records": total,
            "caller_pass_biallelic_snvs": snvs,
            "germline_phased_vcf": str(longphase_dir / "germline_phased_merged.vcf.gz"),
            "normal_bam": option(command, "-b"),
            "tumor_bam": option(command, "--tumor-bam-file"),
            "reference": option(command, "-r"),
            "historical_tagged_bam": str(tagged_bam),
            "historical_truth_vcf": option(command, "--truth-vcf"),
            "historical_truth_bed": option(command, "--truth-bed"),
        }
        required = [sample["caller_raw_vcf"], sample["caller_pass_vcf"], sample["germline_phased_vcf"],
                    sample["normal_bam"], sample["tumor_bam"], sample["reference"]]
        missing = [path for path in required if not path or not Path(path).is_file()]
        if missing:
            raise RuntimeError(f"{sample['sample']} missing production inputs: {missing}")
        samples.append(sample)
    output = {
        "schema_version": "1.0",
        "manifest_id": "20260711_longphase_s_production_tag_sidecar",
        "dataset_count": len(samples),
        "biological_sample_count": len({item["biological_id"] for item in samples}),
        "contract": {
            "tagging_scope": "genome-wide; no truth VCF/BED flags",
            "somatic_input_role": "ClairS paired FILTER=PASS input; all record keys preserved",
            "output_mode": "read-tag sidecar plus recalibrated VCF; BAM payload discarded to /dev/null",
            "read_tags": [".", "1", "2", "3", "4", "1-1", "2-1", "1-2", "2-2"],
            "mapq": 20,
            "tag_supplementary": True,
        },
        "samples": samples,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"PRODUCTION MANIFEST: {len(samples)}/{len(samples)} inputs valid -> {args.output}")


if __name__ == "__main__":
    main()
