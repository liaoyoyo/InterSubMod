#!/usr/bin/env python3
"""Prepare legacy one-dataset backbone probes; excluded from the authoritative validation."""

import argparse
import json
import subprocess
from pathlib import Path

import pysam


LPS_VCF = Path("/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged_sc.vcf")
HC_BED = Path("/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed")
DEEPSOMATIC_VCF = Path("/big8_disk/data/HCC1395/ONT_Dorado/DeepSomatic_v1_8_0/deepsomatic_HCC1395.vcf.gz")


def run(command, log):
    proc = subprocess.run(command, text=True, capture_output=True, check=False)
    log.append({"command": command, "exit_code": proc.returncode,
                "stdout": proc.stdout[-4000:], "stderr": proc.stderr[-4000:]})
    if proc.returncode:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(command)}\n{proc.stderr}")


def count_snvs(path):
    counts = {"total": 0, "autosomal": 0, "out_of_scope": 0}
    autosomes = {f"chr{i}" for i in range(1, 23)}
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            if len(record.ref) != 1 or len(record.alts or []) != 1 or len(record.alts[0]) != 1:
                continue
            counts["total"] += 1
            counts["autosomal" if record.contig in autosomes else "out_of_scope"] += 1
    return counts


def one_sample_manifest(base, sample, vcf, label, source):
    meta = next(dict(item) for item in base["samples"] if item["sample"] == sample)
    meta["somatic_vcf"] = str(vcf.resolve())
    meta["backbone_sensitivity_label"] = label
    meta["backbone_source"] = source
    return {
        "schema_version": "2.0",
        "manifest_id": f"20260710_backbone_sensitivity_{label}",
        "analysis_scope": "chr1-22",
        "dataset_count": 1,
        "biological_sample_count": 1,
        "backbone_contract": "Sensitivity input only; not the primary operational universe",
        "samples": [meta],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-manifest", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--allow-legacy-partial", action="store_true",
                        help="Acknowledge that these one-dataset probes cannot satisfy the current full-scope gate")
    args = parser.parse_args()
    if not args.allow_legacy_partial:
        raise SystemExit(
            "legacy partial workflow is not authoritative; use "
            "prepare_clairs_backbone_sensitivity_manifest.py for the full seven-dataset comparison"
        )
    if args.output_root.exists():
        raise SystemExit(f"immutable output root exists: {args.output_root}")
    args.output_root.mkdir(parents=True)
    base = json.loads(args.base_manifest.read_text(encoding="utf-8"))
    base_hcc = Path(next(x["somatic_vcf"] for x in base["samples"] if x["sample"] == "HCC1395"))
    jobs = [
        ("longphase_s_pass", "HCC1395", LPS_VCF, [], "LongPhase-S recalibrated FILTER=PASS"),
        ("clairs_pass_in_seqc2_hc", "HCC1395", base_hcc, ["-R", str(HC_BED)], "ClairS paired PASS restricted to SEQC2 HC regions"),
        ("deepsomatic_pass", "HCC1395_DORADO", DEEPSOMATIC_VCF, [], "DeepSomatic v1.8.0 FILTER=PASS"),
    ]
    provenance = []
    for label, sample, source_vcf, extra, description in jobs:
        for required in (source_vcf, *([HC_BED] if extra else [])):
            if not required.is_file():
                raise SystemExit(f"missing sensitivity input: {required}")
        job_dir = args.output_root / label
        job_dir.mkdir()
        output_vcf = job_dir / f"{label}.vcf.gz"
        logs = []
        command = ["bcftools", "view", "-f", "PASS", "-m2", "-M2", "-v", "snps", *extra,
                   "-Oz", "-o", str(output_vcf), str(source_vcf)]
        run(command, logs)
        run(["bcftools", "index", "-c", str(output_vcf)], logs)
        manifest = one_sample_manifest(base, sample, output_vcf, label, description)
        manifest_path = job_dir / "input_manifest.json"
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        counts = count_snvs(output_vcf)
        (job_dir / "prepare_commands.json").write_text(json.dumps(logs, indent=2) + "\n", encoding="utf-8")
        provenance.append({"label": label, "sample": sample, "source": str(source_vcf),
                           "output_vcf": str(output_vcf), "manifest": str(manifest_path), "counts": counts})
        print(f"{label}: {counts} -> {output_vcf}")
    (args.output_root / "preparation_summary.json").write_text(
        json.dumps({"schema_version": "2.0", "jobs": provenance}, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8")


if __name__ == "__main__":
    main()
