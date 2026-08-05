#!/usr/bin/env python3
"""Prepare a layered-v2 manifest backed by production LongPhase-S tag sidecars."""

import argparse
import hashlib
import json
from pathlib import Path


def path_fingerprint(path):
    path = Path(path)
    logical = path.lstat()
    resolved = path.resolve()
    target = resolved.stat()
    return {
        "logical_path": str(path),
        "logical_size_bytes": logical.st_size,
        "logical_mtime_ns": logical.st_mtime_ns,
        "resolved_path": str(resolved),
        "resolved_size_bytes": target.st_size,
        "resolved_mtime_ns": target.st_mtime_ns,
    }


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_hash_manifest(path):
    values = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, raw_path = line.split(None, 1)
        artifact = Path(raw_path.strip()).resolve()
        if artifact in values and values[artifact] != digest:
            raise RuntimeError(f"conflicting production artifact hashes: {artifact}")
        values[artifact] = digest
    return values


def verify_recorded_artifact(path, hashes):
    artifact = Path(path).resolve()
    expected = hashes.get(artifact)
    if expected is None:
        raise RuntimeError(f"production artifact is absent from closeout hash manifest: {artifact}")
    if not artifact.is_file() or sha256(artifact) != expected:
        raise RuntimeError(f"production artifact changed after closeout: {artifact}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-manifest", required=True, type=Path)
    parser.add_argument("--longphase-manifest", required=True, type=Path)
    parser.add_argument("--production-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    base = json.loads(args.base_manifest.read_text(encoding="utf-8"))
    lps = json.loads(args.longphase_manifest.read_text(encoding="utf-8"))
    production_manifest_path = args.production_root / "input_manifest.json"
    verification_path = args.production_root / "verification_summary.json"
    success_path = args.production_root / "_SUCCESS"
    closeout_path = args.production_root / "closeout" / "production_closeout.json"
    if (not production_manifest_path.is_file() or not verification_path.is_file()
            or not success_path.is_file() or not closeout_path.is_file()):
        raise RuntimeError("production root is incomplete: frozen manifest, verification, closeout, or _SUCCESS missing")
    production_manifest = json.loads(production_manifest_path.read_text(encoding="utf-8"))
    verification = json.loads(verification_path.read_text(encoding="utf-8"))
    closeout = json.loads(closeout_path.read_text(encoding="utf-8"))
    success = json.loads(success_path.read_text(encoding="utf-8"))
    if production_manifest != lps:
        raise RuntimeError("production frozen input manifest differs from requested LongPhase manifest")
    expected_samples = {item["sample"] for item in base["samples"]}
    verified_samples = {
        Path(item.get("sidecar", "")).parent.name for item in verification.get("samples", [])
    }
    if ({item["sample"] for item in lps["samples"]} != expected_samples
            or verified_samples != expected_samples):
        raise RuntimeError("base/LongPhase/verification sample sets differ")
    if (verification.get("dataset_count") != len(expected_samples)
            or verification.get("n_pass") != len(expected_samples)
            or verification.get("all_pass") is not True):
        raise RuntimeError("production LongPhase root is not all-pass")
    if (closeout.get("status") != "PASS" or closeout.get("dataset_count") != len(expected_samples)
            or closeout.get("n_pass") != len(expected_samples)
            or closeout.get("truth_flags") is not False
            or closeout.get("tree_backbone") != "LongPhase-S _sc.vcf FILTER=PASS"):
        raise RuntimeError("production closeout does not satisfy the strict seven-dataset contract")
    artifacts_manifest = args.production_root / "closeout" / "artifacts.final.sha256"
    if (success.get("status") != "SUCCESS"
            or Path(success.get("closeout_receipt", "")).resolve() != closeout_path.resolve()
            or success.get("closeout_receipt_sha256") != sha256(closeout_path)
            or Path(success.get("artifacts_manifest", "")).resolve() != artifacts_manifest.resolve()
            or not artifacts_manifest.is_file()
            or success.get("artifacts_manifest_sha256") != sha256(artifacts_manifest)):
        raise RuntimeError("production _SUCCESS receipt does not match the closeout artifacts")
    production_hashes = load_hash_manifest(artifacts_manifest)
    lps_by_sample = {item["sample"]: item for item in lps["samples"]}
    samples = []
    for original in base["samples"]:
        sample = original["sample"]
        source = lps_by_sample[sample]
        sample_dir = args.production_root / "samples" / sample
        item = dict(original)
        item.update({
            "tumor_bam": source["tumor_bam"],
            "tumor_bam_role": "unaltered aligned ONT BAM; HP/PS joined from exact alignment sidecar",
            "normal_bam": source["normal_bam"],
            "tumor_bam_fingerprint": path_fingerprint(source["tumor_bam"]),
            "caller_raw_vcf": source["caller_raw_vcf"],
            "longphase_input_vcf": source["caller_pass_vcf"],
            "longphase_input_vcf_role": "ClairS paired FILTER=PASS supplied to LongPhase-S",
            "read_tag_sidecar": str(sample_dir / f"{sample}.read_tags.tsv.gz"),
            "read_tag_sidecar_index": str(sample_dir / f"{sample}.read_tags.tsv.gz.tbi"),
            "read_tag_validation": str(sample_dir / "sidecar_validation.json"),
            "longphase_input_inventory": str(sample_dir / "input_files.tsv"),
            "longphase_production_closeout": str(closeout_path),
            "longphase_recalibrated_all_vcf": str(sample_dir / f"{sample}.longphase_s.recalibrated.all.vcf.gz"),
            "longphase_recalibrated_pass_vcf": str(sample_dir / f"{sample}.longphase_s.recalibrated.pass.vcf.gz"),
            "somatic_vcf": str(sample_dir / f"{sample}.longphase_s.recalibrated.pass.vcf.gz"),
            "somatic_vcf_role": "LongPhase-S recalibrated FILTER=PASS; primary sSNV tree universe",
            "longphase_tagging_scope": "production genome-wide; no truth VCF/BED flags",
            "historical_truth_restricted_tagged_bam": original["tumor_bam"],
        })
        required = [item["tumor_bam"], f"{item['tumor_bam']}.bai", item["caller_raw_vcf"],
                    item["longphase_input_vcf"], item["somatic_vcf"], item["read_tag_sidecar"],
                    item["read_tag_sidecar_index"],
                    item["read_tag_validation"], item["longphase_input_inventory"],
                    item["longphase_recalibrated_all_vcf"],
                    item["longphase_recalibrated_pass_vcf"]]
        missing = [path for path in required if not Path(path).is_file()]
        if missing:
            raise RuntimeError(f"{sample} clean layered inputs missing: {missing}")
        production_outputs = [
            item["read_tag_sidecar"], item["read_tag_sidecar_index"], item["read_tag_validation"],
            item["longphase_input_inventory"], item["longphase_recalibrated_all_vcf"],
            f"{item['longphase_recalibrated_all_vcf']}.csi", item["longphase_recalibrated_pass_vcf"],
            f"{item['longphase_recalibrated_pass_vcf']}.csi", item["longphase_input_vcf"],
        ]
        for artifact in production_outputs:
            verify_recorded_artifact(artifact, production_hashes)
        validation = json.loads(Path(item["read_tag_validation"]).read_text(encoding="utf-8"))
        if (validation.get("pass") is not True or validation.get("region") != "all"
                or validation.get("duplicate_exact_alignment_rows") != 0
                or validation.get("duplicate_exact_alignment_conflicts") != 0):
            raise RuntimeError(f"{sample} production tag validation did not pass strict all-genome identity gates")
        samples.append(item)
    output = {
        "schema_version": "2.1",
        "task_type": "B_COMPREHENSIVE_VALIDATION",
        "manifest_id": "20260711_layered_v2_clean_production_tags",
        "analysis_scope": "all raw ClairS records ledger; chr1-22 LongPhase-S recalibrated PASS sSNV reconstruction",
        "tree_input_contract": "longphase_recalibrated_PASS",
        "dataset_count": len(samples),
        "biological_sample_count": len({item["biological_id"] for item in samples}),
        "tag_contract": {
            "source": "LongPhase-S production tag sidecar",
            "tree_backbone": "LongPhase-S _sc.vcf FILTER=PASS",
            "longphase_filtering_policy": "production_default_filter",
            "truth_flags": False,
            "alignment_identity": "QNAME+chrom+start+end+FLAG+BLAKE2b8(CIGAR)",
            "HP_vocabulary": [".", "1", "2", "3", "4", "1-1", "2-1", "1-2", "2-2"],
            "PS_preserved": True,
            "bam_identity_locked": True,
            "production_closeout": str(closeout_path),
            "production_closeout_sha256": sha256(closeout_path),
            "production_success": str(success_path),
            "production_success_sha256": sha256(success_path),
            "production_artifacts_manifest": str(artifacts_manifest),
            "production_artifacts_manifest_sha256": sha256(artifacts_manifest),
        },
        "samples": samples,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"CLEAN LAYERED MANIFEST: {len(samples)}/{len(samples)} valid -> {args.output}")


if __name__ == "__main__":
    main()
