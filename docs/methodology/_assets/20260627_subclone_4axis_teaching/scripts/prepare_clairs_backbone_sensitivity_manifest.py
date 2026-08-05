#!/usr/bin/env python3
"""Derive the legacy v2.1 ClairS tree sensitivity manifest.

For raw-all layered v3, rebuild from production with
``prepare_clean_layered_manifest_v3.py --tree-input-contract clairs_FILTER_PASS_sensitivity``
so all six VCF roles and the producer-bound LongPhase-S PASS artifact remain explicit.
"""

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--clean-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"immutable sensitivity manifest exists: {args.output}")
    canonical = json.loads(args.clean_manifest.read_text(encoding="utf-8"))
    if (canonical.get("schema_version") != "2.1"
            or canonical.get("tree_input_contract") != "longphase_recalibrated_PASS"
            or canonical.get("tag_contract", {}).get("tree_backbone") != "LongPhase-S _sc.vcf FILTER=PASS"
            or canonical.get("dataset_count") != 7):
        raise SystemExit("source is not a full clean canonical LongPhase-S manifest")
    samples = []
    for source in canonical["samples"]:
        item = dict(source)
        item.update({
            "somatic_vcf": source["longphase_input_vcf"],
            "somatic_vcf_role": "ClairS paired FILTER=PASS tree-input sensitivity; not canonical",
            "backbone_sensitivity_label": "clairs_PASS_input",
            "canonical_tree_vcf": source["longphase_recalibrated_pass_vcf"],
        })
        if not Path(item["somatic_vcf"]).is_file():
            raise SystemExit(f"{item['sample']} ClairS sensitivity input missing: {item['somatic_vcf']}")
        samples.append(item)
    tag_contract = dict(canonical["tag_contract"])
    tag_contract["tree_backbone"] = "ClairS PASS sensitivity"
    output = dict(canonical)
    output.update({
        "task_type": "B_BACKBONE_SENSITIVITY",
        "manifest_id": "20260711_clairs_PASS_tree_backbone_sensitivity_7datasets",
        "analysis_scope": "chr1-22 x 7 datasets; ClairS PASS tree input using the same clean LongPhase-S HP/PS tags",
        "tree_input_contract": "clairs_PASS_input",
        "scope_flag": "FULL_7_DATASET_SENSITIVITY_NOT_CANONICAL",
        "canonical_manifest": str(args.clean_manifest.resolve()),
        "canonical_manifest_sha256": sha256(args.clean_manifest),
        "tag_contract": tag_contract,
        "samples": samples,
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"CLAIRS BACKBONE SENSITIVITY MANIFEST: 7/7 -> {args.output}")


if __name__ == "__main__":
    main()
