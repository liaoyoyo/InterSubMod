#!/usr/bin/env python3
"""Audit that every tree and focal analysis uses the same latest LongPhase-S PASS run."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


TOPIC_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_MANIFEST = TOPIC_ROOT / "results" / "all_ssnv_input_manifest.json"


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_file(path: Path) -> Path:
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(path)
    return path


def identity(record: dict[str, Any]) -> tuple[str, int, str]:
    path = require_file(Path(record["path"]).resolve())
    expected = record.get("identity", record).get("sha256")
    if not isinstance(expected, str) or len(expected) != 64:
        raise RuntimeError(f"Missing SHA-256 identity for {path}")
    observed = sha256(path)
    return str(path), path.stat().st_size, observed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-manifest", type=Path, default=DEFAULT_INPUT_MANIFEST)
    parser.add_argument(
        "--output",
        type=Path,
        default=TOPIC_ROOT / "results" / "latest_tree_input_contract_audit.json",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite existing audit: {args.output}")

    input_manifest_path = require_file(args.input_manifest.resolve())
    focal = json.loads(input_manifest_path.read_text(encoding="utf-8"))
    layered_record = focal.get("layered_root") or {}
    layered_path = require_file(Path(layered_record["path"]).resolve())
    layered_sha = sha256(layered_path)
    layered = json.loads(layered_path.read_text(encoding="utf-8"))

    focal_rows = {row["sample"]: row for row in focal.get("samples", [])}
    layered_rows = {row["sample"]: row for row in layered.get("samples", [])}
    expected_samples = [row["sample"] for row in focal.get("samples", [])]
    sample_rows: list[dict[str, Any]] = []

    for sample in expected_samples:
        focal_row = focal_rows[sample]
        layered_row = layered_rows.get(sample)
        if layered_row is None:
            raise RuntimeError(f"Layered manifest is missing sample {sample}")
        somatic = layered_row["somatic"]
        tree = somatic["tree_vcf"]
        recalibrated_pass = somatic["longphase_recalibrated_pass_vcf"]
        tree_path, tree_size, tree_sha = identity(tree)
        pass_path, pass_size, pass_sha = identity(recalibrated_pass)
        focal_vcf_path = require_file(Path(focal_row["all_ssnv_vcf"]["path"]).resolve())
        focal_vcf_sha = sha256(focal_vcf_path)

        sidecar_path = str(Path(focal_row["latest_read_tag_sidecar"]["path"]).resolve())
        layered_sidecar_path = str(Path(layered_row["read_tags"]["sidecar"]["path"]).resolve())
        raw_bam_path = str(Path(focal_row["raw_alignment"]["path"]).resolve())
        layered_bam_path = str(Path(layered_row["alignment_payload"]["path"]).resolve())
        checks = {
            "backbone_role_is_latest_longphase_s_pass": (
                somatic.get("backbone_role") == "longphase_s_recalibrated_filter_pass"
            ),
            "tree_and_recalibrated_pass_paths_equal": tree_path == pass_path,
            "tree_and_recalibrated_pass_sizes_equal": tree_size == pass_size,
            "tree_and_recalibrated_pass_recorded_and_observed_sha_equal": (
                tree_sha == pass_sha
                and tree_sha == tree["identity"]["sha256"]
                and pass_sha == recalibrated_pass["identity"]["sha256"]
            ),
            "focal_subset_vcf_hash_matches_frozen_manifest": (
                focal_vcf_sha == focal_row["all_ssnv_vcf"]["sha256"]
            ),
            "focal_manifest_checks_all_pass": all(focal_row.get("checks", {}).values()),
            "raw_bam_path_matches_layered_manifest": raw_bam_path == layered_bam_path,
            "latest_sidecar_path_matches_layered_manifest": sidecar_path == layered_sidecar_path,
            "embedded_bam_hp_ps_ignored": (
                layered_row["alignment_payload"].get("embedded_hp_ps_policy") == "ignore"
            ),
            "external_sidecar_is_tag_authority": (
                layered_row["read_tags"].get("authority") == "external_coordinate_sidecar"
                and layered_row["read_tags"].get("fallback_policy") == "forbidden"
            ),
        }
        sample_rows.append(
            {
                "sample": sample,
                "biological_id": focal_row["biological_id"],
                "all_ssnv_count": focal_row["counts"]["all_ssnv"],
                "tree_vcf": tree_path,
                "tree_vcf_sha256": tree_sha,
                "all_ssnv_focal_subset_vcf": str(focal_vcf_path),
                "all_ssnv_focal_subset_vcf_sha256": focal_vcf_sha,
                "raw_bam": raw_bam_path,
                "latest_hp_ps_sidecar": sidecar_path,
                "checks": checks,
                "pass": all(checks.values()),
            }
        )

    top_checks = {
        "input_manifest_pass": focal.get("pass") is True,
        "layered_manifest_sha_matches_frozen_record": layered_sha == layered_record.get("sha256"),
        "same_seven_datasets": (
            len(expected_samples) == 7 and set(expected_samples) == set(layered_rows)
        ),
        "producer_seven_of_seven_pass": (
            layered.get("production_summary", {}).get("all_pass") is True
            and layered.get("production_summary", {}).get("passed_dataset_count") == 7
        ),
        "tree_input_contract_is_latest_longphase_s_pass": (
            layered.get("analysis_contract", {}).get("tree_input_contract")
            == "longphase_s_recalibrated_FILTER_PASS"
        ),
        "read_tag_contract_requires_exact_external_join": (
            layered.get("analysis_contract", {}).get("read_tag_mode") == "external_sidecar"
            and layered.get("analysis_contract", {}).get("require_exact_join") is True
            and layered.get("analysis_contract", {}).get("embedded_tag_policy") == "ignore"
        ),
        "focal_scope_is_469849": focal.get("totals", {}).get("all_ssnv") == 469849,
    }
    payload = {
        "schema_name": "intersubmod.latest_tree_input_contract_audit",
        "schema_version": "1.0.0",
        "created_at_utc": now_utc(),
        "input_manifest": {
            "path": str(input_manifest_path),
            "sha256": sha256(input_manifest_path),
        },
        "layered_manifest": {"path": str(layered_path), "sha256": layered_sha},
        "scope": focal.get("scope"),
        "totals": focal.get("totals"),
        "top_level_checks": top_checks,
        "samples": sample_rows,
        "interpretation": (
            "Every layered tree consumes the same-run LongPhase-S recalibrated FILTER=PASS VCF. "
            "The 469,849 focal universe is its frozen chr1-22 biallelic-sSNV PASS subset. Raw BAM "
            "embedded HP/PS is ignored; latest HP/PS is supplied only by the same-run sidecar."
        ),
        "limitations": [
            "This audit proves file and manifest identity, not that every LongPhase-S PASS call is a true somatic mutation.",
            "The full all-site read-level HP/PS terminal join is audited by the downstream 469,849-site screen receipt.",
        ],
    }
    payload["pass"] = all(top_checks.values()) and all(row["pass"] for row in sample_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "samples": len(sample_rows),
                "all_ssnv": focal["totals"]["all_ssnv"],
                "pass": payload["pass"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
