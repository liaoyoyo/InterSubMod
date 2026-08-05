#!/usr/bin/env python3
"""Screen all frozen LongPhase-S PASS sSNVs for focal-ALT methylation multigroups.

Clustering is deliberately truth-blind. Truth labels and layered reconstruction
context are joined only after the methylation screen has completed for a site.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import multiprocessing
import os
import re
import sys
from collections import Counter, deque
from concurrent.futures import Future, ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator

import numpy as np
import pysam


TOPIC_ROOT = Path(__file__).resolve().parents[1]
FP_SCRIPT_ROOT = (
    TOPIC_ROOT.parent / "20260715_single_fp_alt_multicluster_subclone_validation" / "scripts"
)
for module_root in (Path(__file__).resolve().parent, FP_SCRIPT_ROOT):
    if str(module_root) not in sys.path:
        sys.path.insert(0, str(module_root))

import focal_alt_cluster_lib as F  # noqa: E402
import latest_tag_join as TAGS  # noqa: E402


DATASETS = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]
EXPECTED_ALL_SSNV = 469_849
OUTPUT_SCHEMA_VERSION = "1.2.0"
TRUTH_LABELS = {"TP", "FP", "UNASSESSED"}
SITE_PATTERN = re.compile(r"^(chr(?:[1-9]|1[0-9]|2[0-2]))_(\d+)$")
MODAL_ASSIGNMENT_ARI_THRESHOLD = 0.8
SCREEN_CONTRACT = "phylo-v4.1_column_null95_modal_K10_RNULL40_min_group3"
M1_STABILITY_GATE_CONTRACT = (
    "coarse_ng>=2 AND not unstable AND modal_assignment_ari_min>=0.8"
)
CLUSTERING_FORBIDDEN_FIELDS = {
    "biological_id",
    "truth_label",
    "ssnv_branch",
    "ledger_selected_for_read_census",
    "component_id",
    "component_size",
    "selected_group_id",
    "selected_group_size",
    "pooled_ref_reads",
    "pooled_alt_reads",
    "layered_family_coverage",
    "layered_alt_families",
}
PRIOR_SCREEN_THRESHOLDS = {
    "minimum_group_size": 3,
    "minimum_between_within_ratio": 1.3,
    "rnull": 40,
    "modal_confidence": 0.7,
    "hidden_heterogeneity_fraction": 0.30,
    "seeds": 10,
    "coarse_null_percentile": 95.0,
    "fine_null_percentile": 90.0,
    "modal_assignment_ari_min": MODAL_ASSIGNMENT_ARI_THRESHOLD,
    "null_mode": "column",
    "empirical_alpha": None,
}

_PHYLO_PARALLEL_WORKERS = 1
_PHYLO_PARALLEL_MIN_READS = 0
_PARALLEL_PHYLO_DISTANCE: np.ndarray | None = None
_PARALLEL_PHYLO_METHYLATION: np.ndarray | None = None
STRICT_NOT_RUN_REASON = (
    "NOT_EVALUATED_AT_M1_SCREEN: R1 strict methyl-partition robustness is downstream-only and "
    "requires a G2 multi-marker molecular-haplotype base candidate; not evaluated is not failure"
)

VariantKey = tuple[str, int, str, str]
PositionKey = tuple[str, int]


CATEGORICAL_ASSOCIATIONS = ("hp_exact", "hp_family", "strand")
CONTINUOUS_ASSOCIATIONS = ("start", "end", "length", "mapq", "cpg_called")
ASSOCIATION_FIELDS = [
    field
    for prefix in CATEGORICAL_ASSOCIATIONS
    for field in (f"{prefix}_v", f"{prefix}_p_perm", f"{prefix}_n", f"{prefix}_aligned")
] + [
    field
    for prefix in CONTINUOUS_ASSOCIATIONS
    for field in (f"{prefix}_eta2", f"{prefix}_p_perm", f"{prefix}_n", f"{prefix}_aligned")
]

SITE_FIELDS = [
    "dataset",
    "sample",
    "biological_id",
    "truth_label",
    "chrom",
    "pos",
    "ref",
    "alt",
    "qual",
    "filter",
    "caller_gt",
    "caller_gq",
    "caller_dp",
    "caller_af",
    "caller_ad",
    "normal_af",
    "normal_dp",
    "normal_ad",
    "caller_original_filter",
    "longphase_filter_transition",
    "ssnv_branch",
    "ledger_selected_for_read_census",
    "component_id",
    "component_size",
    "selected_group_id",
    "selected_group_size",
    "pooled_ref_reads",
    "pooled_alt_reads",
    "layered_family_coverage",
    "layered_alt_families",
    "region_dir",
    "screen_contract",
    "m1_stability_gate_contract",
    "screen_global_fdr_calibrated",
    "screen_rnull",
    "screen_min_group_size",
    "latest_tag_join_status",
    "latest_tag_rows_fetched",
    "latest_tag_rows_eligible",
    "latest_tag_reads_joined",
    "latest_tag_ps_present",
    "latest_tag_projection_multimatch_reads",
    "source_hp_replaced_reads",
    "n_reads_total",
    "n_alt_raw",
    "n_alt_matrix",
    "n_alt_after_peel",
    "n_alt_peeled",
    "alt_readset_sha256",
    "alt_hp_counts",
    "alt_hp_family_counts",
    "alt_ps_counts",
    "alt_strand_counts",
    "phase_anchored_fraction_raw",
    "explicit_second_haplotype_fraction_raw",
    "ambiguous_tag_fraction_raw",
    "m1_evaluable",
    "m1_not_evaluable_reason",
    "analysis_status",
    "forced_k",
    "forced_silhouette",
    "coarse_ng",
    "modal_fraction",
    "fine_ng",
    "n_other",
    "n_outlier",
    "unstable",
    "hidden_heterogeneity",
    "modal_assignment_pair_count",
    "modal_assignment_ari_median",
    "modal_assignment_ari_min",
    "modal_assignment_all_pairs_ari_ge_0_8",
    "stable_null_multigroup",
    "cluster_sizes",
    "seed_group_counts",
    "hp_axis_confound",
    "technical_axis_confound",
    "residual_unexplained_multigroup",
    "phase_anchored_fraction_core",
    "phase_anchored_robust_epigenetic_candidate",
    "evidence_tier",
    "strict_methyl_partition_robustness_status",
    "strict_methyl_partition_robustness_not_evaluable_reason",
    "strict_confirm_status",
    "strict_confirm_candidate",
    "strict_confirm_candidate_is_formal_r1_claim",
    "strict_confirm_reason",
    *ASSOCIATION_FIELDS,
]


_SIDECAR_CACHE: dict[tuple[str, str], pysam.TabixFile] = {}


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def json_text(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def scalar(value: Any) -> Any:
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, np.generic):
        return value.item()
    return value


def sample_value(call: pysam.libcbcf.VariantRecordSample, key: str) -> Any:
    try:
        return scalar(call.get(key))
    except (KeyError, ValueError):
        return None


def optional_int(value: str | None) -> int | None:
    return int(value) if value not in {None, ""} else None


def open_text(path: Path) -> Any:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def variant_key(chrom: str, pos: int | str, ref: str, alt: str) -> VariantKey:
    return chrom, int(pos), ref.upper(), alt.upper()


def verify_manifest_artifact(record: dict[str, Any], require_hash: bool = False) -> Path:
    path = Path(record["path"])
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    if "size_bytes" in record and path.stat().st_size != int(record["size_bytes"]):
        raise RuntimeError(f"Frozen artifact size drift: {path}")
    if require_hash and "sha256" not in record:
        raise RuntimeError(f"Frozen artifact lacks required SHA-256: {path}")
    if "sha256" in record and sha256(path) != record["sha256"]:
        raise RuntimeError(f"Frozen artifact hash drift: {path}")
    if "sha256" not in record and "mtime_ns" in record and path.stat().st_mtime_ns != int(record["mtime_ns"]):
        raise RuntimeError(f"Frozen unhashed artifact mtime drift: {path}")
    return path


def validate_prior_screen_contract() -> None:
    observed = {
        "minimum_group_size": F.MIN_SIZE,
        "minimum_between_within_ratio": F.SEP_MIN,
        "rnull": F.RNULL,
        "modal_confidence": F.MODAL_CONFIDENCE,
        "hidden_heterogeneity_fraction": F.HIDDEN_HETEROGENEITY_FRACTION,
    }
    expected = {key: PRIOR_SCREEN_THRESHOLDS[key] for key in observed}
    if observed != expected:
        raise RuntimeError(f"Prior phylo-v4.1 threshold drift: observed={observed} expected={expected}")


def stable_null_multigroup_gate(phylo: dict[str, Any]) -> bool:
    """Require a multigroup modal solution and robust assignments across modal seeds."""
    try:
        coarse_ng = int(phylo["coarse_ng"])
        ari_min = float(phylo["modal_assignment_ari_min"])
    except (KeyError, TypeError, ValueError):
        return False
    unstable = phylo.get("unstable")
    return bool(
        coarse_ng >= 2
        and unstable is not None
        and not bool(unstable)
        and math.isfinite(ari_min)
        and ari_min >= MODAL_ASSIGNMENT_ARI_THRESHOLD
    )


def apply_m1_screen_flags(output: dict[str, Any], phylo: dict[str, Any]) -> bool:
    stable = stable_null_multigroup_gate(phylo)
    output["stable_null_multigroup"] = stable
    # Deprecated selection alias retained for existing consumers; this is not an R1 claim.
    output["strict_confirm_candidate"] = stable
    output["strict_confirm_candidate_is_formal_r1_claim"] = False
    return stable


def validate_completed_sample_run(sample_root: Path, entry: dict[str, Any]) -> dict[str, Any]:
    receipt_path = sample_root / "run_receipt.json"
    if not receipt_path.exists() or receipt_path.stat().st_size == 0:
        raise FileNotFoundError(f"Completed InterSubMod run receipt is missing: {receipt_path}")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    sample = entry["sample"]
    expected = int(entry["counts"]["all_ssnv"])
    validation = receipt.get("validation") or {}
    exact_counts = {
        "expected_vcf_sites": expected,
        "reads_files": expected,
        "bernoulli_matrix_files": expected,
        "methylation_files": expected,
    }
    for field, expected_value in exact_counts.items():
        if int(validation.get(field, -1)) != expected_value:
            raise RuntimeError(
                f"{sample} run receipt count mismatch for {field}: "
                f"{validation.get(field)!r} != {expected_value}"
            )
    summary_rows = int(validation.get("summary_rows", -1))
    if not 0 < summary_rows <= expected:
        raise RuntimeError(f"{sample} invalid run receipt summary row count: {summary_rows}")
    receipt_sidecar = (receipt.get("latest_read_tag_sidecar") or {}).get("path")
    expected_sidecar = entry["latest_read_tag_sidecar"]["path"]
    gates = {
        "schema": receipt.get("schema_name") == "intersubmod.all_ssnv_site_run",
        "sample": receipt.get("sample") == sample,
        "pass": receipt.get("pass") is True,
        "exit_code": receipt.get("exit_code") == 0,
        "validation_pass": validation.get("pass") is True,
        "vcf_sha256": receipt.get("input_vcf_sha256") == entry["all_ssnv_vcf"]["sha256"],
        "output_dir": Path(receipt.get("output_dir", "")).resolve() == sample_root.resolve(),
        "latest_sidecar": Path(receipt_sidecar or "").resolve() == Path(expected_sidecar).resolve(),
    }
    failed = [name for name, passed in gates.items() if not passed]
    if failed:
        raise RuntimeError(f"{sample} completed-run receipt failed gates: {failed}")
    return {
        "receipt": artifact(receipt_path),
        "validation": {
            key: validation.get(key)
            for key in (*exact_counts, "summary_rows", "summary_rows_not_emitted")
        },
    }


def load_truth_labels(path: Path, sample: str) -> dict[VariantKey, str]:
    labels: dict[VariantKey, str] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "chrom", "pos", "ref", "alt", "truth_label"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing truth manifest columns: {sorted(missing)}")
        for row in reader:
            if row["sample"] != sample:
                raise RuntimeError(f"Truth manifest sample mismatch: expected={sample} observed={row['sample']}")
            label = row["truth_label"]
            if label not in TRUTH_LABELS:
                raise RuntimeError(f"Unexpected truth label {label!r} in {path}")
            key = variant_key(row["chrom"], row["pos"], row["ref"], row["alt"])
            if key in labels:
                raise RuntimeError(f"Duplicate truth manifest key: {key}")
            labels[key] = label
    return labels


def ledger_row_context(row: dict[str, str]) -> dict[str, Any]:
    family_coverage = json.loads(row.get("family_coverage_json") or "{}")
    return {
        "caller_original_filter": row.get("caller_filter") or None,
        "longphase_filter_transition": row.get("longphase_filter_transition") or None,
        "ssnv_branch": row.get("ssnv_branch") or None,
        "ledger_selected_for_read_census": row.get("selected_for_read_census") or None,
        "component_id": row.get("component_id") or None,
        "component_size": optional_int(row.get("component_size")),
        "selected_group_id": row.get("selected_group_id") or None,
        "selected_group_size": optional_int(row.get("selected_group_size")),
        "pooled_ref_reads": optional_int(row.get("pooled_ref_reads")),
        "pooled_alt_reads": optional_int(row.get("pooled_alt_reads")),
        "layered_family_coverage": json_text(family_coverage),
        "layered_alt_families": json_text(
            sorted(
                family
                for family, counts in family_coverage.items()
                if isinstance(counts, list) and len(counts) >= 2 and counts[1] > 0
            )
        ),
    }


def load_ledger_context(path: Path, sample: str) -> dict[VariantKey, dict[str, Any]]:
    context: dict[VariantKey, dict[str, Any]] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "sample",
            "chrom",
            "pos",
            "ref",
            "alt",
            "longphase_recalibrated_filter",
            "ssnv_branch",
            "component_id",
            "selected_group_id",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing ledger columns: {sorted(missing)}")
        for row in reader:
            if row["longphase_recalibrated_filter"] != "PASS":
                continue
            if row["sample"] != sample:
                raise RuntimeError(f"Ledger sample mismatch: expected={sample} observed={row['sample']}")
            key = variant_key(row["chrom"], row["pos"], row["ref"], row["alt"])
            if key in context:
                raise RuntimeError(f"Duplicate PASS ledger key: {key}")
            context[key] = ledger_row_context(row)
    return context


def variant_metadata(record: pysam.libcbcf.VariantRecord) -> dict[str, Any]:
    if not record.alts or len(record.alts) != 1:
        raise RuntimeError(f"Expected one ALT allele: {record}")
    filters = list(record.filter.keys())
    if filters != ["PASS"]:
        raise RuntimeError(f"All-sSNV screen received non-PASS record: {record.chrom}:{record.pos} {filters}")
    call = next(iter(record.samples.values())) if record.samples else None
    return {
        "chrom": record.chrom,
        "pos": int(record.pos),
        "ref": record.ref.upper(),
        "alt": record.alts[0].upper(),
        "qual": record.qual,
        "filter": ";".join(filters),
        "caller_gt": json_text(sample_value(call, "GT")) if call else None,
        "caller_gq": sample_value(call, "GQ") if call else None,
        "caller_dp": sample_value(call, "DP") if call else None,
        "caller_af": sample_value(call, "AF") if call else None,
        "caller_ad": json_text(sample_value(call, "AD")) if call else None,
        "normal_af": sample_value(call, "NAF") if call else None,
        "normal_dp": sample_value(call, "NDP") if call else None,
        "normal_ad": json_text(sample_value(call, "NAD")) if call else None,
    }


def parse_site_key(reads_path: Path) -> PositionKey:
    match = SITE_PATTERN.fullmatch(reads_path.parents[2].name)
    if not match:
        raise ValueError(f"Cannot derive focal site from {reads_path}")
    return match.group(1), int(match.group(2))


def scan_site_outputs(sample_root: Path) -> dict[PositionKey, dict[str, str]]:
    if not sample_root.exists():
        raise FileNotFoundError(f"Completed InterSubMod sample output is missing: {sample_root}")
    outputs: dict[PositionKey, dict[str, str]] = {}
    for reads_path in sample_root.rglob("reads.tsv"):
        key = parse_site_key(reads_path)
        if key in outputs:
            raise RuntimeError(f"Duplicate InterSubMod output site: {sample_root.name} {key}")
        region_dir = reads_path.parents[1]
        distance_path = region_dir / "distance" / "BERNOULLI" / "matrix.csv"
        methylation_path = region_dir / "methylation" / "methylation.csv"
        for required in (reads_path, distance_path, methylation_path):
            if not required.exists() or required.stat().st_size == 0:
                raise FileNotFoundError(required)
        outputs[key] = {
            "reads_path": str(reads_path),
            "distance_path": str(distance_path),
            "methylation_path": str(methylation_path),
            "region_dir": str(region_dir),
        }
    if not outputs:
        raise RuntimeError(f"No per-site InterSubMod outputs under {sample_root}")
    return outputs


def join_posthoc_metadata(
    variant: dict[str, Any],
    truth_labels: dict[VariantKey, str],
    ledger: dict[VariantKey, dict[str, Any]],
    biological_id: str,
) -> dict[str, Any]:
    key = variant_key(variant["chrom"], variant["pos"], variant["ref"], variant["alt"])
    if key not in truth_labels:
        raise KeyError(f"Variant missing truth-manifest label: {key}")
    if key not in ledger:
        raise KeyError(f"Variant missing latest PASS ledger context: {key}")
    return {
        "biological_id": biological_id,
        "truth_label": truth_labels[key],
        **variant,
        **ledger[key],
    }


def validate_focal_key_contract(
    sample: str,
    expected: int,
    observed_keys: set[VariantKey],
    truth_labels: dict[VariantKey, str],
    ledger: dict[VariantKey, dict[str, Any]],
    observed_positions: set[PositionKey],
    outputs: dict[PositionKey, dict[str, str]],
) -> dict[str, int | bool]:
    truth_keys = set(truth_labels)
    ledger_keys = set(ledger)
    output_positions = set(outputs)
    if len(observed_keys) != expected:
        raise RuntimeError(f"{sample} VCF record count mismatch: {len(observed_keys)} != {expected}")
    if observed_keys != truth_keys:
        raise RuntimeError(
            f"{sample} VCF/truth key sets differ: "
            f"missing_truth={sorted(observed_keys - truth_keys)[:3]} "
            f"extra_truth={sorted(truth_keys - observed_keys)[:3]}"
        )
    missing_ledger = observed_keys - ledger_keys
    if missing_ledger:
        raise RuntimeError(
            f"{sample} focal VCF keys are missing from the LongPhase-S PASS ledger: "
            f"n={len(missing_ledger)} examples={sorted(missing_ledger)[:3]}"
        )
    if observed_positions != output_positions:
        raise RuntimeError(
            f"{sample} VCF/output position sets differ: "
            f"missing_outputs={sorted(observed_positions - output_positions)[:3]} "
            f"extra_outputs={sorted(output_positions - observed_positions)[:3]}"
        )
    return {
        "expected_focal_keys": expected,
        "observed_focal_keys": len(observed_keys),
        "truth_keys": len(truth_keys),
        "longphase_pass_ledger_keys": len(ledger_keys),
        "longphase_pass_ledger_extra_context_keys": len(ledger_keys - observed_keys),
        "all_focal_keys_in_longphase_pass_ledger": True,
        "output_positions": len(output_positions),
        "pass": True,
    }


def iter_sample_tasks(entry: dict[str, Any], intersubmod_root: Path) -> Iterator[dict[str, Any]]:
    sample = entry["sample"]
    expected = int(entry["counts"]["all_ssnv"])
    vcf_path = verify_manifest_artifact(entry["all_ssnv_vcf"], require_hash=True)
    truth_path = verify_manifest_artifact(entry["site_manifest"], require_hash=True)
    ledger_path = verify_manifest_artifact(entry["site_ledger"], require_hash=True)
    sidecar_path = verify_manifest_artifact(entry["latest_read_tag_sidecar"])
    sidecar_index = verify_manifest_artifact(entry["latest_read_tag_sidecar_index"])
    validate_completed_sample_run(intersubmod_root / sample, entry)
    truth_labels = load_truth_labels(truth_path, sample)
    ledger = load_ledger_context(ledger_path, sample)
    outputs = scan_site_outputs(intersubmod_root / sample)
    if len(truth_labels) != expected or len(outputs) != expected:
        raise RuntimeError(
            f"{sample} full-scope count mismatch expected={expected} truth={len(truth_labels)} "
            f"ledger_pass_context={len(ledger)} outputs={len(outputs)}"
        )

    observed_keys: set[VariantKey] = set()
    observed_positions: set[PositionKey] = set()
    with pysam.VariantFile(str(vcf_path)) as source:
        for record in source:
            metadata = variant_metadata(record)
            key = variant_key(metadata["chrom"], metadata["pos"], metadata["ref"], metadata["alt"])
            position_key = (metadata["chrom"], int(metadata["pos"]))
            if key in observed_keys or position_key in observed_positions:
                raise RuntimeError(f"Duplicate all-sSNV VCF output identity: {sample} {key}")
            if position_key not in outputs:
                raise KeyError(f"VCF site missing InterSubMod output: {sample} {position_key}")
            observed_keys.add(key)
            observed_positions.add(position_key)
            posthoc = join_posthoc_metadata(metadata, truth_labels, ledger, entry["biological_id"])
            yield {
                "screen": {
                    "sample": sample,
                    "chrom": metadata["chrom"],
                    "pos": int(metadata["pos"]),
                    "sidecar_path": str(sidecar_path),
                    "sidecar_index": str(sidecar_index),
                    **outputs[position_key],
                },
                "posthoc": posthoc,
            }
    validate_focal_key_contract(
        sample,
        expected,
        observed_keys,
        truth_labels,
        ledger,
        observed_positions,
        outputs,
    )


def get_sidecar(path: str, index: str) -> pysam.TabixFile:
    key = (path, index)
    if key not in _SIDECAR_CACHE:
        _SIDECAR_CACHE[key] = pysam.TabixFile(path, index=index)
    return _SIDECAR_CACHE[key]


def close_sidecar_cache() -> None:
    for handle in _SIDECAR_CACHE.values():
        handle.close()
    _SIDECAR_CACHE.clear()


def load_latest_joined_reads(
    reads_path: Path,
    sidecar_path: str,
    sidecar_index: str,
    chrom: str,
    pos: int,
) -> tuple[dict[str, dict[str, Any]], dict[str, int]]:
    lookup = TAGS.fetch_site_lookup(get_sidecar(sidecar_path, sidecar_index), chrom, pos)
    joined = TAGS.join_read_rows(reads_path, lookup)
    reads: dict[str, dict[str, Any]] = {}
    source_hp_replaced = 0
    ps_present = 0
    projection_multimatch = 0
    for source_row, latest, full_identity_count in joined:
        if full_identity_count != 1:
            raise RuntimeError(
                "Latest LongPhase-S projection is not a unique full alignment identity: "
                f"{source_row['read_name']} {source_row['chr']}:{source_row['start']}-{source_row['end']} "
                f"count={full_identity_count}"
            )
        read_id = source_row["read_id"]
        if read_id in reads:
            raise RuntimeError(f"Duplicate reads.tsv read_id: {reads_path} {read_id}")
        row: dict[str, Any] = dict(source_row)
        row["source_hp"] = TAGS.normalize_hp(source_row.get("hp"))
        row["hp"] = latest.hp
        row["ps"] = latest.ps
        row["latest_tag_full_identity_count"] = full_identity_count
        source_hp_replaced += row["source_hp"] != latest.hp
        ps_present += latest.ps is not None
        projection_multimatch += full_identity_count > 1
        reads[read_id] = row
    return reads, {
        "latest_tag_rows_fetched": lookup.rows_fetched,
        "latest_tag_rows_eligible": lookup.rows_eligible,
        "latest_tag_reads_joined": len(joined),
        "latest_tag_ps_present": ps_present,
        "latest_tag_projection_multimatch_reads": projection_multimatch,
        "source_hp_replaced_reads": source_hp_replaced,
    }


def readset_digest(rows: list[dict[str, Any]]) -> str:
    identities = sorted(
        f"{row['read_name']}|{row['chr']}|{row['start']}|{row['end']}|{row['strand']}" for row in rows
    )
    return hashlib.sha256("\n".join(identities).encode()).hexdigest()


def flatten_association(prefix: str, result: dict[str, Any], output: dict[str, Any]) -> None:
    for key, value in result.items():
        output[f"{prefix}_{key}"] = value


def configure_phylo_parallel(workers: int, min_reads: int) -> None:
    """Configure optional seed-level parallelism inside an outer site worker."""
    if workers < 1 or min_reads < 0:
        raise ValueError("phylo parallel workers must be positive and min reads nonnegative")
    if workers > 1 and min_reads < 2 * F.MIN_SIZE:
        raise ValueError("phylo parallel min reads must be at least the clustering minimum")
    global _PHYLO_PARALLEL_WORKERS, _PHYLO_PARALLEL_MIN_READS
    _PHYLO_PARALLEL_WORKERS = workers
    _PHYLO_PARALLEL_MIN_READS = min_reads


def _parallel_phylo_label(spec: tuple[str, int, int, float, int, str, float | None, float]) -> dict[str, Any]:
    kind, index, seed, null_pct, rnull, null_mode, empirical_alpha, min_valid_null_fraction = spec
    if _PARALLEL_PHYLO_DISTANCE is None or _PARALLEL_PHYLO_METHYLATION is None:
        raise RuntimeError("parallel phylo arrays are not initialized")
    trace: list[dict[str, Any]] = []
    labels = F.phylo_label(
        _PARALLEL_PHYLO_DISTANCE,
        _PARALLEL_PHYLO_METHYLATION,
        np.random.default_rng(seed),
        null_pct=null_pct,
        rnull=rnull,
        null_mode=null_mode,
        empirical_alpha=empirical_alpha,
        min_valid_null_fraction=min_valid_null_fraction,
        trace=trace,
    )
    return {"kind": kind, "index": index, "labels": labels, "trace": trace}


def analyze_phylo_parallel(
    distance: np.ndarray,
    methylation: np.ndarray,
    *,
    workers: int,
    base_seed: int = 20260622,
    seeds: int = 10,
    rnull: int = F.RNULL,
    null_mode: str = "column",
    empirical_alpha: float | None = None,
    min_valid_null_fraction: float = 0.8,
) -> dict[str, Any]:
    """Run independent coarse/fine phylo labels in forked workers without changing RNG streams."""
    if workers < 2:
        return F.analyze_phylo(
            distance,
            methylation,
            base_seed=base_seed,
            seeds=seeds,
            rnull=rnull,
            null_mode=null_mode,
            empirical_alpha=empirical_alpha,
            min_valid_null_fraction=min_valid_null_fraction,
        )
    if seeds < 1:
        raise ValueError("seeds must be positive")
    if "fork" not in multiprocessing.get_all_start_methods():
        raise RuntimeError("seed-level phylo parallelism requires the Linux fork start method")

    specs = [
        (
            "coarse",
            seed_index,
            base_seed + seed_index * 101,
            95.0,
            rnull,
            null_mode,
            empirical_alpha,
            min_valid_null_fraction,
        )
        for seed_index in range(seeds)
    ]
    specs.append(
        (
            "fine",
            0,
            base_seed,
            90.0,
            rnull,
            null_mode,
            empirical_alpha,
            min_valid_null_fraction,
        )
    )

    global _PARALLEL_PHYLO_DISTANCE, _PARALLEL_PHYLO_METHYLATION
    if _PARALLEL_PHYLO_DISTANCE is not None or _PARALLEL_PHYLO_METHYLATION is not None:
        raise RuntimeError("nested parallel phylo invocation is not supported")
    _PARALLEL_PHYLO_DISTANCE = distance
    _PARALLEL_PHYLO_METHYLATION = methylation
    try:
        context = multiprocessing.get_context("fork")
        with context.Pool(processes=min(workers, len(specs))) as pool:
            completed = pool.map(_parallel_phylo_label, specs, chunksize=1)
    finally:
        _PARALLEL_PHYLO_DISTANCE = None
        _PARALLEL_PHYLO_METHYLATION = None

    coarse = sorted(
        (run for run in completed if run["kind"] == "coarse"), key=lambda run: int(run["index"])
    )
    fine_runs = [run for run in completed if run["kind"] == "fine"]
    if len(coarse) != seeds or len(fine_runs) != 1:
        raise RuntimeError("parallel phylo run cardinality mismatch")
    coarse_runs = [run["labels"] for run in coarse]
    coarse_counts = [F.number_of_groups(labels) for labels in coarse_runs]
    frequencies = Counter(coarse_counts)
    modal_groups, modal_count = frequencies.most_common(1)[0]
    representative_index = next(index for index, count in enumerate(coarse_counts) if count == modal_groups)
    representative = coarse_runs[representative_index]
    fine_labels = fine_runs[0]["labels"]
    modal_fraction = modal_count / seeds
    modal_indices = [index for index, count in enumerate(coarse_counts) if count == modal_groups]
    assignment_aris = [
        float(F.adjusted_rand_score(coarse_runs[first], coarse_runs[second]))
        for offset, first in enumerate(modal_indices)
        for second in modal_indices[offset + 1 :]
    ]
    if assignment_aris:
        assignment_ari_median = float(np.median(assignment_aris))
        assignment_ari_min = min(assignment_aris)
    else:
        assignment_ari_median = 1.0
        assignment_ari_min = 1.0
    return {
        "coarse_ng": modal_groups,
        "modal_fraction": modal_fraction,
        "fine_ng": F.number_of_groups(fine_labels),
        "n_other": representative.count("other"),
        "n_outlier": representative.count("outlier"),
        "unstable": modal_fraction < F.MODAL_CONFIDENCE,
        "ng_min": min(coarse_counts),
        "ng_max": max(coarse_counts),
        "seed_group_counts": coarse_counts,
        "modal_assignment_pair_count": len(assignment_aris),
        "modal_assignment_ari_median": assignment_ari_median,
        "modal_assignment_ari_min": assignment_ari_min,
        "modal_assignment_all_pairs_ari_ge_0_8": assignment_ari_min >= 0.8,
        "hidden_heterogeneity": representative.count("other") / len(representative)
        > F.HIDDEN_HETEROGENEITY_FRACTION,
        "coarse_labels": representative,
        "fine_labels": fine_labels,
        "coarse_split_trace": coarse[representative_index]["trace"],
        "fine_split_trace": fine_runs[0]["trace"],
    }


def analyze_phylo_dispatch(distance: np.ndarray, methylation: np.ndarray) -> dict[str, Any]:
    if _PHYLO_PARALLEL_WORKERS > 1 and distance.shape[0] >= _PHYLO_PARALLEL_MIN_READS:
        return analyze_phylo_parallel(
            distance,
            methylation,
            workers=_PHYLO_PARALLEL_WORKERS,
            seeds=10,
            rnull=F.RNULL,
            null_mode="column",
            empirical_alpha=None,
        )
    return F.analyze_phylo(
        distance,
        methylation,
        seeds=10,
        rnull=F.RNULL,
        null_mode="column",
        empirical_alpha=None,
    )


def stable_cluster_detail(
    screen: dict[str, Any],
    reads: dict[str, dict[str, Any]],
    core_ids: list[str],
    core_labels: list[str],
    all_after_peel_ids: list[str],
    all_after_peel_labels: list[str],
    associations: dict[str, Any],
    phylo: dict[str, Any],
) -> dict[str, Any]:
    if not stable_null_multigroup_gate(phylo):
        raise RuntimeError("Stable assignment detail requested for a site that does not pass the M1 gate")
    region_dir = Path(screen["region_dir"])
    primary_artifacts = {
        "reads": artifact(region_dir / "reads" / "reads.tsv"),
        "distance_matrix": artifact(region_dir / "distance" / "BERNOULLI" / "matrix.csv"),
        "methylation_matrix": artifact(region_dir / "methylation" / "methylation.csv"),
    }
    core_reads = []
    for read_id, label in zip(core_ids, core_labels):
        row = reads[read_id]
        core_reads.append(
            {
                "read_id": read_id,
                "read_name": row["read_name"],
                "chrom": row["chr"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "mapq": int(row["mapq"]),
                "strand": row["strand"],
                "latest_hp": row["hp"],
                "latest_ps": row["ps"],
                "latest_tag_full_identity_count": int(row["latest_tag_full_identity_count"]),
                "label": label,
            }
        )
    return {
        "schema_name": "intersubmod.all_ssnv_stable_multigroup_read_assignment",
        "schema_version": "1.0.0",
        "screen_contract": SCREEN_CONTRACT,
        "m1_stability_gate_contract": M1_STABILITY_GATE_CONTRACT,
        "artifact_identity_contract": "sha256_size_path_v1",
        "primary_artifacts": primary_artifacts,
        "dataset": screen["sample"],
        "sample": screen["sample"],
        "chrom": screen["chrom"],
        "pos": int(screen["pos"]),
        "region_dir": screen["region_dir"],
        "core_reads": core_reads,
        "read_ids": core_ids,
        "read_names": [reads[read_id]["read_name"] for read_id in core_ids],
        "labels": core_labels,
        "latest_hp": [reads[read_id]["hp"] for read_id in core_ids],
        "latest_ps": [reads[read_id]["ps"] for read_id in core_ids],
        "associations": associations,
        "coarse_split_trace": phylo["coarse_split_trace"],
        "fine_split_trace": phylo["fine_split_trace"],
        "all_after_peel_read_ids": all_after_peel_ids,
        "coarse_labels_all_after_peel": all_after_peel_labels,
        "strict_confirm_status": "NOT_RUN",
        "strict_confirm_candidate": True,
        "strict_confirm_candidate_is_formal_r1_claim": False,
        "strict_confirm_reason": STRICT_NOT_RUN_REASON,
        "strict_methyl_partition_robustness_status": "NOT_EVALUATED_AT_M1_SCREEN",
        "strict_methyl_partition_robustness_not_evaluable_reason": STRICT_NOT_RUN_REASON,
    }


def reject_clustering_label_leakage(screen: dict[str, Any]) -> None:
    leaked = sorted(CLUSTERING_FORBIDDEN_FIELDS.intersection(screen))
    leaked.extend(
        sorted(
            field
            for field in screen
            if field.startswith(("truth_", "cooccurrence_", "partner_")) and field not in leaked
        )
    )
    if leaked:
        raise RuntimeError(f"Truth/cooccurrence fields are forbidden in clustering input: {leaked}")


def initial_screen_row(screen: dict[str, Any], reads: dict[str, dict[str, Any]], join_stats: dict[str, int]) -> dict[str, Any]:
    alt_rows = [
        row for row in reads.values() if F.is_tumor(row["is_tumor"]) and row["alt_support"] == "ALT"
    ]
    phase_tags = {"1-1", "2-1", "1-2", "2-2"}
    second_tags = {"1-2", "2-2"}
    ambiguous_tags = {"", ".", "0", "3", "4"}
    return {
        "dataset": screen["sample"],
        "sample": screen["sample"],
        "chrom": screen["chrom"],
        "pos": int(screen["pos"]),
        "region_dir": screen["region_dir"],
        "screen_contract": SCREEN_CONTRACT,
        "m1_stability_gate_contract": M1_STABILITY_GATE_CONTRACT,
        "screen_global_fdr_calibrated": False,
        "screen_rnull": F.RNULL,
        "screen_min_group_size": F.MIN_SIZE,
        "latest_tag_join_status": "PASS",
        **join_stats,
        "n_reads_total": len(reads),
        "n_alt_raw": len(alt_rows),
        "n_alt_matrix": 0,
        "n_alt_after_peel": 0,
        "n_alt_peeled": 0,
        "alt_readset_sha256": readset_digest(alt_rows),
        "alt_hp_counts": json_text(Counter(row["hp"] for row in alt_rows)),
        "alt_hp_family_counts": json_text(Counter(F.hp_family(row["hp"]) for row in alt_rows)),
        "alt_ps_counts": json_text(Counter(str(row["ps"]) if row["ps"] is not None else "." for row in alt_rows)),
        "alt_strand_counts": json_text(Counter(row["strand"] for row in alt_rows)),
        "phase_anchored_fraction_raw": (
            sum(row["hp"] in phase_tags for row in alt_rows) / len(alt_rows) if alt_rows else None
        ),
        "explicit_second_haplotype_fraction_raw": (
            sum(row["hp"] in second_tags for row in alt_rows) / len(alt_rows) if alt_rows else None
        ),
        "ambiguous_tag_fraction_raw": (
            sum(row["hp"] in ambiguous_tags for row in alt_rows) / len(alt_rows) if alt_rows else None
        ),
        "m1_evaluable": False,
        "m1_not_evaluable_reason": "INSUFFICIENT_FOCAL_ALT_READS",
        "analysis_status": "insufficient_alt_reads",
        "forced_k": 0,
        "forced_silhouette": None,
        "coarse_ng": 0,
        "modal_fraction": None,
        "fine_ng": 0,
        "n_other": 0,
        "n_outlier": 0,
        "unstable": None,
        "hidden_heterogeneity": None,
        "modal_assignment_pair_count": None,
        "modal_assignment_ari_median": None,
        "modal_assignment_ari_min": None,
        "modal_assignment_all_pairs_ari_ge_0_8": None,
        "stable_null_multigroup": False,
        "cluster_sizes": "{}",
        "seed_group_counts": "[]",
        "hp_axis_confound": False,
        "technical_axis_confound": False,
        "residual_unexplained_multigroup": False,
        "phase_anchored_fraction_core": None,
        "phase_anchored_robust_epigenetic_candidate": False,
        "evidence_tier": "E0_NOT_EVALUABLE",
        "strict_methyl_partition_robustness_status": "NOT_EVALUATED_AT_M1_SCREEN",
        "strict_methyl_partition_robustness_not_evaluable_reason": STRICT_NOT_RUN_REASON,
        "strict_confirm_status": "NOT_RUN",
        "strict_confirm_candidate": False,
        "strict_confirm_candidate_is_formal_r1_claim": False,
        "strict_confirm_reason": STRICT_NOT_RUN_REASON,
    }


def screen_site(screen: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """Run the truth-blind methylation screen using only read/matrix/tag inputs."""
    reject_clustering_label_leakage(screen)
    reads, join_stats = load_latest_joined_reads(
        Path(screen["reads_path"]),
        screen["sidecar_path"],
        screen["sidecar_index"],
        screen["chrom"],
        int(screen["pos"]),
    )
    output = initial_screen_row(screen, reads, join_stats)
    distance_ids, distance = F.load_matrix(Path(screen["distance_path"]))
    methylation_ids, methylation = F.load_matrix(Path(screen["methylation_path"]))
    if len(distance_ids) != len(set(distance_ids)) or len(methylation_ids) != len(set(methylation_ids)):
        raise RuntimeError("Duplicate read IDs in InterSubMod matrix")
    if distance.shape != (len(distance_ids), len(distance_ids)):
        raise RuntimeError(f"Invalid C++ BERNOULLI matrix shape: {distance.shape}")
    if methylation.shape[0] != len(methylation_ids):
        raise RuntimeError(f"Invalid methylation matrix row count: {methylation.shape}")
    distance_index = {read_id: index for index, read_id in enumerate(distance_ids)}
    methylation_index = {read_id: index for index, read_id in enumerate(methylation_ids)}
    alt_ids = [
        read_id
        for read_id in distance_ids
        if read_id in reads
        and read_id in methylation_index
        and F.is_tumor(reads[read_id]["is_tumor"])
        and reads[read_id]["alt_support"] == "ALT"
    ]
    output["n_alt_matrix"] = len(alt_ids)
    if len(alt_ids) < 2 * F.MIN_SIZE:
        output["m1_not_evaluable_reason"] = "INSUFFICIENT_MATRIX_JOINED_FOCAL_ALT_READS"
        return output, None

    distance_rows = [distance_index[read_id] for read_id in alt_ids]
    sub_distance = distance[np.ix_(distance_rows, distance_rows)]
    kept = F.peel_complete(sub_distance)
    kept_ids = [alt_ids[index] for index in kept]
    output["n_alt_after_peel"] = len(kept_ids)
    output["n_alt_peeled"] = len(alt_ids) - len(kept_ids)
    if len(kept_ids) < 2 * F.MIN_SIZE:
        output["m1_not_evaluable_reason"] = "INCOMPLETE_DISTANCE_BELOW_MINIMUM"
        output["analysis_status"] = "incomplete_distance_below_minimum"
        return output, None

    sub_distance = sub_distance[np.ix_(kept, kept)]
    sub_methylation = methylation[[methylation_index[read_id] for read_id in kept_ids]]
    forced = F.forced_silhouette_split(sub_distance.copy())
    phylo = analyze_phylo_dispatch(sub_distance, sub_methylation)
    output["forced_k"] = forced["k"]
    output["forced_silhouette"] = forced["silhouette"]
    for key in (
        "coarse_ng",
        "modal_fraction",
        "fine_ng",
        "n_other",
        "n_outlier",
        "unstable",
        "hidden_heterogeneity",
        "modal_assignment_pair_count",
        "modal_assignment_ari_median",
        "modal_assignment_ari_min",
        "modal_assignment_all_pairs_ari_ge_0_8",
    ):
        output[key] = phylo[key]
    output["seed_group_counts"] = json_text(phylo["seed_group_counts"])
    output["m1_evaluable"] = True
    output["m1_not_evaluable_reason"] = ""
    output["analysis_status"] = "evaluable"
    output["evidence_tier"] = "E1_EVALUABLE_OR_SILHOUETTE_ONLY"
    apply_m1_screen_flags(output, phylo)

    labels = phylo["coarse_labels"]
    core_indices = [index for index, label in enumerate(labels) if label not in {"other", "outlier"}]
    core_labels = [labels[index] for index in core_indices]
    core_ids = [kept_ids[index] for index in core_indices]
    output["cluster_sizes"] = json_text(Counter(core_labels))
    if not output["stable_null_multigroup"]:
        return output, None

    output["evidence_tier"] = "E2_STABLE_NULL_MULTIGROUP"
    core_rows = [reads[read_id] for read_id in core_ids]
    cpg_called = [int(np.isfinite(sub_methylation[index]).sum()) for index in core_indices]
    seed = F.stable_seed(screen["sample"], screen["chrom"], int(screen["pos"]))
    associations = {
        "hp_exact": F.categorical_permutation_association(
            [row["hp"] for row in core_rows], core_labels, seed + 1
        ),
        "hp_family": F.categorical_permutation_association(
            [F.hp_family(row["hp"]) for row in core_rows], core_labels, seed + 2
        ),
        "strand": F.categorical_permutation_association(
            [row["strand"] for row in core_rows], core_labels, seed + 3
        ),
        "start": F.continuous_permutation_association(
            [float(row["start"]) for row in core_rows], core_labels, seed + 4
        ),
        "end": F.continuous_permutation_association(
            [float(row["end"]) for row in core_rows], core_labels, seed + 5
        ),
        "length": F.continuous_permutation_association(
            [float(row["end"]) - float(row["start"]) for row in core_rows], core_labels, seed + 6
        ),
        "mapq": F.continuous_permutation_association(
            [float(row["mapq"]) for row in core_rows], core_labels, seed + 7
        ),
        "cpg_called": F.continuous_permutation_association(cpg_called, core_labels, seed + 8),
    }
    for name, association in associations.items():
        flatten_association(name, association, output)
    output["hp_axis_confound"] = associations["hp_exact"]["aligned"] or associations["hp_family"]["aligned"]
    output["technical_axis_confound"] = any(
        associations[name]["aligned"]
        for name in ("strand", "start", "end", "length", "mapq", "cpg_called")
    )
    output["residual_unexplained_multigroup"] = not output["hp_axis_confound"] and not output[
        "technical_axis_confound"
    ]
    if output["residual_unexplained_multigroup"]:
        output["evidence_tier"] = "E3_UNEXPLAINED_AFTER_MEASURED_AXES"
    phase_tags = {"1-1", "2-1", "1-2", "2-2"}
    core_phase_fraction = sum(row["hp"] in phase_tags for row in core_rows) / len(core_rows)
    core_sizes = Counter(core_labels)
    output["phase_anchored_fraction_core"] = core_phase_fraction
    output["phase_anchored_robust_epigenetic_candidate"] = bool(
        output["residual_unexplained_multigroup"]
        and len(core_rows) >= 10
        and phylo["modal_fraction"] == 1.0
        and core_sizes
        and min(core_sizes.values()) >= 5
        and core_phase_fraction >= 0.80
    )
    if output["phase_anchored_robust_epigenetic_candidate"]:
        output["evidence_tier"] = "E4_PHASE_ANCHORED_ROBUST_EPIGENETIC_CANDIDATE"
    detail = stable_cluster_detail(
        screen,
        reads,
        core_ids,
        core_labels,
        kept_ids,
        labels,
        associations,
        phylo,
    )
    return output, detail


def attach_posthoc(
    screen_row: dict[str, Any],
    detail: dict[str, Any] | None,
    posthoc: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    identity_fields = {"chrom", "pos"}
    for field in identity_fields:
        if field in screen_row and field in posthoc and screen_row[field] != posthoc[field]:
            raise RuntimeError(f"Screen/posthoc identity mismatch for {field}")
    row = {**posthoc, **screen_row}
    if detail is not None:
        detail = dict(detail)
        detail["posthoc"] = {
            key: posthoc.get(key)
            for key in (
                "biological_id",
                "truth_label",
                "ref",
                "alt",
                "ssnv_branch",
                "component_id",
                "component_size",
                "selected_group_id",
                "selected_group_size",
            )
        }
    return row, detail


def process_site_task(task: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any] | None]:
    screen = task["screen"]
    try:
        screen_row, detail = screen_site(screen)
        return attach_posthoc(screen_row, detail, task["posthoc"])
    except Exception as error:
        raise RuntimeError(
            f"{screen['sample']}:{screen['chrom']}:{screen['pos']} screening failed: {error}"
        ) from error


def process_chunk(chunk: list[dict[str, Any]]) -> list[tuple[dict[str, Any], dict[str, Any] | None]]:
    return [process_site_task(task) for task in chunk]


def chunked(values: Iterable[dict[str, Any]], size: int) -> Iterator[list[dict[str, Any]]]:
    if size < 1:
        raise ValueError("Chunk size must be positive")
    chunk: list[dict[str, Any]] = []
    for value in values:
        chunk.append(value)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def bounded_chunk_results(
    executor: ProcessPoolExecutor,
    chunks: Iterable[list[dict[str, Any]]],
    max_pending: int,
) -> Iterator[list[tuple[dict[str, Any], dict[str, Any] | None]]]:
    """Submit at most max_pending chunk futures and emit results in input order."""
    if max_pending < 1:
        raise ValueError("max_pending must be positive")
    iterator = iter(chunks)
    pending: deque[Future[list[tuple[dict[str, Any], dict[str, Any] | None]]]] = deque()

    def fill() -> None:
        while len(pending) < max_pending:
            try:
                pending.append(executor.submit(process_chunk, next(iterator)))
            except StopIteration:
                break

    fill()
    while pending:
        yield pending.popleft().result()
        fill()


def wilson(successes: int, total: int, z: float = 1.959963984540054) -> list[float | None]:
    if total <= 0:
        return [None, None]
    proportion = successes / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    radius = z * math.sqrt(proportion * (1 - proportion) / total + z * z / (4 * total * total)) / denominator
    return [max(0.0, center - radius), min(1.0, center + radius)]


class ScreenSummary:
    def __init__(self) -> None:
        self.n_sites = 0
        self.status: Counter[str] = Counter()
        self.evidence: Counter[str] = Counter()
        self.branches: Counter[str] = Counter()
        self.evaluable = 0
        self.stable = 0
        self.residual = 0
        self.high_threshold = 0
        self.strict_candidates = 0
        self.strict_not_run = 0
        self.latest_tag_join_status: Counter[str] = Counter()
        self.latest_tag_rows_fetched = 0
        self.latest_tag_rows_eligible = 0
        self.latest_tag_reads_joined = 0
        self.latest_tag_ps_present = 0
        self.latest_tag_projection_multimatch_reads = 0
        self.source_hp_replaced_reads = 0
        self.reads_tsv_site_rows = 0
        self.sites_with_zero_reads = 0

    def add(self, row: dict[str, Any]) -> None:
        join_status = str(row["latest_tag_join_status"])
        joined = int(row["latest_tag_reads_joined"])
        reads_tsv_rows = int(row["n_reads_total"])
        tag_counts = {
            "latest_tag_rows_fetched": int(row["latest_tag_rows_fetched"]),
            "latest_tag_rows_eligible": int(row["latest_tag_rows_eligible"]),
            "latest_tag_reads_joined": joined,
            "latest_tag_ps_present": int(row["latest_tag_ps_present"]),
            "latest_tag_projection_multimatch_reads": int(
                row["latest_tag_projection_multimatch_reads"]
            ),
            "source_hp_replaced_reads": int(row["source_hp_replaced_reads"]),
        }
        if join_status != "PASS":
            raise RuntimeError(f"Latest HP/PS join did not pass: {join_status}")
        if any(value < 0 for value in tag_counts.values()) or reads_tsv_rows < 0:
            raise RuntimeError("Latest HP/PS join counts must be nonnegative")
        if joined != reads_tsv_rows:
            raise RuntimeError(
                "Latest HP/PS join did not consume every reads.tsv row: "
                f"joined={joined} reads_tsv_rows={reads_tsv_rows}"
            )
        if not (
            tag_counts["latest_tag_rows_fetched"]
            >= tag_counts["latest_tag_rows_eligible"]
            >= joined
        ):
            raise RuntimeError("Latest HP/PS sidecar fetched/eligible/joined counts are inconsistent")
        for field in (
            "latest_tag_ps_present",
            "latest_tag_projection_multimatch_reads",
            "source_hp_replaced_reads",
        ):
            if tag_counts[field] > joined:
                raise RuntimeError(f"Latest HP/PS count exceeds joined reads: {field}")
        self.n_sites += 1
        self.status[row["analysis_status"]] += 1
        self.evidence[row["evidence_tier"]] += 1
        self.branches[row.get("ssnv_branch") or "NA"] += 1
        self.evaluable += row["analysis_status"] == "evaluable"
        self.stable += bool(row["stable_null_multigroup"])
        self.residual += bool(row["residual_unexplained_multigroup"])
        self.high_threshold += bool(row["phase_anchored_robust_epigenetic_candidate"])
        self.strict_candidates += bool(row["strict_confirm_candidate"])
        self.strict_not_run += row["strict_confirm_status"] == "NOT_RUN"
        self.latest_tag_join_status[join_status] += 1
        self.latest_tag_rows_fetched += tag_counts["latest_tag_rows_fetched"]
        self.latest_tag_rows_eligible += tag_counts["latest_tag_rows_eligible"]
        self.latest_tag_reads_joined += joined
        self.latest_tag_ps_present += tag_counts["latest_tag_ps_present"]
        self.latest_tag_projection_multimatch_reads += tag_counts[
            "latest_tag_projection_multimatch_reads"
        ]
        self.source_hp_replaced_reads += tag_counts["source_hp_replaced_reads"]
        self.reads_tsv_site_rows += reads_tsv_rows
        self.sites_with_zero_reads += reads_tsv_rows == 0

    def latest_tag_payload(self) -> dict[str, Any]:
        all_sites_pass = self.latest_tag_join_status == Counter({"PASS": self.n_sites})
        every_row_joined = self.latest_tag_reads_joined == self.reads_tsv_site_rows
        return {
            "authoritative_tag_source": "same_run_LongPhase_S_external_HP_PS_sidecar",
            "embedded_reads_tsv_hp_used_for_analysis": False,
            "join_occurs_before_focal_ALT_selection": True,
            "counting_unit": "site_read_observation_not_globally_unique_read",
            "n_sites": self.n_sites,
            "join_status_counts": dict(self.latest_tag_join_status),
            "all_sites_pass": all_sites_pass,
            "n_reads_tsv_site_rows": self.reads_tsv_site_rows,
            "n_exact_hp_ps_site_read_joins": self.latest_tag_reads_joined,
            "every_reads_tsv_row_joined": every_row_joined,
            "n_ps_present_site_read_joins": self.latest_tag_ps_present,
            "n_source_hp_replaced_site_read_joins": self.source_hp_replaced_reads,
            "n_sidecar_rows_fetched_site_observations": self.latest_tag_rows_fetched,
            "n_sidecar_rows_eligible_site_observations": self.latest_tag_rows_eligible,
            "n_projection_multimatch_site_reads": self.latest_tag_projection_multimatch_reads,
            "all_projection_identities_unique": self.latest_tag_projection_multimatch_reads == 0,
            "n_sites_with_zero_reads_tsv_rows": self.sites_with_zero_reads,
            "pass": all_sites_pass
            and every_row_joined
            and self.latest_tag_projection_multimatch_reads == 0,
        }

    def payload(self) -> dict[str, Any]:
        return {
            "n_sites": self.n_sites,
            "status_counts": dict(self.status),
            "evidence_tier_counts": dict(self.evidence),
            "ssnv_branch_counts": dict(self.branches),
            "n_evaluable": self.evaluable,
            "evaluable_fraction_all": self.evaluable / self.n_sites if self.n_sites else None,
            "n_stable_null_multigroup": self.stable,
            "stable_fraction_evaluable": self.stable / self.evaluable if self.evaluable else None,
            "stable_fraction_evaluable_wilson95": wilson(self.stable, self.evaluable),
            "n_residual_unexplained_multigroup": self.residual,
            "residual_fraction_evaluable": self.residual / self.evaluable if self.evaluable else None,
            "n_phase_anchored_robust_epigenetic_candidate": self.high_threshold,
            "high_threshold_fraction_evaluable": (
                self.high_threshold / self.evaluable if self.evaluable else None
            ),
            "n_strict_confirm_candidates": self.strict_candidates,
            "n_legacy_strict_confirm_candidate_alias": self.strict_candidates,
            "strict_confirm_status_counts": {"NOT_RUN": self.strict_not_run},
            "latest_hp_ps_terminal_join_audit": self.latest_tag_payload(),
        }


def validate_manifest(manifest: dict[str, Any]) -> None:
    if manifest.get("schema_name") != "intersubmod.all_ssnv_focal_alt_input_manifest":
        raise RuntimeError("Unexpected all-sSNV input manifest schema")
    if manifest.get("pass") is not True:
        raise RuntimeError("All-sSNV input manifest is not passing")
    samples = manifest.get("samples", [])
    if [entry.get("sample") for entry in samples] != DATASETS:
        raise RuntimeError("All-sSNV manifest dataset order/content drift")
    if int(manifest["totals"]["all_ssnv"]) != EXPECTED_ALL_SSNV:
        raise RuntimeError("All-sSNV manifest total drift")


def create_output_dir(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing result directory: {path}")
    path.mkdir(parents=True, exist_ok=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=TOPIC_ROOT / "results" / "all_ssnv_input_manifest.json",
    )
    parser.add_argument("--intersubmod-root", type=Path, default=None)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=TOPIC_ROOT / "results" / "all_ssnv_focal_alt_multigroup_v1",
    )
    parser.add_argument("--sample", action="append", choices=DATASETS)
    parser.add_argument("--workers", type=int, default=max(1, min(32, (os.cpu_count() or 4) - 4)))
    parser.add_argument("--chunk-size", type=int, default=128)
    parser.add_argument("--max-pending-chunks", type=int, default=0)
    parser.add_argument("--progress-every", type=int, default=1000)
    parser.add_argument(
        "--phylo-parallel-workers",
        type=int,
        default=1,
        help="Fork workers for independent coarse/fine runs at very high-depth sites (1 disables)",
    )
    parser.add_argument(
        "--phylo-parallel-min-reads",
        type=int,
        default=0,
        help="Minimum peeled focal-ALT reads for seed-level parallelism",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_at = now_utc()
    if args.workers < 1 or args.chunk_size < 1 or args.max_pending_chunks < 0 or args.progress_every < 1:
        raise ValueError("workers, chunk-size, and progress-every must be positive; max-pending-chunks must be nonnegative")
    configure_phylo_parallel(args.phylo_parallel_workers, args.phylo_parallel_min_reads)
    logical_cpus = os.cpu_count() or 1
    maximum_processes = (
        args.workers * args.phylo_parallel_workers
        if args.phylo_parallel_workers > 1
        else args.workers
    )
    if maximum_processes > logical_cpus:
        raise ValueError(
            "outer workers plus nested phylo workers exceed logical CPUs: "
            f"{maximum_processes} > {logical_cpus}"
        )
    validate_prior_screen_contract()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    validate_manifest(manifest)
    selected = set(args.sample or DATASETS)
    entries = [entry for entry in manifest["samples"] if entry["sample"] in selected]
    expected_sites = sum(int(entry["counts"]["all_ssnv"]) for entry in entries)
    full_scope = selected == set(DATASETS)
    if full_scope and expected_sites != EXPECTED_ALL_SSNV:
        raise RuntimeError(f"Full-scope selected total drift: {expected_sites}")
    intersubmod_root = args.intersubmod_root or (Path(manifest["workspace"]) / "intersubmod_all_ssnv_v1")
    completion_receipts = {
        entry["sample"]: validate_completed_sample_run(intersubmod_root / entry["sample"], entry)
        for entry in entries
    }
    create_output_dir(args.output_dir)

    site_path = args.output_dir / "all_ssnv_site_results.tsv.gz"
    assignment_path = args.output_dir / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
    summary_path = args.output_dir / "all_ssnv_summary.json"
    run_manifest_path = args.output_dir / "run_manifest.json"
    pooled = ScreenSummary()
    by_sample = {entry["sample"]: ScreenSummary() for entry in entries}
    by_truth = {label: ScreenSummary() for label in sorted(TRUTH_LABELS)}
    by_biological_id = {
        biological_id: ScreenSummary()
        for biological_id in sorted({entry["biological_id"] for entry in entries})
    }
    by_ledger_branch: dict[str, ScreenSummary] = {}
    processed = 0
    assignments = 0

    def task_stream() -> Iterator[dict[str, Any]]:
        for entry in entries:
            yield from iter_sample_tasks(entry, intersubmod_root)

    max_pending = args.max_pending_chunks or max(1, args.workers * 2)
    with gzip.open(site_path, "wt", encoding="utf-8", newline="") as site_handle, gzip.open(
        assignment_path, "wt", encoding="utf-8", newline=""
    ) as assignment_handle:
        writer = csv.DictWriter(site_handle, SITE_FIELDS, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        chunks = chunked(task_stream(), args.chunk_size)
        with ProcessPoolExecutor(
            max_workers=max(1, args.workers),
            initializer=configure_phylo_parallel,
            initargs=(args.phylo_parallel_workers, args.phylo_parallel_min_reads),
        ) as executor:
            for results in bounded_chunk_results(executor, chunks, max_pending):
                for row, detail in results:
                    writer.writerow(row)
                    if detail is not None:
                        assignment_handle.write(json.dumps(detail, ensure_ascii=False, separators=(",", ":")) + "\n")
                        assignments += 1
                    pooled.add(row)
                    by_sample[row["sample"]].add(row)
                    by_truth[row["truth_label"]].add(row)
                    by_biological_id[row["biological_id"]].add(row)
                    branch = row.get("ssnv_branch") or "NA"
                    if branch not in by_ledger_branch:
                        by_ledger_branch[branch] = ScreenSummary()
                    by_ledger_branch[branch].add(row)
                    processed += 1
                    if processed % args.progress_every == 0 or processed == expected_sites:
                        print(
                            f"processed={processed}/{expected_sites} stable={pooled.stable} "
                            f"assignments={assignments}",
                            flush=True,
                        )

    if processed != expected_sites:
        raise RuntimeError(f"Processed site count mismatch: {processed} != {expected_sites}")
    if assignments != pooled.stable:
        raise RuntimeError(f"Stable assignment count mismatch: {assignments} != {pooled.stable}")
    latest_tag_audit = pooled.latest_tag_payload()
    if not latest_tag_audit["pass"]:
        raise RuntimeError("Pooled latest HP/PS terminal join audit failed")
    summary = {
        "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_screen",
        "schema_version": OUTPUT_SCHEMA_VERSION,
        "created_at_utc": now_utc(),
        "status": "EXECUTION_PASS",
        "pass_semantics": "execution_integrity_only_not_scientific_confirmation",
        "scope": {
            "full_469849": full_scope,
            "selected_datasets": [entry["sample"] for entry in entries],
            "selected_samples": [entry["sample"] for entry in entries],
            "expected_sites": expected_sites,
            "processed_sites": processed,
        },
        "population": "LongPhase-S recalibrated FILTER=PASS chr1-22 biallelic sSNVs",
        "clustering_contract": {
            "truth_blind": True,
            "cooccurrence_blind": True,
            "distance_source": "exact existing C++ BERNOULLI matrices",
            "screen": SCREEN_CONTRACT,
            "screen_contract_semantics": "legacy_algorithm_identity_for_downstream_compatibility",
            "m1_stability_gate_contract": M1_STABILITY_GATE_CONTRACT,
            "prior_screen_thresholds": PRIOR_SCREEN_THRESHOLDS,
            "global_fdr_calibrated": False,
            "strict_methyl_partition_robustness_status": "NOT_EVALUATED_AT_M1_SCREEN",
            "strict_confirm_status_legacy_alias": "NOT_RUN",
            "stable_null_multigroup_basis": M1_STABILITY_GATE_CONTRACT,
            "strict_confirm_candidate_legacy_alias_basis": "same_as_stable_null_multigroup",
            "strict_confirm_candidate_is_formal_r1_claim": False,
        },
        "pooled_site_weighted": pooled.payload(),
        "per_dataset": {sample: stats.payload() for sample, stats in by_sample.items()},
        "per_sample": {sample: stats.payload() for sample, stats in by_sample.items()},
        "posthoc_truth_strata": {label: by_truth[label].payload() for label in sorted(TRUTH_LABELS)},
        "posthoc_biological_id_strata": {
            biological_id: stats.payload() for biological_id, stats in by_biological_id.items()
        },
        "posthoc_ledger_branch_strata": {
            branch: stats.payload() for branch, stats in sorted(by_ledger_branch.items())
        },
        "n_stable_assignment_records": assignments,
        "latest_hp_ps_terminal_join_audit": latest_tag_audit,
        "interpretation_guardrail": (
            "This is an M1 screen, not a genome-wide FDR procedure. A focal-ALT stable methyl "
            "multigroup supports ALT-carrier read-level methyl heterogeneity only, not a genetic clone, "
            "cell fraction, or lineage topology. R1 strict methyl-partition robustness is not evaluated "
            "at this stage; not evaluated is not failure."
        ),
        "pass": processed == expected_sites
        and assignments == pooled.stable
        and latest_tag_audit["pass"],
    }
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    finished_at = now_utc()
    run_manifest = {
        "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_run_manifest",
        "schema_version": OUTPUT_SCHEMA_VERSION,
        "created_at_utc": finished_at,
        "status": "EXECUTION_PASS",
        "pass_semantics": "execution_integrity_only_not_scientific_confirmation",
        "started_at_utc": started_at,
        "finished_at_utc": finished_at,
        "command": sys.argv,
        "input_manifest": artifact(args.manifest),
        "intersubmod_root": str(intersubmod_root.resolve()),
        "completed_dataset_runs": completion_receipts,
        "completed_sample_runs": completion_receipts,
        "source_code": {
            "analyzer": artifact(Path(__file__)),
            "focal_alt_cluster_lib": artifact(Path(F.__file__)),
            "latest_tag_join": artifact(Path(TAGS.__file__)),
            "claim_contract_v2": artifact(TOPIC_ROOT / "claim-contract-v2.md"),
        },
        "execution": {
            "workers": args.workers,
            "chunk_size": args.chunk_size,
            "max_pending_chunks": max_pending,
            "per_site_future_submission": False,
            "phylo_parallel_workers": args.phylo_parallel_workers,
            "phylo_parallel_min_reads": args.phylo_parallel_min_reads,
            "maximum_processes_by_contract": maximum_processes,
            "phylo_parallel_semantics": (
                "independent coarse/fine runs only; per-run RNG seed and null sequence unchanged"
            ),
            "selected_datasets": [entry["sample"] for entry in entries],
            "selected_samples": [entry["sample"] for entry in entries],
        },
        "contracts": {
            "truth_and_cooccurrence_enter_clustering": False,
            "latest_hp_ps_join": "same-run sidecar projected join before ALT selection; missing/conflict hard fail",
            "screen_global_fdr_calibrated": False,
            "screen_contract_semantics": "legacy_algorithm_identity_for_downstream_compatibility",
            "m1_stability_gate_contract": M1_STABILITY_GATE_CONTRACT,
            "strict_methyl_partition_robustness_status": "NOT_EVALUATED_AT_M1_SCREEN",
            "strict_confirm_status_legacy_alias": "NOT_RUN",
            "stable_null_multigroup_basis": M1_STABILITY_GATE_CONTRACT,
            "strict_confirm_candidate_legacy_alias_basis": "same_as_stable_null_multigroup",
            "strict_confirm_candidate_is_formal_r1_claim": False,
            "prior_screen_thresholds": PRIOR_SCREEN_THRESHOLDS,
            "existing_results_overwritten": False,
            "latest_hp_ps_terminal_join": latest_tag_audit,
        },
        "outputs": {
            "site_results": artifact(site_path),
            "stable_assignments": artifact(assignment_path),
            "summary": artifact(summary_path),
        },
        "counts": {
            "expected_sites": expected_sites,
            "processed_sites": processed,
            "stable_assignment_records": assignments,
            "reads_tsv_site_rows": latest_tag_audit["n_reads_tsv_site_rows"],
            "exact_hp_ps_site_read_joins": latest_tag_audit[
                "n_exact_hp_ps_site_read_joins"
            ],
            "ps_present_site_read_joins": latest_tag_audit[
                "n_ps_present_site_read_joins"
            ],
            "source_hp_replaced_site_read_joins": latest_tag_audit[
                "n_source_hp_replaced_site_read_joins"
            ],
        },
        "pass": summary["pass"],
    }
    run_manifest_path.write_text(
        json.dumps(run_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "output_dir": str(args.output_dir.resolve()),
                "summary": str(summary_path.resolve()),
                "processed_sites": processed,
                "stable_assignments": assignments,
                "pass": summary["pass"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
