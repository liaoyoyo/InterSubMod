#!/usr/bin/env python3
"""Merge a complete canonical prefix with one complete replacement sample run."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping

import analyze_all_ssnv_focal_alt_multigroup as A


SCHEMA_NAME = "intersubmod.all_ssnv_focal_alt_multigroup_recovery_merge"
EQUIVALENCE_SCHEMA = "intersubmod.phylo_parallel_exact_equivalence"
SOURCE_LOCK_SCHEMA = "intersubmod.source_locked_screen_run"
BOOL_FIELDS = {
    "stable_null_multigroup",
    "residual_unexplained_multigroup",
    "phase_anchored_robust_epigenetic_candidate",
    "strict_confirm_candidate",
}


def parse_bool(value: Any, field: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "1"}:
        return True
    if normalized in {"false", "0"}:
        return False
    raise RuntimeError(f"Invalid boolean for {field}: {value!r}")


def summary_row(row: Mapping[str, Any]) -> dict[str, Any]:
    normalized = dict(row)
    for field in BOOL_FIELDS:
        normalized[field] = parse_bool(normalized[field], field)
    return normalized


def site_key(row: Mapping[str, Any]) -> tuple[str, str, int, str, str]:
    return (
        str(row["sample"]),
        str(row["chrom"]),
        int(row["pos"]),
        str(row["ref"]).upper(),
        str(row["alt"]).upper(),
    )


def assignment_key(row: Mapping[str, Any]) -> tuple[str, str, int, str, str]:
    posthoc = row.get("posthoc")
    if not isinstance(posthoc, Mapping):
        raise RuntimeError("Stable assignment is missing posthoc identity")
    return (
        str(row["sample"]),
        str(row["chrom"]),
        int(row["pos"]),
        str(posthoc["ref"]).upper(),
        str(posthoc["alt"]).upper(),
    )


def logical_digest_update(digest: Any, payload: str) -> None:
    digest.update(payload.encode("utf-8"))
    digest.update(b"\n")


def parse_utc(value: Any, label: str) -> datetime:
    if not isinstance(value, str) or not value.strip():
        raise RuntimeError(f"Missing timestamp for {label}")
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise RuntimeError(f"Naive timestamp for {label}")
    return parsed.astimezone(timezone.utc)


def validate_current_artifact(reference: Any, label: str) -> Path:
    if not isinstance(reference, Mapping):
        raise RuntimeError(f"Missing artifact reference for {label}")
    path = Path(str(reference.get("path", "")))
    if not path.is_file() or A.artifact(path) != dict(reference):
        raise RuntimeError(f"Artifact identity mismatch for {label}")
    return path


def iter_prefix_sites(path: Path, replacement_sample: str) -> Iterator[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != A.SITE_FIELDS:
            raise RuntimeError("Prefix site header differs from the canonical screen schema")
        for row in reader:
            if row["sample"] == replacement_sample:
                break
            yield row


def iter_prefix_assignments(path: Path, replacement_sample: str) -> Iterator[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("sample") == replacement_sample:
                break
            if row.get("sample") not in A.DATASETS:
                raise RuntimeError(f"Unexpected prefix assignment sample at line {line_number}")
            yield row


def iter_complete_sites(path: Path, expected_sample: str) -> Iterator[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != A.SITE_FIELDS:
            raise RuntimeError("Replacement site header differs from the canonical screen schema")
        for row in reader:
            if row["sample"] != expected_sample:
                raise RuntimeError("Replacement site file contains another sample")
            yield row


def iter_complete_assignments(path: Path, expected_sample: str) -> Iterator[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("sample") != expected_sample:
                raise RuntimeError(f"Replacement assignment sample drift at line {line_number}")
            yield row


def validate_prefix(
    prefix_dir: Path,
    manifest_path: Path,
    expected_samples: list[str],
    expected_sites: int,
) -> tuple[Path, Path, dict[str, Any], dict[str, Any], dict[str, Any]]:
    site_path = prefix_dir / "all_ssnv_site_results.tsv.gz"
    assignment_path = prefix_dir / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
    summary_path = prefix_dir / "all_ssnv_summary.json"
    receipt_path = prefix_dir / "run_manifest.json"
    source_lock_path = prefix_dir / "source_lock_receipt.json"
    for path in (site_path, assignment_path, summary_path, receipt_path, source_lock_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    source_lock = json.loads(source_lock_path.read_text(encoding="utf-8"))
    if summary.get("schema_name") != "intersubmod.all_ssnv_focal_alt_multigroup_screen":
        raise RuntimeError("Prefix summary schema mismatch")
    if receipt.get("schema_name") != "intersubmod.all_ssnv_focal_alt_multigroup_run_manifest":
        raise RuntimeError("Prefix run-manifest schema mismatch")
    if summary.get("pass") is not True or receipt.get("pass") is not True:
        raise RuntimeError("Prefix run is not passing")
    scope = summary.get("scope") or {}
    if scope.get("full_469849") is not False:
        raise RuntimeError("Prefix summary must be an explicit subset run")
    if scope.get("selected_samples") != expected_samples:
        raise RuntimeError("Prefix selected-sample order mismatch")
    if int(scope.get("processed_sites", -1)) != expected_sites:
        raise RuntimeError("Prefix summary site count mismatch")
    outputs = receipt.get("outputs") or {}
    for name, path in {
        "site_results": site_path,
        "stable_assignments": assignment_path,
        "summary": summary_path,
    }.items():
        if outputs.get(name) != A.artifact(path):
            raise RuntimeError(f"Prefix receipt artifact mismatch: {name}")
    if receipt.get("input_manifest") != A.artifact(manifest_path):
        raise RuntimeError("Prefix input-manifest identity mismatch")
    execution = receipt.get("execution") or {}
    if execution.get("selected_samples") != expected_samples:
        raise RuntimeError("Prefix execution scope mismatch")
    if int(execution.get("phylo_parallel_workers", 1)) != 1:
        raise RuntimeError("Prefix must be the pinned serial execution")
    if source_lock.get("schema_name") != SOURCE_LOCK_SCHEMA or source_lock.get("pass") is not True:
        raise RuntimeError("Prefix source-lock receipt is not passing")
    if source_lock.get("scope") != expected_samples:
        raise RuntimeError("Prefix source-lock scope mismatch")
    if source_lock.get("child_run_manifest") != A.artifact(receipt_path):
        raise RuntimeError("Prefix source-lock child receipt mismatch")
    if source_lock.get("source_identity_before") != source_lock.get("source_identity_after"):
        raise RuntimeError("Prefix source identities changed during execution")
    checks = source_lock.get("checks") or {}
    if not checks or not all(value is True for value in checks.values()):
        raise RuntimeError("Prefix source-lock checks are not all passing")
    return site_path, assignment_path, summary, receipt, source_lock


def validate_equivalence(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_name") != EQUIVALENCE_SCHEMA or payload.get("pass") is not True:
        raise RuntimeError("Serial/parallel equivalence receipt is not passing")
    if payload.get("scope") != "algorithm_and_real_nested_high_read_fixture":
        raise RuntimeError("Serial/parallel equivalence scope mismatch")
    checks = payload.get("checks") or {}
    required = {
        "pinned_analyzer_sha256_exact",
        "source_identity_unchanged",
        "synthetic_full_payload_exact",
        "real_nested_full_payload_exact",
        "real_nested_parallel_triggered",
    }
    if set(checks) != required or not all(checks.values()):
        raise RuntimeError("Serial/parallel equivalence checks are incomplete")
    real = payload.get("real_fixture") or {}
    if real.get("full_row_and_assignment_payload_exact") is not True:
        raise RuntimeError("Real nested equivalence payload is not exact")
    if int(real.get("n_alt_after_peel", 0)) < int(real.get("parallel_min_reads", 1)):
        raise RuntimeError("Real equivalence fixture did not trigger parallel dispatch")
    source = ((payload.get("inputs") or {}).get("source_identity_before") or {}).get(
        "analyzer"
    )
    validate_current_artifact(source, "equivalence analyzer")
    return payload


def validate_replacement(
    replacement_dir: Path,
    manifest_path: Path,
    replacement_sample: str,
    expected_sites: int,
) -> tuple[Path, Path, dict[str, Any], dict[str, Any]]:
    site_path = replacement_dir / "all_ssnv_site_results.tsv.gz"
    assignment_path = replacement_dir / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
    summary_path = replacement_dir / "all_ssnv_summary.json"
    receipt_path = replacement_dir / "run_manifest.json"
    for path in (site_path, assignment_path, summary_path, receipt_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if summary.get("schema_name") != "intersubmod.all_ssnv_focal_alt_multigroup_screen":
        raise RuntimeError("Replacement summary schema mismatch")
    if receipt.get("schema_name") != "intersubmod.all_ssnv_focal_alt_multigroup_run_manifest":
        raise RuntimeError("Replacement run-manifest schema mismatch")
    if summary.get("pass") is not True or receipt.get("pass") is not True:
        raise RuntimeError("Replacement run is not passing")
    scope = summary.get("scope") or {}
    if scope.get("full_469849") is not False:
        raise RuntimeError("Replacement summary must be an explicit single-sample run")
    if scope.get("selected_samples") != [replacement_sample]:
        raise RuntimeError("Replacement summary selected sample mismatch")
    if int(scope.get("processed_sites", -1)) != expected_sites:
        raise RuntimeError("Replacement summary site count mismatch")
    outputs = receipt.get("outputs") or {}
    expected_outputs = {
        "site_results": site_path,
        "stable_assignments": assignment_path,
        "summary": summary_path,
    }
    for name, path in expected_outputs.items():
        if outputs.get(name) != A.artifact(path):
            raise RuntimeError(f"Replacement receipt artifact mismatch: {name}")
    if receipt.get("input_manifest") != A.artifact(manifest_path):
        raise RuntimeError("Replacement receipt input manifest mismatch")
    execution = receipt.get("execution") or {}
    if execution.get("selected_samples") != [replacement_sample]:
        raise RuntimeError("Replacement execution scope mismatch")
    if int(execution.get("phylo_parallel_workers", 1)) <= 1:
        raise RuntimeError("Replacement did not record seed-level parallel execution")
    source_code = receipt.get("source_code")
    if not isinstance(source_code, Mapping):
        raise RuntimeError("Replacement source-code inventory is missing")
    started = parse_utc(receipt.get("started_at_utc"), "replacement start")
    finished = parse_utc(receipt.get("finished_at_utc"), "replacement finish")
    if started > finished:
        raise RuntimeError("Replacement execution time window is reversed")
    source_paths = {
        role: validate_current_artifact(
            source_code.get(role), f"replacement {role}"
        )
        for role in (
            "analyzer",
            "focal_alt_cluster_lib",
            "latest_tag_join",
            "claim_contract_v2",
        )
    }
    for role, source_path in source_paths.items():
        if datetime.fromtimestamp(source_path.stat().st_mtime, timezone.utc) > started:
            raise RuntimeError(
                f"Replacement {role} was modified after execution started"
            )
    if source_paths["analyzer"].stat().st_mode & 0o222:
        raise RuntimeError("Replacement analyzer is not a read-only pinned source file")
    return site_path, assignment_path, summary, receipt


def build_summary(
    manifest: Mapping[str, Any],
    pooled: A.ScreenSummary,
    by_sample: Mapping[str, A.ScreenSummary],
    by_truth: Mapping[str, A.ScreenSummary],
    by_biological: Mapping[str, A.ScreenSummary],
    by_branch: Mapping[str, A.ScreenSummary],
    assignments: int,
) -> dict[str, Any]:
    tag_audit = pooled.latest_tag_payload()
    expected_sites = int(manifest["totals"]["all_ssnv"])
    return {
        "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_screen",
        "schema_version": A.OUTPUT_SCHEMA_VERSION,
        "created_at_utc": A.now_utc(),
        "status": "EXECUTION_PASS",
        "pass_semantics": "execution_integrity_only_not_scientific_confirmation",
        "scope": {
            "full_469849": True,
            "selected_datasets": A.DATASETS,
            "selected_samples": A.DATASETS,
            "expected_sites": expected_sites,
            "processed_sites": pooled.n_sites,
        },
        "population": "LongPhase-S recalibrated FILTER=PASS chr1-22 biallelic sSNVs",
        "clustering_contract": {
            "truth_blind": True,
            "cooccurrence_blind": True,
            "distance_source": "exact existing C++ BERNOULLI matrices",
            "screen": A.SCREEN_CONTRACT,
            "screen_contract_semantics": "legacy_algorithm_identity_for_downstream_compatibility",
            "m1_stability_gate_contract": A.M1_STABILITY_GATE_CONTRACT,
            "prior_screen_thresholds": A.PRIOR_SCREEN_THRESHOLDS,
            "global_fdr_calibrated": False,
            "strict_methyl_partition_robustness_status": "NOT_EVALUATED_AT_M1_SCREEN",
            "strict_confirm_status_legacy_alias": "NOT_RUN",
            "stable_null_multigroup_basis": A.M1_STABILITY_GATE_CONTRACT,
            "strict_confirm_candidate_legacy_alias_basis": "same_as_stable_null_multigroup",
            "strict_confirm_candidate_is_formal_r1_claim": False,
        },
        "pooled_site_weighted": pooled.payload(),
        "per_dataset": {sample: by_sample[sample].payload() for sample in A.DATASETS},
        "per_sample": {sample: by_sample[sample].payload() for sample in A.DATASETS},
        "posthoc_truth_strata": {label: by_truth[label].payload() for label in sorted(A.TRUTH_LABELS)},
        "posthoc_biological_id_strata": {
            biological: by_biological[biological].payload() for biological in sorted(by_biological)
        },
        "posthoc_ledger_branch_strata": {
            branch: by_branch[branch].payload() for branch in sorted(by_branch)
        },
        "n_stable_assignment_records": assignments,
        "latest_hp_ps_terminal_join_audit": tag_audit,
        "interpretation_guardrail": (
            "This is an M1 screen, not a genome-wide FDR procedure. A focal-ALT stable methyl "
            "multigroup supports ALT-carrier read-level methyl heterogeneity only, not a genetic clone, "
            "cell fraction, or lineage topology. R1 strict methyl-partition robustness is not evaluated "
            "at this stage; not evaluated is not failure."
        ),
        "recovery_merge": {
            "replacement_sample": "HCC1954",
            "method_change": False,
            "execution_change": "independent_seed_runs_parallelized_only",
        },
        "pass": pooled.n_sites == expected_sites and assignments == pooled.stable and tag_audit["pass"],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--intersubmod-root", type=Path, required=True)
    parser.add_argument("--prefix-dir", type=Path, required=True)
    parser.add_argument("--replacement-dir", type=Path, required=True)
    parser.add_argument("--replacement-sample", choices=A.DATASETS, required=True)
    parser.add_argument("--equivalence-receipt", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_at = A.now_utc()
    source_code_start = {
        "recovery_merger": A.artifact(Path(__file__)),
        "current_analyzer": A.artifact(Path(A.__file__)),
        "focal_alt_cluster_lib": A.artifact(Path(A.F.__file__)),
        "latest_tag_join": A.artifact(Path(A.TAGS.__file__)),
        "claim_contract_v2": A.artifact(A.TOPIC_ROOT / "claim-contract-v2.md"),
    }
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    A.validate_manifest(manifest)
    if args.replacement_sample != A.DATASETS[-1]:
        raise RuntimeError("Recovery merge currently requires replacement of the final canonical sample")
    expected_by_sample = {
        entry["sample"]: int(entry["counts"]["all_ssnv"]) for entry in manifest["samples"]
    }
    prefix_samples = A.DATASETS[:-1]
    prefix_expected_sites = sum(expected_by_sample[sample] for sample in prefix_samples)
    prefix_sites, prefix_assignments, prefix_summary, prefix_receipt, prefix_source_lock = (
        validate_prefix(
            args.prefix_dir,
            args.manifest,
            prefix_samples,
            prefix_expected_sites,
        )
    )
    replacement_sites, replacement_assignments, replacement_summary, replacement_receipt = (
        validate_replacement(
            args.replacement_dir,
            args.manifest,
            args.replacement_sample,
            expected_by_sample[args.replacement_sample],
        )
    )
    equivalence = validate_equivalence(args.equivalence_receipt)
    prefix_sources = prefix_receipt.get("source_code") or {}
    replacement_sources = replacement_receipt.get("source_code") or {}
    equivalence_analyzer = (
        ((equivalence.get("inputs") or {}).get("source_identity_before") or {}).get(
            "analyzer"
        )
        or {}
    )
    analyzer_hashes = {
        str((prefix_sources.get("analyzer") or {}).get("sha256", "")),
        str((replacement_sources.get("analyzer") or {}).get("sha256", "")),
        str(equivalence_analyzer.get("sha256", "")),
    }
    if len(analyzer_hashes) != 1 or "" in analyzer_hashes:
        raise RuntimeError("Prefix/replacement/equivalence analyzer SHA-256 mismatch")
    for role in ("focal_alt_cluster_lib", "latest_tag_join", "claim_contract_v2"):
        if prefix_sources.get(role) != replacement_sources.get(role):
            raise RuntimeError(f"Prefix/replacement source identity mismatch: {role}")
    completion_receipts = {
        entry["sample"]: A.validate_completed_sample_run(args.intersubmod_root / entry["sample"], entry)
        for entry in manifest["samples"]
    }
    A.create_output_dir(args.output_dir)

    site_output = args.output_dir / "all_ssnv_site_results.tsv.gz"
    assignment_output = args.output_dir / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
    summary_output = args.output_dir / "all_ssnv_summary.json"
    receipt_output = args.output_dir / "run_manifest.json"
    pooled = A.ScreenSummary()
    by_sample = {sample: A.ScreenSummary() for sample in A.DATASETS}
    by_truth = {label: A.ScreenSummary() for label in sorted(A.TRUTH_LABELS)}
    biological_ids = sorted({entry["biological_id"] for entry in manifest["samples"]})
    by_biological = {biological: A.ScreenSummary() for biological in biological_ids}
    by_branch: dict[str, A.ScreenSummary] = {}
    sample_counts: Counter[str] = Counter()
    seen_sites: set[tuple[str, str, int, str, str]] = set()
    stable_sites: set[tuple[str, str, int, str, str]] = set()
    assignment_keys: set[tuple[str, str, int, str, str]] = set()
    prefix_site_digest = hashlib.sha256()
    prefix_assignment_digest = hashlib.sha256()

    def consume_site(row: dict[str, str], writer: csv.DictWriter, *, prefix: bool) -> None:
        key = site_key(row)
        if key in seen_sites:
            raise RuntimeError(f"Duplicate merged site key: {key}")
        seen_sites.add(key)
        sample_counts[row["sample"]] += 1
        normalized = summary_row(row)
        pooled.add(normalized)
        by_sample[row["sample"]].add(normalized)
        by_truth[row["truth_label"]].add(normalized)
        by_biological[row["biological_id"]].add(normalized)
        branch = row.get("ssnv_branch") or "NA"
        by_branch.setdefault(branch, A.ScreenSummary()).add(normalized)
        if normalized["stable_null_multigroup"]:
            stable_sites.add(key)
        writer.writerow(row)
        if prefix:
            logical_digest_update(prefix_site_digest, "\t".join(row[field] for field in A.SITE_FIELDS))

    with gzip.open(site_output, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, A.SITE_FIELDS, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        for row in iter_prefix_sites(prefix_sites, args.replacement_sample):
            consume_site(row, writer, prefix=True)
        expected_prefix = sum(expected_by_sample[sample] for sample in A.DATASETS[:-1])
        if len(seen_sites) != expected_prefix:
            raise RuntimeError(f"Canonical prefix count mismatch: {len(seen_sites)} != {expected_prefix}")
        if any(sample_counts[sample] != expected_by_sample[sample] for sample in A.DATASETS[:-1]):
            raise RuntimeError("Canonical prefix per-sample counts are incomplete")
        for row in iter_complete_sites(replacement_sites, args.replacement_sample):
            consume_site(row, writer, prefix=False)

    def consume_assignment(row: dict[str, Any], handle: Any, *, prefix: bool) -> None:
        key = assignment_key(row)
        if key in assignment_keys:
            raise RuntimeError(f"Duplicate merged assignment key: {key}")
        assignment_keys.add(key)
        encoded = json.dumps(row, ensure_ascii=False, separators=(",", ":"))
        handle.write(encoded + "\n")
        if prefix:
            logical_digest_update(prefix_assignment_digest, encoded)

    with gzip.open(assignment_output, "wt", encoding="utf-8") as handle:
        for row in iter_prefix_assignments(prefix_assignments, args.replacement_sample):
            consume_assignment(row, handle, prefix=True)
        for row in iter_complete_assignments(replacement_assignments, args.replacement_sample):
            consume_assignment(row, handle, prefix=False)

    expected_total = int(manifest["totals"]["all_ssnv"])
    if len(seen_sites) != expected_total or sample_counts != Counter(expected_by_sample):
        raise RuntimeError("Merged site scope does not match the frozen manifest")
    if assignment_keys != stable_sites:
        raise RuntimeError(
            "Merged stable/assignment key-set mismatch: "
            f"stable={len(stable_sites)} assignments={len(assignment_keys)}"
        )
    summary = build_summary(
        manifest,
        pooled,
        by_sample,
        by_truth,
        by_biological,
        by_branch,
        len(assignment_keys),
    )
    summary["recovery_merge"].update(
        {
            "prefix_source_locked": True,
            "replacement_read_only_pinned_source": True,
            "serial_parallel_exact_equivalence_pass": True,
            "pinned_analyzer_sha256": next(iter(analyzer_hashes)),
        }
    )
    if summary["pass"] is not True:
        raise RuntimeError("Merged summary execution integrity failed")
    source_code_end = {
        "recovery_merger": A.artifact(Path(__file__)),
        "current_analyzer": A.artifact(Path(A.__file__)),
        "focal_alt_cluster_lib": A.artifact(Path(A.F.__file__)),
        "latest_tag_join": A.artifact(Path(A.TAGS.__file__)),
        "claim_contract_v2": A.artifact(A.TOPIC_ROOT / "claim-contract-v2.md"),
    }
    if source_code_start != source_code_end:
        raise RuntimeError("Recovery merger source identity changed during execution")
    summary_output.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    finished_at = A.now_utc()
    receipt = {
        "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_run_manifest",
        "schema_version": A.OUTPUT_SCHEMA_VERSION,
        "created_at_utc": finished_at,
        "status": "EXECUTION_PASS",
        "pass_semantics": "execution_integrity_only_not_scientific_confirmation",
        "started_at_utc": started_at,
        "finished_at_utc": finished_at,
        "command": sys.argv,
        "input_manifest": A.artifact(args.manifest),
        "intersubmod_root": str(args.intersubmod_root.resolve()),
        "completed_dataset_runs": completion_receipts,
        "completed_sample_runs": completion_receipts,
        "source_code": source_code_end,
        "execution": {
            "workers": 0,
            "chunk_size": 0,
            "max_pending_chunks": 0,
            "per_site_future_submission": False,
            "selected_datasets": A.DATASETS,
            "selected_samples": A.DATASETS,
            "recovery_merge": True,
        },
        "recovery": {
            "schema_name": SCHEMA_NAME,
            "replacement_sample": args.replacement_sample,
            "prefix_samples": prefix_samples,
            "prefix_site_logical_sha256": prefix_site_digest.hexdigest(),
            "prefix_assignment_logical_sha256": prefix_assignment_digest.hexdigest(),
            "prefix_site_count": sum(sample_counts[sample] for sample in A.DATASETS[:-1]),
            "prefix_source_artifacts": {
                "sites": A.artifact(prefix_sites),
                "assignments": A.artifact(prefix_assignments),
                "summary": A.artifact(args.prefix_dir / "all_ssnv_summary.json"),
                "receipt": A.artifact(args.prefix_dir / "run_manifest.json"),
                "source_lock_receipt": A.artifact(
                    args.prefix_dir / "source_lock_receipt.json"
                ),
            },
            "prefix_summary_scope": prefix_summary["scope"],
            "prefix_execution": prefix_receipt["execution"],
            "prefix_source_lock_checks": prefix_source_lock["checks"],
            "prefix_source_files_complete_closed_run": True,
            "prefix_read_reached_eof": True,
            "replacement_summary": A.artifact(args.replacement_dir / "all_ssnv_summary.json"),
            "replacement_receipt": A.artifact(args.replacement_dir / "run_manifest.json"),
            "replacement_execution": replacement_receipt["execution"],
            "replacement_summary_scope": replacement_summary["scope"],
            "method_change": False,
            "execution_change": "independent_seed_runs_parallelized_only",
            "serial_parallel_exact_equivalence_required": True,
            "serial_parallel_exact_equivalence_receipt": A.artifact(
                args.equivalence_receipt
            ),
            "serial_parallel_exact_equivalence_checks": equivalence["checks"],
            "pinned_analyzer_sha256": next(iter(analyzer_hashes)),
            "prefix_replacement_source_dependencies_exact": True,
            "recovery_source_identity_unchanged_during_merge": True,
        },
        "contracts": {
            "truth_and_cooccurrence_enter_clustering": False,
            "latest_hp_ps_join": "same-run sidecar projected join before ALT selection; missing/conflict hard fail",
            "screen_global_fdr_calibrated": False,
            "screen_contract_semantics": "legacy_algorithm_identity_for_downstream_compatibility",
            "m1_stability_gate_contract": A.M1_STABILITY_GATE_CONTRACT,
            "strict_methyl_partition_robustness_status": "NOT_EVALUATED_AT_M1_SCREEN",
            "strict_confirm_status_legacy_alias": "NOT_RUN",
            "stable_null_multigroup_basis": A.M1_STABILITY_GATE_CONTRACT,
            "strict_confirm_candidate_legacy_alias_basis": "same_as_stable_null_multigroup",
            "strict_confirm_candidate_is_formal_r1_claim": False,
            "prior_screen_thresholds": A.PRIOR_SCREEN_THRESHOLDS,
            "existing_results_overwritten": False,
            "latest_hp_ps_terminal_join": pooled.latest_tag_payload(),
        },
        "outputs": {
            "site_results": A.artifact(site_output),
            "stable_assignments": A.artifact(assignment_output),
            "summary": A.artifact(summary_output),
        },
        "counts": {
            "expected_sites": expected_total,
            "processed_sites": pooled.n_sites,
            "stable_assignment_records": len(assignment_keys),
            "reads_tsv_site_rows": pooled.reads_tsv_site_rows,
            "exact_hp_ps_site_read_joins": pooled.latest_tag_reads_joined,
            "ps_present_site_read_joins": pooled.latest_tag_ps_present,
            "source_hp_replaced_site_read_joins": pooled.source_hp_replaced_reads,
        },
        "pass": True,
    }
    receipt_output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output_dir": str(args.output_dir.resolve()),
                "processed_sites": pooled.n_sites,
                "stable_assignments": len(assignment_keys),
                "prefix_site_logical_sha256": prefix_site_digest.hexdigest(),
                "pass": True,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
