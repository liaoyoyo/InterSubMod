#!/usr/bin/env python3
"""Fail-closed validation for the 2026-08-13 research handoff package."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from typing import Any, Callable
from urllib.parse import unquote, urlsplit


PACKAGE_RELATIVE_PATH = Path("docs/handoff/20260813_完整研究資料與軟體交接_01")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
MARKDOWN_LINK_RE = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
EXTERNAL_URL_RE = re.compile(r"(?:https?:)?//", re.IGNORECASE)
EVIDENCE_STATUS_ENUM = {
    "AUTHORITY",
    "VALIDATED_DERIVED",
    "PARTIAL",
    "HISTORICAL",
    "INVALIDATED",
    "IN_PROGRESS",
}
SCOPE_ENUM = {"FULL", "PARTIAL", "DEMO"}
FINALITY_ENUM = {"FINAL_FOR_SCOPE", "NON_FINAL", "SUPERSEDED"}
REGENERATION_SEMANTIC_PREFIXES = ("REPLAY_ONLY:", "VERIFY_ONLY:", "NOT_REGENERABLE_FROM_HANDOFF:")


class ContractError(ValueError):
    """Raised when a handoff contract is not satisfied."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def declared_source_path(raw_source: object) -> tuple[Path | None, str]:
    """Resolve an evidence source without pretending unavailable private repos exist.

    InterSubMod sources are always required to resolve against this checkout.  A
    LongLineage source is verified only when the caller explicitly supplies the
    clean preview checkout through ``INTERSUBMOD_HANDOFF_LONGLINEAGE_ROOT``.
    """

    require(isinstance(raw_source, str) and raw_source, "evidence source_path must be non-empty")
    if raw_source.startswith("InterSubMod/"):
        return repository_root() / raw_source.removeprefix("InterSubMod/"), "REQUIRED_LOCAL_SOURCE"
    if raw_source.startswith("LongLineage/"):
        external_root = os.environ.get("INTERSUBMOD_HANDOFF_LONGLINEAGE_ROOT")
        if external_root:
            return Path(external_root).expanduser().resolve() / raw_source.removeprefix("LongLineage/"), "EXPLICIT_EXTERNAL_SOURCE"
        return None, "EXTERNAL_SOURCE_NOT_MOUNTED"
    return None, "NON_FILE_SOURCE_DESCRIPTION"


def load_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def package_path(package_root: Path, relative: str, label: str, *, base: Path | None = None) -> Path:
    require(isinstance(relative, str) and relative.strip() != "", f"{label} must be a non-empty path")
    candidate = ((base or package_root) / relative).resolve()
    try:
        candidate.relative_to(package_root.resolve())
    except ValueError as error:
        raise ContractError(f"{label} escapes the package: {relative}") from error
    return candidate


def check_markdown_links(package_root: Path) -> dict[str, Any]:
    markdown_files = sorted(package_root.rglob("*.md"))
    require(markdown_files, "package contains no Markdown files")
    checked = 0
    missing: list[str] = []
    for markdown in markdown_files:
        text = markdown.read_text(encoding="utf-8")
        for match in MARKDOWN_LINK_RE.finditer(text):
            raw_target = match.group(1).strip()
            if raw_target.startswith("<") and raw_target.endswith(">"):
                raw_target = raw_target[1:-1].strip()
            if not raw_target or raw_target.startswith("#") or re.match(r"https?://", raw_target, re.IGNORECASE):
                continue
            parsed = urlsplit(raw_target)
            if parsed.scheme or raw_target.startswith("/"):
                continue
            target_text = unquote(parsed.path)
            if not target_text:
                continue
            target = (markdown.parent / target_text).resolve()
            checked += 1
            if not target.exists():
                missing.append(f"{markdown.relative_to(package_root)} -> {raw_target}")
    require(not missing, "missing relative Markdown links: " + "; ".join(missing))
    return {"markdown_files": len(markdown_files), "relative_links_checked": checked, "missing": 0}


def check_json_parse(package_root: Path) -> dict[str, Any]:
    json_files = sorted(package_root.rglob("*.json"))
    require(json_files, "package contains no JSON files")
    errors: list[str] = []
    for path in json_files:
        try:
            load_json(path)
        except (OSError, UnicodeError, json.JSONDecodeError) as error:
            errors.append(f"{path.relative_to(package_root)}: {error}")
    require(not errors, "unparseable JSON: " + "; ".join(errors))
    return {"json_files": len(json_files), "parse_errors": 0}


def check_evidence_manifest(package_root: Path) -> dict[str, Any]:
    manifest = load_json(package_root / "evidence/EVIDENCE_MANIFEST.json")
    records = manifest.get("records")
    require(isinstance(records, list) and records, "EVIDENCE_MANIFEST records must be a non-empty list")
    evidence_ids = [row.get("evidence_id") for row in records if isinstance(row, dict)]
    require(len(evidence_ids) == len(records), "EVIDENCE_MANIFEST contains a non-object record")
    require(all(isinstance(value, str) and value for value in evidence_ids), "evidence_id must be non-empty")
    require(len(evidence_ids) == len(set(evidence_ids)), "EVIDENCE_MANIFEST evidence_id values are not unique")

    sha256_bound_count = 0
    summary_count = 0
    source_equal_count = 0
    external_source_not_mounted = 0
    size_bound_count = 0
    required_size_bound = {"public_claim_validation_20260813", "full_claim_registry_20260813"}
    for row in records:
        evidence_id = row["evidence_id"]
        if "evidence_status" in row:
            require(
                row["evidence_status"] in EVIDENCE_STATUS_ENUM,
                f"invalid evidence_status enum for {evidence_id}: {row['evidence_status']!r}",
            )
        if "scope" in row:
            require(row["scope"] in SCOPE_ENUM, f"invalid scope enum for {evidence_id}: {row['scope']!r}")
        if "finality" in row:
            require(
                row["finality"] in FINALITY_ENUM,
                f"invalid finality enum for {evidence_id}: {row['finality']!r}",
            )
        path = package_path(package_root, row.get("path"), f"evidence {evidence_id} path")
        require(path.is_file(), f"evidence path does not exist: {row.get('path')}")
        if "size_bytes" in row or evidence_id in required_size_bound:
            expected_size = row.get("size_bytes")
            require(
                isinstance(expected_size, int) and not isinstance(expected_size, bool) and expected_size >= 0,
                f"invalid size_bytes for {evidence_id}",
            )
            require(path.stat().st_size == expected_size, f"evidence size mismatch: {evidence_id}")
            size_bound_count += 1
        copy_status = row.get("copy_status")
        if copy_status in {"EXACT_COPY", "NEW_PACKAGE_RECEIPT", "PACKAGE_RELINKED_COPY"}:
            expected = row.get("sha256")
            require(isinstance(expected, str) and SHA256_RE.fullmatch(expected), f"invalid sha256 for {evidence_id}")
            require(sha256_file(path) == expected, f"evidence sha256 mismatch: {evidence_id}")
            if copy_status == "EXACT_COPY":
                source, source_status = declared_source_path(row.get("source_path"))
                if source is not None:
                    require(source.is_file(), f"declared exact-copy source is missing: {row.get('source_path')}")
                    source_hash = sha256_file(source)
                    require(source_hash == expected, f"exact-copy source hash drift: {evidence_id}")
                    require(source.read_bytes() == path.read_bytes(), f"exact-copy bytes differ from source: {evidence_id}")
                    source_equal_count += 1
                elif source_status == "EXTERNAL_SOURCE_NOT_MOUNTED":
                    external_source_not_mounted += 1
            if copy_status == "PACKAGE_RELINKED_COPY":
                source_sha = row.get("source_sha256")
                require(
                    isinstance(source_sha, str) and SHA256_RE.fullmatch(source_sha),
                    f"invalid original source_sha256 for relinked evidence {evidence_id}",
                )
            sha256_bound_count += 1
        elif copy_status == "SUMMARY_HASH_BOUND":
            source = package_path(package_root, row.get("source_path"), f"evidence {evidence_id} source_path")
            require(source.is_file(), f"summary source path does not exist: {row.get('source_path')}")
            expected = row.get("source_sha256")
            require(
                isinstance(expected, str) and SHA256_RE.fullmatch(expected),
                f"invalid source_sha256 for {evidence_id}",
            )
            require(sha256_file(source) == expected, f"summary source_sha256 mismatch: {evidence_id}")
            summary_count += 1
        else:
            raise ContractError(f"unsupported evidence copy_status for {evidence_id}: {copy_status!r}")
    by_id = {row["evidence_id"]: row for row in records}
    sampled_contracts = {
        "tagged_bam_sampled_identity_baseline_20260711": ("HISTORICAL", "PARTIAL", "NON_FINAL"),
        "tagged_bam_sampled_replay_manifest_20260813": ("VALIDATED_DERIVED", "PARTIAL", "NON_FINAL"),
        "tagged_bam_sampled_identity_replay_20260813": ("VALIDATED_DERIVED", "PARTIAL", "NON_FINAL"),
    }
    for evidence_id, expected in sampled_contracts.items():
        require(evidence_id in by_id, f"missing sampled tagged BAM evidence record: {evidence_id}")
        actual = tuple(by_id[evidence_id].get(key) for key in ("evidence_status", "scope", "finality"))
        require(actual == expected, f"sampled tagged BAM evidence contract drifted for {evidence_id}: {actual!r}")
    reader_id = "reader_acceptance_20260813"
    require(reader_id in by_id, "fresh-reader acceptance evidence record is missing")
    reader_contract = tuple(by_id[reader_id].get(key) for key in ("evidence_status", "scope", "finality"))
    require(
        reader_contract == ("VALIDATED_DERIVED", "FULL", "FINAL_FOR_SCOPE"),
        f"fresh-reader acceptance evidence contract drifted: {reader_contract!r}",
    )
    return {
        "records": len(records),
        "sha256_bound_records": sha256_bound_count,
        "summary_hash_bound": summary_count,
        "exact_copy_source_equal": source_equal_count,
        "external_exact_copy_sources_not_mounted": external_source_not_mounted,
        "size_bound_records": size_bound_count,
        "missing": 0,
        "hash_mismatch": 0,
        "reader_acceptance_records": 1,
    }


def check_dataset_registry(package_root: Path) -> dict[str, Any]:
    records = load_json(package_root / "registries/dataset_registry.json")
    require(isinstance(records, list), "dataset_registry must be a JSON array")
    technical_ids = [row.get("technical_dataset_id") for row in records if isinstance(row, dict)]
    require(len(technical_ids) == len(records), "dataset_registry contains a non-object record")
    require(len(records) == 7, f"expected 7 technical datasets, found {len(records)}")
    require(all(isinstance(value, str) and value for value in technical_ids), "technical_dataset_id must be non-empty")
    require(len(technical_ids) == len(set(technical_ids)), "technical_dataset_id values are not unique")
    biological_ids = {row.get("biological_id") for row in records}
    require(None not in biological_ids and "" not in biological_ids, "biological_id must be non-empty")
    require(len(biological_ids) == 6, f"expected 6 biological samples, found {len(biological_ids)}")
    return {"technical_datasets": len(records), "biological_samples": len(biological_ids)}


def check_run_registry(package_root: Path) -> dict[str, Any]:
    registry = load_json(package_root / "registries/run_registry.json")
    records = registry.get("records")
    require(isinstance(records, list), "run_registry records must be a list")
    require(len(records) == 51, f"expected 51 physical runs, found {len(records)}")
    physical_ids = [row.get("physical_run_id") for row in records if isinstance(row, dict)]
    physical_paths = [row.get("physical_path") for row in records if isinstance(row, dict)]
    require(len(physical_ids) == len(records), "run_registry contains a non-object record")
    require(all(isinstance(value, str) and value for value in physical_ids), "physical_run_id must be non-empty")
    require(all(isinstance(value, str) and value for value in physical_paths), "physical_path must be non-empty")
    require(len(set(physical_ids)) == 51, "physical_run_id values are not unique")
    require(len(set(physical_paths)) == 51, "physical_path values are not unique")
    reconciliation = registry.get("reconciliation")
    require(isinstance(reconciliation, dict), "run_registry reconciliation must be an object")
    require(reconciliation.get("physical_directories_total") == 51, "run reconciliation total is not 51")
    require(
        reconciliation.get("current_physical_directories", 0)
        + reconciliation.get("pending_archive_physical_directories", 0)
        == 51,
        "run reconciliation current + pending archive does not equal 51",
    )
    return {"physical_runs": 51, "unique_physical_run_ids": 51, "unique_physical_paths": 51}


def check_artifact_registry(package_root: Path) -> dict[str, Any]:
    records = load_json(package_root / "registries/artifact_registry.json")
    require(isinstance(records, list), "artifact_registry must be a JSON array")
    require(len(records) == 36, f"expected 36 artifacts, found {len(records)}")
    artifact_ids = [row.get("artifact_id") for row in records if isinstance(row, dict)]
    require(len(artifact_ids) == len(records), "artifact_registry contains a non-object record")
    require(all(isinstance(value, str) and value for value in artifact_ids), "artifact_id must be non-empty")
    require(len(set(artifact_ids)) == 36, "artifact_id values are not unique")

    final = [row for row in records if row.get("finality") == "FINAL_FOR_SCOPE"]
    require(len(final) == 20, f"expected 20 FINAL_FOR_SCOPE artifacts, found {len(final)}")
    final_status = Counter(row.get("evidence_status") for row in final)
    require(
        final_status == Counter({"AUTHORITY": 19, "VALIDATED_DERIVED": 1}),
        f"unexpected FINAL_FOR_SCOPE evidence statuses: {dict(final_status)}",
    )
    provenance = [row for row in final if row.get("evidence_status") == "VALIDATED_DERIVED"]
    require(
        len(provenance) == 1 and provenance[0].get("artifact_id") == "adjudication.frozen_binary_source",
        "the sole provenance-only FINAL_FOR_SCOPE artifact is not frozen binary source adjudication",
    )
    for row in final:
        artifact_id = row["artifact_id"]
        require(isinstance(row.get("producer"), str) and row["producer"].strip(), f"final artifact lacks producer: {artifact_id}")
        require(
            isinstance(row.get("sha256"), str) and SHA256_RE.fullmatch(row["sha256"]),
            f"final artifact lacks valid sha256: {artifact_id}",
        )
        require(
            isinstance(row.get("producer_commit"), str) and row["producer_commit"].strip(),
            f"final artifact lacks producer_commit provenance: {artifact_id}",
        )
        require(
            isinstance(row.get("inputs"), list) and len(row["inputs"]) > 0,
            f"final artifact lacks input provenance: {artifact_id}",
        )
        command = row.get("regeneration_command")
        require(
            isinstance(command, str) and command.startswith(REGENERATION_SEMANTIC_PREFIXES),
            f"final artifact lacks typed replay/regeneration semantics: {artifact_id}",
        )

    tagged_bams = [row for row in records if row.get("artifact_type") == "longphase_s_tagged_bam"]
    require(len(tagged_bams) == 14, f"expected 14 tagged BAM artifacts, found {len(tagged_bams)}")
    for row in tagged_bams:
        require(row.get("scope") == "PARTIAL", f"tagged BAM is not PARTIAL: {row.get('artifact_id')}")
        require(row.get("evidence_status") == "PARTIAL", f"tagged BAM evidence is not PARTIAL: {row.get('artifact_id')}")
        require(row.get("finality") == "NON_FINAL", f"tagged BAM is not NON_FINAL: {row.get('artifact_id')}")

    paired_full = [row for row in tagged_bams if row.get("mode") == "paired_full"]
    paired_pileup = [row for row in tagged_bams if row.get("mode") == "paired_pileup"]
    require(len(paired_full) == 7 and len(paired_pileup) == 7, "tagged BAM inventory is not 7 paired_full + 7 paired_pileup")
    require(
        sum(row.get("size_bytes", -1) for row in paired_full) == 1_840_983_466_353,
        "paired_full tagged BAM byte total drifted",
    )
    require(
        sum(row.get("size_bytes", -1) for row in tagged_bams) == 3_709_322_840_333,
        "14-object tagged BAM byte total drifted",
    )

    baseline_path = package_root / "evidence/canonical_longphase_tagged_bam_sampled_identity_v1.json"
    replay_manifest_path = package_root / "evidence/tagged_bam_sampled_replay_manifest.json"
    replay_path = package_root / "evidence/tagged_bam_sampled_identity_replay_20260813.json"
    baseline = load_json(baseline_path)
    replay_manifest = load_json(replay_manifest_path)
    replay = load_json(replay_path)
    require(
        baseline.get("schema_name") == "intersubmod.canonical_longphase_bam_immutability"
        and baseline.get("schema_version") == "1.0.0"
        and baseline.get("dataset_count") == 7,
        "tagged BAM sampled baseline contract mismatch",
    )
    require(
        baseline.get("identity_set_sha256")
        == "ce6c63d42e3f334d6847a1a6d3e46ead165b59a03197acb098319be5c67fcf90",
        "tagged BAM sampled baseline set identity drifted",
    )
    baseline_samples = baseline.get("samples")
    require(isinstance(baseline_samples, list) and len(baseline_samples) == 7, "sampled baseline must contain seven rows")
    artifacts_by_dataset = {row.get("technical_dataset_id"): row for row in paired_full}
    require(len(artifacts_by_dataset) == 7, "paired_full technical dataset IDs are not unique")
    for sample_row in baseline_samples:
        require(isinstance(sample_row, dict), "sampled baseline contains a non-object row")
        sample = sample_row.get("sample")
        identity = sample_row.get("identity")
        require(sample in artifacts_by_dataset and isinstance(identity, dict), f"sampled baseline dataset is not registered: {sample!r}")
        artifact = artifacts_by_dataset[sample]
        locations = artifact.get("machine_locations")
        require(isinstance(locations, list) and len(locations) == 1, f"paired_full location contract mismatch: {sample}")
        require(identity.get("requested_path") == locations[0].get("path"), f"sampled baseline path mismatch: {sample}")
        require(identity.get("size_bytes") == artifact.get("size_bytes"), f"sampled baseline size mismatch: {sample}")
        chunks = identity.get("chunks")
        require(
            isinstance(chunks, list)
            and [chunk.get("label") for chunk in chunks if isinstance(chunk, dict)] == ["first", "middle", "last"]
            and all(isinstance(chunk.get("sha256"), str) and SHA256_RE.fullmatch(chunk["sha256"]) for chunk in chunks),
            f"sampled baseline chunk contract mismatch: {sample}",
        )

    require(
        replay_manifest.get("schema_name") == "intersubmod.tagged_bam_sampled_replay_manifest"
        and replay_manifest.get("scope") == "HISTORICAL_PRE_FIX_PAIRED_FULL_SAMPLED_IDENTITY_ONLY"
        and replay_manifest.get("evidence_status") == "PARTIAL"
        and replay_manifest.get("finality") == "NON_FINAL",
        "tagged BAM replay manifest overstates its scope or finality",
    )
    manifest_samples = replay_manifest.get("samples")
    require(isinstance(manifest_samples, list) and len(manifest_samples) == 7, "tagged BAM replay manifest must contain seven rows")
    manifest_map = {
        row.get("sample"): row.get("tumor_bam")
        for row in manifest_samples
        if isinstance(row, dict)
    }
    require(
        manifest_map == {
            dataset: artifact["machine_locations"][0]["path"]
            for dataset, artifact in artifacts_by_dataset.items()
        },
        "tagged BAM replay manifest does not exactly join the paired_full artifact registry",
    )
    require(
        replay.get("schema_name") == "intersubmod.canonical_longphase_bam_immutability_verification"
        and replay.get("schema_version") == "1.0.0"
        and replay.get("baseline_sha256") == sha256_file(baseline_path)
        and replay.get("manifest_sha256") == sha256_file(replay_manifest_path)
        and replay.get("dataset_count") == 7
        and replay.get("all_match") is True,
        "tagged BAM sampled replay contract is not a seven-object PASS",
    )
    comparisons = replay.get("comparisons")
    require(isinstance(comparisons, list) and len(comparisons) == 7, "tagged BAM replay must contain seven comparisons")
    require(
        {row.get("sample") for row in comparisons if isinstance(row, dict)} == set(artifacts_by_dataset)
        and all(
            isinstance(row, dict) and row.get("match") is True and row.get("differing_fields") == []
            for row in comparisons
        ),
        "tagged BAM replay contains a mismatch or incomplete comparison",
    )
    return {
        "artifacts": 36,
        "unique_artifact_ids": 36,
        "final_for_scope": 20,
        "science_authority_final": 19,
        "provenance_only_validated_derived_final": 1,
        "final_with_producer_commit_inputs_and_typed_replay_semantics": 20,
        "tagged_bam_partial_non_final": 14,
        "tagged_bam_inventory_bytes": 3_709_322_840_333,
        "paired_full_sampled_identity": {
            "objects": 7,
            "bytes": 1_840_983_466_353,
            "identity_set_sha256": baseline["identity_set_sha256"],
            "replay_match": 7,
            "full_file_sha256": "NOT_COMPUTED",
            "evidence_status": "PARTIAL",
            "finality": "NON_FINAL",
        },
    }


def check_claim_registries(package_root: Path) -> dict[str, Any]:
    pointer = load_json(package_root / "registries/claim_registry.json")
    full_path = package_path(
        package_root,
        pointer.get("records_uri"),
        "claim records_uri",
        base=package_root / "registries",
    )
    require(full_path.is_file(), f"full claim registry does not exist: {pointer.get('records_uri')}")
    expected_hash = pointer.get("records_sha256")
    require(isinstance(expected_hash, str) and SHA256_RE.fullmatch(expected_hash), "claim registry pointer sha256 is invalid")
    require(sha256_file(full_path) == expected_hash, "full claim registry hash does not match pointer")
    full = load_json(full_path)
    claims = full.get("claims")
    require(isinstance(claims, list), "full claim registry claims must be a list")
    require(pointer.get("records_count") == 158, "embedded claim pointer records_count is not 158")
    require(full.get("counts", {}).get("claims") == 158, "full claim registry count is not 158")
    require(len(claims) == 158, f"full claim registry contains {len(claims)} claims, expected 158")
    claim_ids = [row.get("claim_id") for row in claims if isinstance(row, dict)]
    require(len(claim_ids) == len(claims), "full claim registry contains a non-object claim")
    require(len(set(claim_ids)) == 158, "full claim registry claim_id values are not unique")
    full_counts = dict(full.get("counts", {}))
    full_counts.pop("claims", None)
    require(pointer.get("counts") == full_counts, "embedded and full claim count distributions differ")
    for category, distribution in full_counts.items():
        require(isinstance(distribution, dict), f"claim count category is not an object: {category}")
        require(sum(distribution.values()) == 158, f"claim count category does not sum to 158: {category}")
    expected_gates = {
        "P0_SOURCE_READY": "PASS",
        "SOURCE_READY": "BLOCKED",
        "PUBLICATION_READY": "BLOCKED",
        "RELEASE_READY": "BLOCKED",
    }
    require(pointer.get("gates") == expected_gates, "embedded release gate statuses are not fail-closed")
    for gate, expected_status in expected_gates.items():
        require(full.get("gates", {}).get(gate, {}).get("status") == expected_status, f"full claim gate differs: {gate}")

    # The exact-copy receipt must also bind this embedded registry and the
    # current QA matrix.  A matching stale receipt is not freshness evidence.
    validation = load_json(package_root / "evidence/public_claim_validation_receipt.json")
    require(
        validation.get("schema_name") == "intersubmod.public_claim_validation_receipt"
        and validation.get("schema_version") == "2.0.0",
        "public claim validation receipt is not freshness-bound schema 2.0.0",
    )
    require(validation.get("verdict") == "PASS", "public claim validation receipt is not PASS")
    require(
        isinstance(validation.get("publication_status"), str)
        and validation["publication_status"].startswith("BLOCKED_")
        and validation.get("release_status") == "BLOCKED",
        "public claim validation receipt promotes publication or release",
    )
    require(
        validation.get("scope") == "LOCAL_SOURCE_PLUS_C108_LIVE_RECEIPT_PLUS_CHROMIUM_QA",
        "public claim validation receipt scope is stale or incomplete",
    )
    claim_contract = validation.get("claim_registry_contract")
    require(isinstance(claim_contract, dict), "public claim receipt lacks claim_registry_contract")
    require(claim_contract.get("verdict") == "PASS", "claim registry contract is not PASS")
    require(claim_contract.get("errors") == [], "claim registry contract contains errors")
    require(claim_contract.get("registry_sha256") == expected_hash, "validation receipt binds a different claim registry")
    require(claim_contract.get("public_source_files") == 34, "validation receipt does not bind 34 public sources")
    browser_contract = validation.get("browser_qa_contract")
    require(isinstance(browser_contract, dict), "public claim receipt lacks browser_qa_contract")
    require(browser_contract.get("verdict") == "PASS", "browser QA contract is not PASS")
    require(browser_contract.get("errors") == [], "browser QA contract contains errors")
    require(
        (
            browser_contract.get("html_files"),
            browser_contract.get("standalone_svg_files"),
            browser_contract.get("browser_cases"),
        )
        == (21, 21, 84),
        "validation receipt does not bind the exact 21 HTML / 21 SVG / 84 browser-case matrix",
    )
    runner = validation.get("runner")
    require(isinstance(runner, dict), "public claim validation receipt lacks runner provenance")
    python_match = re.fullmatch(r"(\d+)\.(\d+)(?:\.(\d+))?", str(runner.get("python", "")))
    require(
        python_match is not None and tuple(int(value or 0) for value in python_match.groups()) >= (3, 10, 0),
        "public claim validation receipt was not produced with Python >=3.10",
    )
    require(pointer.get("verified_at") == validation.get("created_at"), "claim pointer and validation receipt timestamps differ")
    return {
        "claims": 158,
        "unique_claim_ids": 158,
        "sha256_match": True,
        "public_source_files": 34,
        "browser_qa_matrix": {"html": 21, "standalone_svg": 21, "cases": 84},
        "validation_receipt_schema": "2.0.0",
        "release_gates": expected_gates,
    }


def check_site_profile_join(package_root: Path) -> dict[str, Any]:
    datasets = load_json(package_root / "registries/dataset_registry.json")
    by_technical_id = {row["technical_dataset_id"]: row for row in datasets}
    alias_registry = load_json(package_root / "registries/dataset_alias_registry.json")
    aliases = alias_registry.get("aliases")
    require(isinstance(aliases, list), "dataset alias registry aliases must be a list")
    alias_map: dict[str, dict[str, Any]] = {}
    for row in aliases:
        require(isinstance(row, dict), "dataset alias registry contains a non-object record")
        key = row.get("site_profile_key")
        require(isinstance(key, str) and key, "site profile alias key must be non-empty")
        require(key not in alias_map, f"duplicate site profile alias: {key}")
        alias_map[key] = row

    # Keep package validation self-contained. The public profile example belongs
    # to the later portable-workflows stack; this path-preflight profile is
    # frozen inside the handoff package itself.
    profile_path = package_root / "machine_profiles/bip7.path-preflight.json"
    profile = load_json(profile_path)
    profile_datasets = profile.get("datasets")
    require(isinstance(profile_datasets, dict) and profile_datasets, "site profile datasets must be a non-empty object")
    alias_joins = 0
    for object_key, row in profile_datasets.items():
        require(isinstance(row, dict), f"site profile dataset is not an object: {object_key}")
        technical_id = row.get("technical_dataset_id")
        require(technical_id in by_technical_id, f"site profile technical_dataset_id is not registered: {technical_id!r}")
        registered = by_technical_id[technical_id]
        require(row.get("biological_id") == registered.get("biological_id"), f"biological_id join mismatch: {object_key}")
        if object_key != technical_id:
            alias = alias_map.get(object_key)
            require(alias is not None, f"site profile object key has no alias registry row: {object_key}")
            require(alias.get("technical_dataset_id") == technical_id, f"alias technical_dataset_id mismatch: {object_key}")
            require(alias.get("biological_id") == row.get("biological_id"), f"alias biological_id mismatch: {object_key}")
            alias_joins += 1
    return {
        "profile": str(profile_path.relative_to(package_root)),
        "profile_datasets": len(profile_datasets),
        "technical_dataset_id_joins": len(profile_datasets),
        "object_key_alias_joins": alias_joins,
    }


class StandaloneHTMLParser(HTMLParser):
    """Small structural parser plus external runtime-resource collector."""

    VOID_TAGS = {"area", "base", "br", "col", "embed", "hr", "img", "input", "link", "meta", "param", "source", "track", "wbr"}
    RUNTIME_ATTRIBUTES = {
        "script": ("src",),
        "link": ("href",),
        "img": ("src", "srcset"),
        "iframe": ("src",),
        "embed": ("src",),
        "object": ("data",),
        "source": ("src", "srcset"),
        "audio": ("src",),
        "video": ("src", "poster"),
        "image": ("href", "xlink:href"),
        "use": ("href", "xlink:href"),
    }

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.stack: list[str] = []
        self.counts: Counter[str] = Counter()
        self.external_runtime_resources: list[str] = []
        self._style_depth = 0
        self._script_depth = 0

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        tag = tag.lower()
        self.counts[tag] += 1
        values = dict(attrs)
        for attribute in self.RUNTIME_ATTRIBUTES.get(tag, ()):
            value = values.get(attribute)
            if value and EXTERNAL_URL_RE.search(value):
                self.external_runtime_resources.append(f"{tag}[{attribute}]={value}")
        style = values.get("style")
        if style and EXTERNAL_URL_RE.search(style):
            self.external_runtime_resources.append(f"{tag}[style]")
        if tag == "style":
            self._style_depth += 1
        if tag == "script":
            self._script_depth += 1
        if tag == "meta" and values.get("http-equiv", "").lower() == "refresh":
            content = values.get("content")
            if content and EXTERNAL_URL_RE.search(content):
                self.external_runtime_resources.append("meta refresh external URL")
        if tag not in self.VOID_TAGS:
            self.stack.append(tag)

    def handle_startendtag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        self.handle_starttag(tag, attrs)
        if tag.lower() not in self.VOID_TAGS:
            self.handle_endtag(tag)

    def handle_endtag(self, tag: str) -> None:
        tag = tag.lower()
        if tag in self.VOID_TAGS:
            return
        if not self.stack or self.stack[-1] != tag:
            current = self.stack[-1] if self.stack else None
            raise ContractError(f"HTML closing tag mismatch: expected {current!r}, found {tag!r}")
        self.stack.pop()
        if tag == "style":
            self._style_depth -= 1
        if tag == "script":
            self._script_depth -= 1

    def handle_data(self, data: str) -> None:
        if self._style_depth and EXTERNAL_URL_RE.search(data):
            self.external_runtime_resources.append("style block external URL")
        if self._script_depth and EXTERNAL_URL_RE.search(data):
            self.external_runtime_resources.append("script block external URL")


def check_standalone_html(package_root: Path) -> dict[str, Any]:
    html_files = sorted(package_root.rglob("*.standalone.html"))
    require(html_files, "package contains no standalone HTML")
    external: list[str] = []
    for path in html_files:
        parser = StandaloneHTMLParser()
        try:
            parser.feed(path.read_text(encoding="utf-8"))
            parser.close()
        except (OSError, UnicodeError, ContractError) as error:
            raise ContractError(f"standalone HTML parse failed for {path.name}: {error}") from error
        require(not parser.stack, f"standalone HTML has unclosed tags in {path.name}: {parser.stack}")
        require(parser.counts["html"] == 1, f"standalone HTML must contain one html element: {path.name}")
        require(parser.counts["head"] == 1, f"standalone HTML must contain one head element: {path.name}")
        require(parser.counts["body"] == 1, f"standalone HTML must contain one body element: {path.name}")
        external.extend(f"{path.name}: {item}" for item in parser.external_runtime_resources)

        normalized = re.sub(r"\s+", " ", path.read_text(encoding="utf-8")).lower()
        require(
            "不可宣稱 production" in normalized or "非 production release" in normalized,
            f"standalone HTML lacks the no-production limitation: {path.name}",
        )
        require(
            re.search(r">\s*0\s*<.*?confirmed cellular subclone", normalized),
            f"standalone HTML lacks confirmed cellular subclone = 0: {path.name}",
        )
        require("88.2579%" in normalized, f"standalone HTML lacks the 88.2579% value: {path.name}")
        require(
            "非 accuracy/prevalence" in normalized and "biological correctness" in normalized,
            f"standalone HTML lacks the 88.2579% interpretation limits: {path.name}",
        )
    require(not external, "external runtime resources found: " + "; ".join(external))
    return {"standalone_html_files": len(html_files), "parse_errors": 0, "external_runtime_resources": 0}


def authority_expected_hashes(manifest: dict[str, Any]) -> dict[str, str]:
    expected: dict[str, str] = {}
    artifacts = manifest.get("artifacts")
    require(isinstance(artifacts, list), "authority manifest artifacts must be a list")
    for row in artifacts:
        expected[row["artifact_id"]] = row["sha256"]
    implementation = manifest.get("implementation")
    require(isinstance(implementation, dict), "authority manifest implementation must be an object")
    frozen_binary = implementation.get("frozen_binary")
    require(isinstance(frozen_binary, dict), "authority manifest frozen_binary must be an object")
    expected["frozen_binary"] = frozen_binary["sha256"]
    snapshots = implementation.get("source_snapshots")
    require(isinstance(snapshots, list), "authority manifest source_snapshots must be a list")
    for index, row in enumerate(snapshots, start=1):
        expected[f"source_snapshot_{index}"] = row["sha256_at_handoff"]
    return expected


def check_authority_replay(package_root: Path) -> dict[str, Any]:
    manifest_path = package_root / "evidence/authority_manifest.json"
    manifest = load_json(manifest_path)
    expected = authority_expected_hashes(manifest)
    require(len(expected) == 19, f"authority manifest resolves to {len(expected)} records, expected 19")
    require(all(SHA256_RE.fullmatch(value) for value in expected.values()), "authority manifest contains invalid sha256")

    receipt = load_json(package_root / "evidence/authority_replay_receipt.json")
    require(receipt.get("manifest_sha256") == sha256_file(manifest_path), "authority replay manifest hash mismatch")
    require(receipt.get("pass") is True, "authority replay receipt is not passing")
    require(receipt.get("total") == 19, "authority replay total is not 19")
    require(receipt.get("tally") == {"MATCH": 19, "MISSING": 0, "HASH_MISMATCH": 0}, "authority replay tally differs")
    results = receipt.get("results")
    require(isinstance(results, list) and len(results) == 19, "authority replay results must contain 19 records")
    by_id: dict[str, dict[str, Any]] = {}
    for row in results:
        require(isinstance(row, dict), "authority replay contains a non-object result")
        artifact_id = row.get("artifact_id")
        require(isinstance(artifact_id, str) and artifact_id, "authority replay artifact_id must be non-empty")
        require(artifact_id not in by_id, f"duplicate authority replay artifact_id: {artifact_id}")
        by_id[artifact_id] = row
    require(set(by_id) == set(expected), "authority replay identities differ from authority manifest")
    for artifact_id, expected_hash in expected.items():
        row = by_id[artifact_id]
        require(row.get("status") == "MATCH", f"authority replay status is not MATCH: {artifact_id}")
        require(row.get("expected_sha256") == expected_hash, f"authority expected sha differs: {artifact_id}")
        require(row.get("actual_sha256") == expected_hash, f"authority actual sha differs: {artifact_id}")
    return {"records": 19, "match": 19, "missing": 0, "hash_mismatch": 0}


def check_release_gate_text(package_root: Path) -> dict[str, Any]:
    index = (package_root / "00_INDEX.md").read_text(encoding="utf-8")
    required_fragments = (
        "## 目前 gate verdict",
        "| Frozen authority bytes | **PASS** | 19/19 SHA-256 MATCH |",
        "| Tag/GitHub Release | **NO-GO** |",
        "目前不得稱為 release-ready",
        "任一 blocking gate存在即不得發布",
    )
    missing = [fragment for fragment in required_fragments if fragment not in index]
    require(not missing, "release gate text is incomplete: " + repr(missing))
    return {"document": "00_INDEX.md", "required_release_gate_fragments": len(required_fragments)}


def check_reader_contract(package_root: Path) -> dict[str, Any]:
    """Guard the copy/paste and denominator issues found by a fresh reader."""

    index = (package_root / "00_INDEX.md").read_text(encoding="utf-8")
    notes = (package_root / "implementation-notes.md").read_text(encoding="utf-8")
    context = (package_root / "ai_context/CONTEXT.md").read_text(encoding="utf-8")
    html = (package_root / "20260813_完整研究交接總覽_01.standalone.html").read_text(encoding="utf-8")
    public_surfaces = "\n".join((index, notes, context, html))

    claim_pointer = load_json(package_root / "registries/claim_registry.json")
    claim_counts = claim_pointer.get("counts", {}).get("by_current_verdict", {})
    require(
        claim_counts == {"CONFIRMED": 69, "CONFIRMED_WITH_LIMITS": 69, "UNVERIFIED": 20},
        f"reader contract claim denominator drifted: {claim_counts!r}",
    )
    unverified_count = claim_counts["UNVERIFIED"]

    require(
        "34/34" not in public_surfaces,
        "P0 readiness must preserve the 33 local guards / 1 bounded About distinction",
    )
    for label, text in (("index", index), ("implementation notes", notes), ("AI context", context), ("standalone HTML", html)):
        require("33/33" in text and "C108" in text, f"{label} lacks the 33 local / 1 About P0 denominator")
        require(
            "bounded live" in text and "CONFIRMED_WITH_LIMITS" in text,
            f"{label} lacks the bounded C108 live verdict",
        )
        require(
            str(unverified_count) in text and "UNVERIFIED" in text,
            f"{label} lacks the current hash-bound UNVERIFIED denominator {unverified_count}",
        )

    require("scripts/site/bootstrap --template" in html, "standalone HTML bootstrap command lacks --template")
    require("--output config/site-profile.local.json" in html, "standalone HTML bootstrap command lacks --output")
    require("scripts/site/bootstrap --profile" not in html, "standalone HTML uses the invalid bootstrap --profile interface")
    for relative in (
        "registries/artifact_registry.json",
        "registries/claim_registry.json",
        "registries/authority_superseded_crosswalk.json",
    ):
        require(relative in context, f"AI context lacks registry prefix: {relative}")

    crosswalk = load_json(package_root / "registries/authority_superseded_crosswalk.json")
    authority = crosswalk.get("immutable_authority", {})
    authority_copy = package_path(
        package_root,
        authority.get("package_path"),
        "immutable authority package_path",
        base=package_root / "registries",
    )
    require(authority_copy.is_file(), "package-local immutable authority copy is missing")
    require(sha256_file(authority_copy) == authority.get("sha256"), "package-local immutable authority hash differs")
    source_crosswalk = next(
        (row for row in crosswalk.get("crosswalks", []) if row.get("crosswalk_id") == "frozen_binary_source_adjudication_20260813"),
        None,
    )
    require(isinstance(source_crosswalk, dict), "frozen binary source crosswalk is missing")
    for field, hash_field in (
        ("adjudication_package_path", "adjudication_sha256"),
        ("receipt_package_path", "receipt_sha256"),
    ):
        local_copy = package_path(
            package_root,
            source_crosswalk.get(field),
            field,
            base=package_root / "registries",
        )
        require(local_copy.is_file(), f"package-local crosswalk evidence is missing: {field}")
        require(sha256_file(local_copy) == source_crosswalk.get(hash_field), f"crosswalk hash differs: {field}")

    require(
        "包可獨立解讀」不等於「包含所有大型payload" in index,
        "index does not explain why frozen payloads remain in the local data plane",
    )
    return {
        "p0_local_guards": 33,
        "p0_about_bounded_live_confirmed_with_limits": 1,
        "claim_verdict_counts": claim_counts,
        "bootstrap_interface": "--template/--output",
        "registry_paths_prefixed": 3,
        "package_local_crosswalk_evidence": 3,
        "raw_authority_payload_policy": "LOCAL_DATA_PLANE",
    }


def check_reader_acceptance_receipt(package_root: Path) -> dict[str, Any]:
    receipt = package_root / "evidence/reader_acceptance_receipt.json"
    require(receipt.is_file(), "fresh-reader acceptance receipt is missing")
    validator_path = repository_root() / "scripts/handoff/validate_reader_acceptance.py"
    spec = importlib.util.spec_from_file_location("handoff_reader_acceptance_validator", validator_path)
    require(spec is not None and spec.loader is not None, "cannot load fresh-reader acceptance validator")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    errors, summary = module.validate(
        receipt,
        package_root,
        require_head=False,
        source_repo=repository_root(),
    )
    require(not errors, "fresh-reader acceptance receipt failed: " + "; ".join(errors))
    require(summary.get("verdict") == "PASS", "fresh-reader acceptance verdict is not PASS")
    return summary


CHECKS: tuple[tuple[str, Callable[[Path], dict[str, Any]]], ...] = (
    ("markdown_links", check_markdown_links),
    ("json_parse", check_json_parse),
    ("evidence_manifest", check_evidence_manifest),
    ("dataset_registry", check_dataset_registry),
    ("run_registry", check_run_registry),
    ("artifact_registry", check_artifact_registry),
    ("claim_registries", check_claim_registries),
    ("site_profile_join", check_site_profile_join),
    ("standalone_html", check_standalone_html),
    ("authority_replay", check_authority_replay),
    ("release_gate_text", check_release_gate_text),
    ("reader_contract", check_reader_contract),
    ("reader_acceptance_receipt", check_reader_acceptance_receipt),
)


def validate_package(package_root: Path) -> dict[str, Any]:
    package_root = package_root.resolve()
    checks: list[dict[str, Any]] = []
    if not package_root.is_dir():
        checks.append({"check_id": "package_root", "status": "FAIL", "error": f"not a directory: {package_root}"})
    else:
        checks.append({"check_id": "package_root", "status": "PASS", "details": {"path": str(package_root)}})
        for check_id, function in CHECKS:
            try:
                details = function(package_root)
                checks.append({"check_id": check_id, "status": "PASS", "details": details})
            except Exception as error:  # Fail closed and preserve all independent gate results.
                checks.append({"check_id": check_id, "status": "FAIL", "error": f"{type(error).__name__}: {error}"})
    passed = all(row["status"] == "PASS" for row in checks)
    return {
        "schema_version": "1.0.0",
        "receipt_type": "handoff_package_validation",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "package_root": str(package_root),
        "pass": passed,
        "tally": {
            "PASS": sum(row["status"] == "PASS" for row in checks),
            "FAIL": sum(row["status"] == "FAIL" for row in checks),
        },
        "checks": checks,
    }


def default_package_root() -> Path:
    return Path(__file__).resolve().parents[2] / PACKAGE_RELATIVE_PATH


def write_receipt(path: Path, receipt: dict[str, Any]) -> None:
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(receipt, ensure_ascii=False, indent=2) + "\n"
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_text(payload, encoding="utf-8")
    os.replace(temporary, path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package_root", nargs="?", type=Path, default=default_package_root())
    parser.add_argument("--receipt", type=Path, help="write the JSON receipt to this explicit path")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    receipt = validate_package(arguments.package_root)
    if arguments.receipt is not None:
        write_receipt(arguments.receipt, receipt)
    json.dump(receipt, sys.stdout, ensure_ascii=False, indent=2)
    sys.stdout.write("\n")
    return 0 if receipt["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
