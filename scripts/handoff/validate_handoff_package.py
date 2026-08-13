#!/usr/bin/env python3
"""Fail-closed validation for the 2026-08-13 research handoff package."""

from __future__ import annotations

import argparse
import hashlib
import json
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
    for row in records:
        evidence_id = row["evidence_id"]
        path = package_path(package_root, row.get("path"), f"evidence {evidence_id} path")
        require(path.is_file(), f"evidence path does not exist: {row.get('path')}")
        copy_status = row.get("copy_status")
        if copy_status in {"EXACT_COPY", "NEW_PACKAGE_RECEIPT", "PACKAGE_RELINKED_COPY"}:
            expected = row.get("sha256")
            require(isinstance(expected, str) and SHA256_RE.fullmatch(expected), f"invalid sha256 for {evidence_id}")
            require(sha256_file(path) == expected, f"evidence sha256 mismatch: {evidence_id}")
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
    return {
        "records": len(records),
        "sha256_bound_records": sha256_bound_count,
        "summary_hash_bound": summary_count,
        "missing": 0,
        "hash_mismatch": 0,
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

    tagged_bams = [row for row in records if row.get("artifact_type") == "longphase_s_tagged_bam"]
    require(len(tagged_bams) == 14, f"expected 14 tagged BAM artifacts, found {len(tagged_bams)}")
    for row in tagged_bams:
        require(row.get("scope") == "PARTIAL", f"tagged BAM is not PARTIAL: {row.get('artifact_id')}")
        require(row.get("evidence_status") == "PARTIAL", f"tagged BAM evidence is not PARTIAL: {row.get('artifact_id')}")
        require(row.get("finality") == "NON_FINAL", f"tagged BAM is not NON_FINAL: {row.get('artifact_id')}")
    return {
        "artifacts": 36,
        "unique_artifact_ids": 36,
        "final_for_scope": 20,
        "science_authority_final": 19,
        "provenance_only_validated_derived_final": 1,
        "tagged_bam_partial_non_final": 14,
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
    return {"claims": 158, "unique_claim_ids": 158, "sha256_match": True, "release_gates": expected_gates}


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

    require("34/34" not in public_surfaces, "P0 readiness incorrectly promotes the external-only C108 claim")
    require("33/33" in index and "C108" in index, "index lacks the 33 local / 1 external P0 denominator")
    require("33/33" in notes and "C108" in notes, "implementation notes lack the P0 denominator boundary")
    require("33/33" in html and "C108" in html, "standalone HTML lacks the P0 denominator boundary")

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
        "p0_external_only_blocked": 1,
        "bootstrap_interface": "--template/--output",
        "registry_paths_prefixed": 3,
        "package_local_crosswalk_evidence": 3,
        "raw_authority_payload_policy": "LOCAL_DATA_PLANE",
    }


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
    path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


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
