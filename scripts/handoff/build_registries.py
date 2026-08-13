#!/usr/bin/env python3
"""Build the 2026-08-13 handoff registries from frozen and live inventories.

The builder intentionally separates logical manifest metadata from physical
directories.  It performs metadata-only scans of large storage roots and never
hashes BAM payloads or recursively sizes TiB-scale trees.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import tempfile
import time
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


# A live inventory must not inherit an earlier run's verification timestamp.
# Set INTERSUBMOD_HANDOFF_VERIFIED_AT only when replaying a declared snapshot.
VERIFY_TIME = os.environ.get(
    "INTERSUBMOD_HANDOFF_VERIFIED_AT",
    datetime.now().astimezone().isoformat(timespec="seconds"),
)
AUTHORITY_CREATED = "2026-08-01T00:00:00+08:00"
GENOME_BUILD = "GRCh38"
EXPECTED_TAGGED_BAM_BYTES = 3_709_322_840_333
TECHNICAL_TO_BIOLOGICAL = {
    "HCC1395": "HCC1395",
    "HCC1395_DORADO": "HCC1395",
    "COLO829": "COLO829",
    "H1437": "H1437",
    "H2009": "H2009",
    "HCC1937": "HCC1937",
    "HCC1954": "HCC1954",
}
EVIDENCE_ENUM = {"AUTHORITY", "VALIDATED_DERIVED", "PARTIAL", "HISTORICAL", "INVALIDATED", "IN_PROGRESS"}
SCOPE_ENUM = {"FULL", "PARTIAL", "DEMO"}
FINALITY_ENUM = {"FINAL_FOR_SCOPE", "NON_FINAL", "SUPERSEDED"}
AVAILABILITY_ENUM = {"GIT", "GITHUB_RELEASE", "LOCAL_CANONICAL", "EXTERNAL_SOURCE", "MISSING"}
REGENERATION_SEMANTIC_PREFIXES = ("REPLAY_ONLY:", "VERIFY_ONLY:", "NOT_REGENERABLE_FROM_HANDOFF:")
AUTHORITY_REPLAY_COMMAND = (
    "REPLAY_ONLY: python3 scripts/handoff/replay_authority.py "
    "--manifest docs/handoff/20260813_完整研究資料與軟體交接_01/evidence/authority_manifest.json "
    "--output-json <NEW_RECEIPT_PATH>; verifies frozen bytes only and does not rerun science"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def iso_from_stat(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(timespec="seconds")


def parse_datetime(value: Any) -> str | None:
    if not value:
        return None
    text = str(value).strip()
    normalized = text.replace(" CST", "+08:00")
    if re.fullmatch(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", normalized):
        normalized = normalized.replace(" ", "T") + "+08:00"
    elif " " in normalized and "T" not in normalized:
        normalized = normalized.replace(" ", "T", 1)
    try:
        parsed = datetime.fromisoformat(normalized)
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.astimezone()
    return parsed.isoformat(timespec="seconds")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_datasets(manifest: Mapping[str, Any], run_rows: Sequence[Mapping[str, str]]) -> list[dict[str, Any]]:
    platforms: dict[str, str] = {}
    for row in run_rows:
        platforms.setdefault(row["sample"], row["platform"])
    datasets = []
    for technical_id in manifest["scope"]["technical_datasets"]:
        datasets.append(
            {
                "technical_dataset_id": technical_id,
                "biological_id": TECHNICAL_TO_BIOLOGICAL[technical_id],
                "platform": platforms.get(technical_id, "UNKNOWN"),
                "genome_build": GENOME_BUILD,
                "scope": "chr1-22",
                "authority_status": "AUTHORITY_SCOPE_MEMBER",
                "known_limits": [
                    "Technical datasets are not independent biological replicates.",
                    (
                        "HCC1395 and HCC1395_DORADO are one biological cell line measured by different technical pipelines."
                        if technical_id in {"HCC1395", "HCC1395_DORADO"}
                        else "No independent cellular lineage truth is registered."
                    ),
                ],
            }
        )
    return datasets


def classify_run(path: Path) -> str:
    name = path.name.lower()
    if (path / "manifest" / "run_context.json").is_file():
        return "MATERIALIZED_MANIFEST_BUNDLE"
    if "complete_matrix" in name:
        return "COMPLETE_MATRIX"
    if "smoke" in name:
        return "SMOKECHECK"
    if "permanova" in name:
        return "PERMANOVA_RERUN"
    if "parallel_test" in name:
        return "PARALLEL_TEST"
    if re.match(r"202604(20|21)", name):
        return "APR_FULL_RERUN"
    return "OTHER"


def infer_mode(path: Path) -> str:
    for mode in ("paired_full", "paired_pileup", "to_pileup"):
        if mode in path.parts:
            return mode
    return "unknown"


def physical_run_dirs(root: Path) -> list[Path]:
    return sorted(path for path in root.glob("*/*/*") if path.is_dir())


def run_context(path: Path) -> tuple[dict[str, Any], str | None]:
    for candidate in (path / "manifest" / "run_context.json", path / "run_context.json", path / "round_context.json"):
        if candidate.is_file():
            try:
                return json.loads(candidate.read_text(encoding="utf-8")), str(candidate)
            except json.JSONDecodeError:
                return {}, str(candidate)
    return {}, None


def pending_path_for_manifest(canonical_path: Path, canonical_root: Path, archive_root: Path) -> Path:
    return archive_root / canonical_path.relative_to(canonical_root)


def build_runs(
    run_rows: Sequence[Mapping[str, str]], canonical_root: Path, archive_root: Path
) -> dict[str, Any]:
    current = physical_run_dirs(canonical_root)
    pending = physical_run_dirs(archive_root)
    all_physical = [(path, "CURRENT") for path in current] + [(path, "PENDING_ARCHIVE") for path in pending]
    physical_by_path = {str(path): (path, storage_state) for path, storage_state in all_physical}

    logical_by_physical: dict[str, Mapping[str, str]] = {}
    for row in run_rows:
        canonical_path = Path(row["canonical_run_path"])
        pending_path = pending_path_for_manifest(canonical_path, canonical_root, archive_root)
        if canonical_path.is_dir():
            selected = canonical_path
        elif pending_path.is_dir():
            selected = pending_path
        else:
            raise ValueError(f"logical run lacks a matching physical directory: {row['run_id']}")
        if str(selected) in logical_by_physical:
            raise ValueError(f"multiple logical rows map to one physical directory: {selected}")
        logical_by_physical[str(selected)] = row

    records: list[dict[str, Any]] = []
    for path, storage_state in all_physical:
        logical = logical_by_physical.get(str(path))
        context, context_path = run_context(path)
        layout = classify_run(path)
        scope = "DEMO" if layout in {"SMOKECHECK", "PARALLEL_TEST"} else (
            "PARTIAL" if layout == "PERMANOVA_RERUN" else "FULL"
        )
        evidence = "HISTORICAL" if layout in {"SMOKECHECK", "PARALLEL_TEST", "PERMANOVA_RERUN"} else "PARTIAL"
        sample_id = str(logical["sample"]) if logical else str(context.get("sample") or path.parts[-3])
        mode = str(logical["canonical_mode"]) if logical else infer_mode(path)
        source_created = parse_datetime(context.get("created_at") or context.get("materialized_at"))
        known_limits = ["Physical directory presence and name do not establish scientific finality."]
        if logical:
            known_limits.extend(
                [
                    f"master manifest completeness_state={logical['completeness_state']}",
                    f"master manifest blocking_reason={logical['blocking_reason']}",
                    "Logical master-manifest metadata was merged into this physical row; it is not a second run.",
                    "tagged_bam_ready=false is a stale logical field and must not override physical artifact inventory.",
                ]
            )
        else:
            known_limits.append("Not present in the 18-row logical master manifest; directory is explicitly unregistered.")
        if storage_state == "PENDING_ARCHIVE":
            known_limits.append("Located under _ARCHIVE_pending_cleanup_202606; archival receipt/finality is not implied.")

        records.append(
            {
                "physical_run_id": f"{storage_state.lower()}::{path.relative_to(canonical_root if storage_state == 'CURRENT' else archive_root)}",
                "run_id": path.name,
                "sample_id": sample_id,
                "mode": mode,
                "layout_class": layout,
                "storage_state": storage_state,
                "physical_path": str(path),
                "logical_manifest_member": logical is not None,
                "logical_manifest_run_id": str(logical["run_id"]) if logical else None,
                "logical_manifest_metadata": dict(logical) if logical else None,
                "run_context_path": context_path,
                "source_created_at": source_created,
                "run_started_at": parse_datetime(context.get("run_started_at") or context.get("started_at")),
                "run_completed_at": parse_datetime(context.get("run_completed_at") or context.get("completed_at")),
                "verified_at": VERIFY_TIME,
                "published_at": None,
                "last_seen_at": VERIFY_TIME,
                "filesystem_mtime": iso_from_stat(path),
                "scope": scope,
                "evidence_status": evidence,
                "finality": "NON_FINAL",
                "availability": "LOCAL_CANONICAL",
                "machine_locations": [
                    {
                        "machine": "bip7",
                        "path": str(path),
                        "last_seen_at": VERIFY_TIME,
                        "location_status": "ARCHIVE_PENDING_UNRECEIPTED" if storage_state == "PENDING_ARCHIVE" else "PRESENT",
                    }
                ],
                "known_limits": known_limits,
            }
        )

    state_logical = Counter(
        row["storage_state"] for row in records if row["logical_manifest_member"]
    )
    state_extra = Counter(
        row["storage_state"] for row in records if not row["logical_manifest_member"]
    )
    reconciliation = {
        "manifest_file_lines": len(run_rows) + 1,
        "manifest_header_lines": 1,
        "logical_manifest_rows": len(run_rows),
        "physical_directories_total": len(records),
        "current_physical_directories": len(current),
        "pending_archive_physical_directories": len(pending),
        "logical_rows_merged_current": state_logical["CURRENT"],
        "logical_rows_merged_pending_archive": state_logical["PENDING_ARCHIVE"],
        "current_physical_unregistered": state_extra["CURRENT"],
        "pending_archive_extras": state_extra["PENDING_ARCHIVE"],
        "explanation_18_19_35": (
            "18 logical data rows + 1 TSV header = 19 file lines; 35 is the current canonical physical-directory count. "
            "The complete physical inventory is 35 current + 16 pending archive = 51."
        ),
        "merge_policy": (
            "A logical row is metadata merged into its matching physical directory by canonical run path or the same "
            "sample/mode/run relative path under _ARCHIVE_pending_cleanup_202606/canonical; it never creates a duplicate record."
        ),
    }
    return {
        "schema_version": "2.0.0",
        "registry_type": "physical_run_registry_with_logical_metadata_merge",
        "verified_at": VERIFY_TIME,
        "claim_ceiling": "Filesystem/run provenance inventory only; no directory name or manifest flag establishes scientific finality.",
        "reconciliation": reconciliation,
        "records": sorted(records, key=lambda row: (row["storage_state"], row["sample_id"], row["mode"], row["run_id"])),
    }


def machine_location(path: Path, status: str = "PRESENT", machine: str = "bip7") -> list[dict[str, Any]]:
    return [{"machine": machine, "path": str(path), "last_seen_at": VERIFY_TIME, "location_status": status}]


def artifact_row(**kwargs: Any) -> dict[str, Any]:
    fields = {
        "artifact_id": None,
        "artifact_type": None,
        "semantic_description": None,
        "sample_id": None,
        "technical_dataset_id": None,
        "mode": None,
        "genome_build": None,
        "producer": None,
        "producer_commit": None,
        "schema_version": "unknown",
        "inputs": [],
        "derived_from": [],
        "validates": [],
        "supersedes": [],
        "used_by": [],
        "created_at": None,
        "generated_at": None,
        "filesystem_mtime": None,
        "verified_at": VERIFY_TIME,
        "scope": "PARTIAL",
        "evidence_status": "PARTIAL",
        "finality": "NON_FINAL",
        "availability": "MISSING",
        "size_bytes": None,
        "sha256": None,
        "license": "UNKNOWN_REVIEW_REQUIRED",
        "claim_ceiling": None,
        "known_limits": [],
        "regeneration_command": None,
        "public_uri": None,
        "machine_locations": [],
    }
    fields.update(kwargs)
    return fields


def tagged_bam_paths(canonical_root: Path) -> list[Path]:
    return sorted(canonical_root.glob("*/*/*/longphase_s/*_tagged.bam"))


def build_artifacts(
    manifest: Mapping[str, Any], authority_path: Path, canonical_root: Path
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    repository_root = authority_path.resolve().parents[3]
    for item in manifest["artifacts"]:
        path = Path(item["path"])
        rows.append(
            artifact_row(
                artifact_id=f"authority.{item['artifact_id']}",
                artifact_type=item["role"],
                semantic_description=f"Frozen 2026-08-01 authority artifact: {item['artifact_id']}",
                genome_build=GENOME_BUILD,
                producer="frozen exact-PS / methylation research pipeline",
                producer_commit=manifest["implementation"]["inter_sub_mod_git_head_at_handoff"],
                schema_version=str(item["schema_version"]),
                inputs=["authority_manifest_20260801"],
                validates=[item["scope"]],
                created_at=iso_from_stat(path) if path.exists() else None,
                generated_at=iso_from_stat(path) if path.exists() else None,
                filesystem_mtime=iso_from_stat(path) if path.exists() else None,
                scope="FULL",
                evidence_status="AUTHORITY",
                finality="FINAL_FOR_SCOPE",
                availability="GIT" if path.exists() and path.resolve().is_relative_to(repository_root) else "LOCAL_CANONICAL",
                size_bytes=path.stat().st_size if path.exists() else None,
                sha256=item["sha256"],
                license="GPL-3.0-or-later for repository content; upstream dataset license tracked separately",
                claim_ceiling=item["claim_boundary"],
                known_limits=list(manifest["known_limitations"]),
                regeneration_command=AUTHORITY_REPLAY_COMMAND,
                machine_locations=machine_location(path, "PRESENT" if path.exists() else "MISSING"),
            )
        )

    implementation = manifest["implementation"]
    binary = implementation["frozen_binary"]
    binary_path = Path(binary["path"])
    rows.append(
        artifact_row(
            artifact_id="authority.frozen_binary",
            artifact_type="executable_binary",
            semantic_description="Frozen exact_ps_topology_af executable used by the authority bundle",
            producer="InterSubMod research CLI linked with LongLineage solver modules",
            producer_commit="InterSubMod 387a101e6a3292e0d7f230ba8d20271c7434972a + LongLineage b979760 (append-only adjudication)",
            schema_version="binary-build-tree-v1",
            inputs=["authority.source_snapshot_1", "authority.source_snapshot_5", "adjudication.frozen_binary_source"],
            derived_from=["authority_manifest_20260801", "adjudication.frozen_binary_source"],
            used_by=["authority.topology_summary", "authority.all7_summary"],
            created_at=iso_from_stat(binary_path) if binary_path.exists() else None,
            generated_at=iso_from_stat(binary_path) if binary_path.exists() else None,
            filesystem_mtime=iso_from_stat(binary_path) if binary_path.exists() else None,
            scope="FULL",
            evidence_status="AUTHORITY",
            finality="FINAL_FOR_SCOPE",
            availability="LOCAL_CANONICAL",
            size_bytes=binary_path.stat().st_size if binary_path.exists() else None,
            sha256=binary["sha256"],
            license="GPL-3.0-or-later; linked source-origin recorded in append-only adjudication",
            claim_ceiling="Frozen binary byte identity and frozen-result reproduction only; not a production LongLineage release.",
            known_limits=[
                "2026-08-01's not_proven status is preserved in the immutable manifest.",
                "2026-08-13 append-only adjudication proves 5/5 solver modules byte-identical to LongLineage b979760.",
                "source_snapshot_4 is a later LongLineage worktree byte snapshot and is explicitly rejected as binary source truth.",
            ],
            regeneration_command=AUTHORITY_REPLAY_COMMAND,
            machine_locations=machine_location(binary_path, "PRESENT" if binary_path.exists() else "MISSING"),
        )
    )
    source_roles = {
        1: ("InterSubMod research CLI compiled into the frozen binary", ["authority.frozen_binary"]),
        2: ("All-seven execution wrapper; not a compiled binary input", ["authority.all7_cohort_receipt"]),
        3: ("Partition preparation script; not a compiled binary input", []),
        4: ("Later LongLineage worktree snapshot rejected as frozen-binary source truth", []),
        5: ("LongLineage parent_mapping module compiled into the frozen binary", ["authority.frozen_binary"]),
    }
    source_content_commits = {
        1: "UNCOMMITTED_AT_20260801; exact bytes later committed as 886fb5e648978f80b5d435360ce1c16432fa6843",
        2: "UNCOMMITTED_AT_20260801; exact bytes later committed as 6583ae1382befd20ff4a4c4ff076cb17020537de",
        3: "UNCOMMITTED_AT_20260801; exact bytes later committed as 23aef3f5527bd2dd9b4b05712414b68f8183baea",
        4: "UNCOMMITTED_AT_20260801; exact bytes later committed as 9ad976b10eb51589a043b55ccc52b6232df9fa3a",
        5: "b9797605cdc859e571c1961c1b4343c40d7c7183",
    }
    for index, source in enumerate(implementation["source_snapshots"], start=1):
        path = Path(source["path"])
        rejected = index == 4
        semantic_role, consumers = source_roles[index]
        rows.append(
            artifact_row(
                artifact_id=f"authority.source_snapshot_{index}",
                artifact_type="registered_source_byte_snapshot",
                semantic_description=f"Immutable 2026-08-01 registered source snapshot {index}: {semantic_role}",
                producer="2026-08-01 handoff inventory",
                producer_commit=source_content_commits[index],
                schema_version="source-snapshot-v1",
                inputs=["authority_manifest_20260801", f"registered filesystem source object: {path}"],
                derived_from=["authority_manifest_20260801"],
                used_by=consumers,
                filesystem_mtime=iso_from_stat(path) if path.exists() else None,
                scope="FULL",
                evidence_status="AUTHORITY",
                finality="FINAL_FOR_SCOPE",
                availability="LOCAL_CANONICAL",
                size_bytes=path.stat().st_size if path.exists() else None,
                sha256=source["sha256_at_handoff"],
                license="GPL-3.0-or-later; cross-repo source-origin recorded separately",
                claim_ceiling=(
                    "Immutable record of bytes observed on 2026-08-01; REJECTED as frozen-binary source truth by adjudication."
                    if rejected
                    else "Immutable registered source bytes; binary-source role is governed by the append-only adjudication."
                ),
                known_limits=(
                    [
                        "The bytes remain part of the immutable 2026-08-01 19-item replay.",
                        "This snapshot includes a 2026-07-26 seed_incumbent change and post-dates the 2026-07-24 binary.",
                        "Do not list this artifact as an input to authority.frozen_binary.",
                    ]
                    if rejected
                    else ["Current mtimes and Git HEAD are not substitutes for the frozen SHA-256."]
                ),
                regeneration_command=(
                    "VERIFY_ONLY: use python3 scripts/handoff/replay_authority.py with the package authority manifest; "
                    "this immutable source snapshot is not regenerated"
                ),
                machine_locations=machine_location(path, "PRESENT" if path.exists() else "MISSING"),
            )
        )

    for path in tagged_bam_paths(canonical_root):
        stat_result = path.stat()
        sample = path.name.removesuffix("_tagged.bam")
        mode = infer_mode(path)
        rows.append(
            artifact_row(
                artifact_id=f"canonical.pre_fix_tagged_bam.{sample.lower()}.{mode}",
                artifact_type="longphase_s_tagged_bam",
                semantic_description=f"Canonical paired {mode} LongPhase-S tagged BAM physical payload (PRE-FIX)",
                sample_id=sample,
                technical_dataset_id=sample,
                mode=mode,
                genome_build=GENOME_BUILD,
                producer="LongPhase-S v1.7.3 --somaticMode historical canonical workflow",
                producer_commit=None,
                schema_version="BAM-with-HP-PS-MM-ML",
                inputs=["external tumor/normal BAM", "historical ClairS candidate VCF/truth-BED-scoped workflow"],
                used_by=["historical ASM and LOH/HP analyses"],
                created_at=iso_from_stat(path),
                generated_at=iso_from_stat(path),
                filesystem_mtime=iso_from_stat(path),
                scope="PARTIAL",
                evidence_status="PARTIAL",
                finality="NON_FINAL",
                availability="LOCAL_CANONICAL",
                size_bytes=stat_result.st_size,
                sha256=None,
                license="Upstream public research dataset license and redistribution review required; payload remains local",
                claim_ceiling="PRE-FIX read-level research input only; never production authority or release evidence.",
                known_limits=[
                    "SHA-256 intentionally NOT_COMPUTED_TIB_PAYLOAD; exact path/bytes/mtime are inventory metadata only.",
                    "Canonical logical manifest tagged_bam_ready=false is stale and conflicts with physical presence.",
                    "6/7 datasets were produced with truth-BED-restricted somatic evidence; cross-sample production use is blocked.",
                    "No BAM payload is included in Git or the handoff package.",
                ],
                regeneration_command=None,
                machine_locations=machine_location(path),
            )
        )

    rows.extend(
        [
            artifact_row(
                artifact_id="adjudication.frozen_binary_source",
                artifact_type="append_only_source_adjudication",
                semantic_description="2026-08-13 append-only correction linking the frozen binary to InterSubMod CLI and LongLineage b979760",
                producer="frozen build-tree cmp/objdump/nm adjudication",
                producer_commit="e3c0889100da88e6056e9268e07f6f94f8e9312b",
                schema_version="1.0.0",
                inputs=["authority_manifest_20260801", "frozen build tree", "LongLineage b979760"],
                derived_from=["authority.source_snapshot_1", "authority.source_snapshot_4", "authority.frozen_binary"],
                validates=["5/5 LongLineage solver modules byte-identical to frozen build objects"],
                supersedes=["status:not_proven_byte_reproducible", "source-role:authority.source_snapshot_4"],
                used_by=["authority.frozen_binary"],
                created_at="2026-08-13T10:49:23+08:00",
                generated_at="2026-08-13T10:49:23+08:00",
                filesystem_mtime=iso_from_stat(
                    repository_root / "docs/provenance/20260813_frozen_binary_source_adjudication_01.md"
                ),
                scope="FULL",
                evidence_status="VALIDATED_DERIVED",
                finality="FINAL_FOR_SCOPE",
                availability="GIT",
                size_bytes=(repository_root / "docs/provenance/20260813_frozen_binary_source_adjudication_01.md").stat().st_size,
                sha256=sha256_file(repository_root / "docs/provenance/20260813_frozen_binary_source_adjudication_01.md"),
                license="GPL-3.0-or-later",
                claim_ceiling="Corrects source identity only; does not alter frozen scientific outputs or make LongLineage production-ready.",
                known_limits=["Append-only correction; the 2026-08-01 authority manifest remains immutable."],
                regeneration_command=(
                    "VERIFY_ONLY: validate docs/provenance/20260813_frozen_binary_bundle_receipt.json and "
                    "docs/provenance/20260813_frozen_binary_source_adjudication_01.md; no science rerun"
                ),
                machine_locations=machine_location(repository_root / "docs/provenance/20260813_frozen_binary_source_adjudication_01.md"),
            ),
            artifact_row(
                artifact_id="legacy.big7_runbook_original",
                artifact_type="legacy_runbook",
                semantic_description="Original InterSubMod_big7_runbook directory documented through 2026-06-01",
                producer="legacy research operations",
                scope="PARTIAL",
                evidence_status="HISTORICAL",
                finality="NON_FINAL",
                availability="MISSING",
                license="UNKNOWN_REVIEW_REQUIRED",
                claim_ceiling="Existence documented through 2026-06-01 only.",
                known_limits=["Original bytes, hashes, deletion or move receipt are unavailable after broad bounded search."],
            ),
            artifact_row(
                artifact_id="claim.invalidated_canonical_19_runs",
                artifact_type="invalidated_claim",
                semantic_description="Legacy claim that master_run_manifest.tsv contains 19 run rows",
                producer="legacy line-count calculation",
                schema_version="claim-v1",
                scope="FULL",
                evidence_status="INVALIDATED",
                finality="SUPERSEDED",
                availability="GIT",
                license="GPL-3.0-or-later",
                claim_ceiling="Invalidated: 19 file lines equals 1 header plus 18 logical rows; physical directories are tracked separately.",
                known_limits=["Do not invent a nineteenth logical run or collapse 51 physical directories into 18 logical records."],
            ),
        ]
    )
    return sorted(rows, key=lambda row: row["artifact_id"])


def mount_info(path: Path) -> dict[str, Any] | None:
    result = subprocess.run(
        ["findmnt", "-T", str(path), "-J", "-o", "TARGET,SOURCE,FSTYPE,OPTIONS"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        timeout=10,
    )
    if result.returncode != 0:
        return None
    filesystems = json.loads(result.stdout).get("filesystems", [])
    return filesystems[0] if filesystems else None


def git_state(path: Path) -> dict[str, Any] | None:
    if not (path / ".git").exists() and not (path / ".git").is_file():
        return None
    try:
        head = subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True, timeout=10).strip()
        branch = subprocess.check_output(["git", "-C", str(path), "branch", "--show-current"], text=True, timeout=10).strip()
        raw = subprocess.check_output(
            ["git", "-C", str(path), "status", "--porcelain=v1", "-z", "--untracked-files=all"],
            timeout=30,
        )
        return {"head": head, "branch": branch, "dirty_entry_count": len([item for item in raw.split(b"\0") if item]), "clean": not raw}
    except (subprocess.SubprocessError, OSError):
        return {"head": None, "branch": None, "dirty_entry_count": None, "clean": None}


def machine_path_record(
    path_id: str,
    path: Path,
    semantic_role: str,
    *,
    expected_machine: str,
    observation_status: str | None = None,
    known_limits: Sequence[str] = (),
) -> dict[str, Any]:
    exists = path.exists()
    actual_host = platform.node()
    mount = mount_info(path if exists else path.parent)
    status = observation_status or ("PRESENT" if exists else "MISSING")
    if exists and expected_machine != actual_host and mount and str(mount.get("fstype", "")).startswith("nfs"):
        status = "NFS_VISIBLE_NOT_HOST_VERIFIED"
    return {
        "path_id": path_id,
        "semantic_role": semantic_role,
        "path": str(path),
        "expected_machine": expected_machine,
        "observed_from_host": actual_host,
        "last_seen_at": VERIFY_TIME,
        "observation_status": status,
        "entry_type": "DIRECTORY" if path.is_dir() else "FILE" if path.is_file() else "MISSING",
        "size_bytes": path.stat().st_size if path.is_file() else None,
        "mtime": iso_from_stat(path) if exists else None,
        "mount": mount,
        "git_state": git_state(path) if path.is_dir() else None,
        "known_limits": list(known_limits),
    }


def build_machine_paths(authority_path: Path) -> dict[str, Any]:
    repo = authority_path.resolve().parents[3]
    specs = [
        ("repo.intersubmod_dirty", Path("/big7_disk/liaoyoyo2001/InterSubMod"), "Active/dirty InterSubMod research repository", "bip7"),
        ("repo.intersubmod_release_worktree", repo, "Clean-baseline-derived handoff worktree", "bip7"),
        ("repo.longlineage_active", Path("/big7_disk/liaoyoyo2001/LongLineage"), "Active LongLineage repository", "bip7"),
        ("repo.longlineage_preview", Path("/big7_disk/liaoyoyo2001/worktrees/ll-preview-20260813"), "LongLineage public-preview worktree", "bip7"),
        ("storage.big7_output", Path("/big7_disk/liaoyoyo2001/big7_disk_output"), "Primary big7 output root", "bip7"),
        ("storage.big7_canonical", Path("/big7_disk/liaoyoyo2001/big7_disk_output/canonical"), "Current canonical physical runs", "bip7"),
        ("storage.big7_synthesis", Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis"), "Research synthesis and observations", "bip7"),
        ("storage.pending_archive", Path("/big7_disk/liaoyoyo2001/big7_disk_output/_ARCHIVE_pending_cleanup_202606"), "Pending-cleanup archive root", "bip7"),
        ("storage.big7_big8_archive", Path("/big7_disk/liaoyoyo2001/big7_disk_output/big8_output_archive"), "Materialized big8 archive", "bip7"),
        ("storage.big7_bip8_archive", Path("/big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive"), "Materialized bip8 archive", "bip7"),
        ("legacy.big8_output", Path("/big8_disk/liaoyoyo2001/InterSubMod_runs/output"), "Legacy big8 InterSubMod output root", "big8"),
        ("legacy.bip8_output", Path("/bip8_disk/liaoyoyo2001/InterSubMod_out/output"), "Legacy bip8 InterSubMod output root", "bip8"),
        ("meeting.root", Path("/big7_disk/liaoyoyo2001/Meeting"), "Meeting notes and presentation workspace", "bip7"),
        ("legacy.runbook", Path("/big7_disk/liaoyoyo2001/InterSubMod_big7_runbook"), "Historical runbook root", "bip7"),
        ("binary.frozen_original", Path("/bip7_disk/liaoyoyo2001/build-exact-af-20260724/bin/exact_ps_topology_af"), "Original frozen exact-PS binary", "bip7"),
        ("binary.frozen_bundle", Path("/bip7_disk/liaoyoyo2001/_frozen/exact_ps_topology_af_20260724"), "Permanent frozen binary/source/build bundle", "bip7"),
    ]
    records = []
    for path_id, path, role, machine in specs:
        limits = []
        if machine == "bip8":
            limits.append("Observed through NFS from bip7; this is not an independent bip8 host receipt.")
        if path_id == "legacy.runbook":
            limits.append("MISSING_SOURCE; reconstructed indexes must never impersonate the original directory.")
        records.append(machine_path_record(path_id, path, role, expected_machine=machine, known_limits=limits))
    return {
        "schema_version": "1.0.0",
        "registry_type": "machine_path_registry",
        "verified_at": VERIFY_TIME,
        "actual_host": platform.node(),
        "claim_ceiling": "Path and mount observations only; NFS visibility is not host-local bip8 validation.",
        "records": records,
    }


def bounded_tree_size(root: Path, *, max_entries: int = 10_000, timeout_seconds: float = 5.0) -> tuple[int | None, str, str | None]:
    """Measure a known-small tree, returning an explicit blocker on bounds."""

    started = time.monotonic()
    entries = 0
    total = 0
    try:
        for current, directories, files in os.walk(root, followlinks=False):
            entries += len(directories) + len(files)
            if entries > max_entries:
                return None, "BLOCKED_SAFETY_BOUND", f"entry count exceeded safe bound {max_entries}"
            if time.monotonic() - started > timeout_seconds:
                return None, "BLOCKED_TIMEOUT", f"metadata walk exceeded {timeout_seconds:.1f}s timeout"
            current_path = Path(current)
            for name in files:
                path = current_path / name
                total += path.lstat().st_size
        return total, "MEASURED", None
    except OSError as error:
        return None, "BLOCKED_IO_ERROR", f"metadata walk failed: {error}"


def build_storage_manifest(machine_registry: Mapping[str, Any]) -> dict[str, Any]:
    storage_ids = {
        "storage.big7_output",
        "storage.big7_canonical",
        "storage.big7_synthesis",
        "storage.pending_archive",
        "storage.big7_big8_archive",
        "storage.big7_bip8_archive",
        "legacy.big8_output",
        "legacy.bip8_output",
        "meeting.root",
        "legacy.runbook",
        "binary.frozen_bundle",
    }
    records = []
    for machine_row in machine_registry["records"]:
        if machine_row["path_id"] not in storage_ids:
            continue
        path = Path(machine_row["path"])
        children: list[dict[str, Any]] = []
        if path.is_dir():
            for child in sorted(path.iterdir(), key=lambda item: item.name):
                children.append(
                    {
                        "name": child.name,
                        "path": str(child),
                        "entry_type": "SYMLINK" if child.is_symlink() else "DIRECTORY" if child.is_dir() else "FILE",
                        "size_bytes": child.lstat().st_size if child.is_file() or child.is_symlink() else None,
                        "mtime": iso_from_stat(child),
                    }
                )
        recursive_size: int | None = None
        size_status = "NOT_APPLICABLE"
        size_blocker: str | None = None
        if path.is_dir() and machine_row["path_id"] == "binary.frozen_bundle":
            recursive_size, size_status, size_blocker = bounded_tree_size(path)
        elif path.is_dir():
            size_status = "NOT_ATTEMPTED_TIB_SCALE_SAFETY"
            size_blocker = (
                "Recursive size intentionally not measured: known TiB payloads or unbounded research trees; "
                "use targeted artifact stat."
            )
        records.append(
            {
                "root_id": machine_row["path_id"],
                "path": machine_row["path"],
                "observation_status": machine_row["observation_status"],
                "last_seen_at": VERIFY_TIME,
                "immediate_child_count": len(children) if path.is_dir() else None,
                "recursive_size_bytes": recursive_size,
                "size_measurement_status": size_status,
                "size_blocker": size_blocker,
                "immediate_children": children,
            }
        )
    return {
        "schema_version": "1.0.0",
        "registry_type": "storage_root_immediate_child_manifest",
        "verified_at": VERIFY_TIME,
        "scan_contract": "Immediate children only; no recursive content hash or du over large roots.",
        "records": records,
    }


def build_authority_crosswalk(authority_path: Path) -> dict[str, Any]:
    repo = authority_path.resolve().parents[3]
    adjudication = repo / "docs/provenance/20260813_frozen_binary_source_adjudication_01.md"
    receipt = repo / "docs/provenance/20260813_frozen_binary_bundle_receipt.json"
    return {
        "schema_version": "1.0.0",
        "registry_type": "authority_superseded_crosswalk",
        "verified_at": VERIFY_TIME,
        "immutable_authority": {
            "path": str(authority_path.resolve()),
            "package_path": "../evidence/authority_manifest.json",
            "sha256": sha256_file(authority_path),
            "as_of_date": "2026-08-01",
            "mutation_policy": "IMMUTABLE_APPEND_ONLY_CORRECTIONS",
        },
        "crosswalks": [
            {
                "crosswalk_id": "frozen_binary_source_adjudication_20260813",
                "subject_artifact_id": "authority.frozen_binary",
                "authority_statement": "relationship_to_current_dirty_sources=not_proven_byte_reproducible",
                "adjudicated_statement": "PROVEN_BYTE_IDENTICAL_AGAINST_LL_b979760 for 5/5 solver modules",
                "relationship": "SUPERSEDES_STATUS_NOT_BYTES",
                "evidence_status": "VALIDATED_DERIVED",
                "adjudication_path": str(adjudication),
                "adjudication_package_path": "../evidence/frozen_binary_source_adjudication.md",
                "adjudication_sha256": sha256_file(adjudication),
                "receipt_path": str(receipt),
                "receipt_package_path": "../evidence/frozen_binary_bundle_receipt.json",
                "receipt_sha256": sha256_file(receipt),
                "claim_ceiling": "Source/build identity correction only; frozen scientific numbers and authority bytes remain unchanged.",
            },
            {
                "crosswalk_id": "source_snapshot_4_role_rejection_20260813",
                "subject_artifact_id": "authority.source_snapshot_4",
                "authority_statement": "Bytes registered in the 2026-08-01 19-item replay.",
                "adjudicated_statement": "Later 2026-07-26 worktree bytes; rejected as source truth for the 2026-07-24 frozen binary.",
                "relationship": "ROLE_INVALIDATED_BYTES_RETAINED",
                "evidence_status": "INVALIDATED",
                "replacement_source": "LongLineage b979760 obligation_bnb plus the five byte-identical solver modules enumerated in the adjudication receipt",
                "claim_ceiling": "Keep hash replay semantics; never connect source_snapshot_4 to authority.frozen_binary as producer input.",
            },
            {
                "crosswalk_id": "canonical_run_count_reconciliation_20260813",
                "subject_artifact_id": "claim.invalidated_canonical_19_runs",
                "authority_statement": "Historical prose called the manifest a 19-run index.",
                "adjudicated_statement": "18 logical rows, 19 TSV lines, 35 current physical directories, 51 including pending archive.",
                "relationship": "SUPERSEDED_COUNT_SEMANTICS",
                "evidence_status": "VALIDATED_DERIVED",
                "claim_ceiling": "Filesystem and registry counting only; not scientific run finality.",
            },
        ],
    }


def validate_registries(
    datasets: Sequence[Mapping[str, Any]],
    run_registry: Mapping[str, Any],
    artifacts: Sequence[Mapping[str, Any]],
    machine_registry: Mapping[str, Any],
    storage_manifest: Mapping[str, Any],
    crosswalk: Mapping[str, Any],
) -> list[str]:
    errors: list[str] = []
    records = run_registry["records"]
    for label, rows, key in (
        ("dataset", datasets, "technical_dataset_id"),
        ("physical run", records, "physical_run_id"),
        ("artifact", artifacts, "artifact_id"),
        ("machine path", machine_registry["records"], "path_id"),
        ("storage root", storage_manifest["records"], "root_id"),
        ("crosswalk", crosswalk["crosswalks"], "crosswalk_id"),
    ):
        ids = [row[key] for row in rows]
        duplicates = sorted(item for item, count in Counter(ids).items() if count > 1)
        if duplicates:
            errors.append(f"duplicate {label} ids: {duplicates}")
    if len(datasets) != 7 or len({row["biological_id"] for row in datasets}) != 6:
        errors.append("dataset scope must remain 7 technical datasets / 6 biological IDs")

    expected_reconciliation = {
        "manifest_file_lines": 19,
        "logical_manifest_rows": 18,
        "physical_directories_total": 51,
        "current_physical_directories": 35,
        "pending_archive_physical_directories": 16,
        "logical_rows_merged_current": 9,
        "logical_rows_merged_pending_archive": 9,
        "current_physical_unregistered": 26,
        "pending_archive_extras": 7,
    }
    for key, expected in expected_reconciliation.items():
        actual = run_registry["reconciliation"].get(key)
        if actual != expected:
            errors.append(f"run reconciliation {key}: expected {expected}, found {actual}")
    if len(records) != 51:
        errors.append(f"run registry must contain 51 physical records, found {len(records)}")
    if sum(bool(row["logical_manifest_member"]) for row in records) != 18:
        errors.append("logical manifest metadata must merge into exactly 18 physical records")

    tagged = [row for row in artifacts if row["artifact_type"] == "longphase_s_tagged_bam"]
    if len(tagged) != 14:
        errors.append(f"expected 14 canonical paired tagged BAM artifacts, found {len(tagged)}")
    if sum(int(row["size_bytes"] or 0) for row in tagged) != EXPECTED_TAGGED_BAM_BYTES:
        errors.append("canonical tagged BAM byte total mismatch")
    for row in tagged:
        if row["sha256"] is not None or row["scope"] != "PARTIAL" or row["finality"] != "NON_FINAL":
            errors.append(f"tagged BAM policy violation: {row['artifact_id']}")
        if row["evidence_status"] != "PARTIAL" or "PRE-FIX" not in row["claim_ceiling"]:
            errors.append(f"tagged BAM claim ceiling violation: {row['artifact_id']}")

    for row in artifacts:
        if row["scope"] not in SCOPE_ENUM or row["evidence_status"] not in EVIDENCE_ENUM:
            errors.append(f"artifact enum violation: {row['artifact_id']}")
        if row["finality"] not in FINALITY_ENUM or row["availability"] not in AVAILABILITY_ENUM:
            errors.append(f"artifact enum violation: {row['artifact_id']}")
        if row["finality"] == "FINAL_FOR_SCOPE":
            if not row["sha256"] or not row["producer"]:
                errors.append(f"final artifact lacks hash or producer: {row['artifact_id']}")
            if not row["producer_commit"] or not row["inputs"]:
                errors.append(f"final artifact lacks producer_commit or input provenance: {row['artifact_id']}")
            command = row["regeneration_command"]
            if not isinstance(command, str) or not command.startswith(REGENERATION_SEMANTIC_PREFIXES):
                errors.append(f"final artifact lacks typed replay/regeneration semantics: {row['artifact_id']}")

    source4 = next(row for row in artifacts if row["artifact_id"] == "authority.source_snapshot_4")
    frozen = next(row for row in artifacts if row["artifact_id"] == "authority.frozen_binary")
    if source4["artifact_id"] in frozen["inputs"] or frozen["artifact_id"] in source4["used_by"]:
        errors.append("source_snapshot_4 must not be connected to frozen binary as source truth")
    for row in records:
        ordered = [
            row.get("source_created_at"),
            row.get("run_started_at"),
            row.get("run_completed_at"),
            row.get("published_at"),
            row.get("last_seen_at"),
            row.get("verified_at"),
        ]
        parsed = [datetime.fromisoformat(value) for value in ordered if value]
        if any(left > right for left, right in zip(parsed, parsed[1:])):
            errors.append(f"illegal run date ordering: {row['physical_run_id']}")
    if not any(row["path_id"] == "legacy.runbook" and row["observation_status"] == "MISSING" for row in machine_registry["records"]):
        errors.append("missing runbook must be explicitly registered as MISSING")
    return errors


def schema_validator_name() -> str:
    try:
        import jsonschema  # type: ignore
    except ImportError:
        executable = shutil.which("jsonschema")
        return f"jsonschema-cli:{executable}" if executable else "UNAVAILABLE"
    return "Draft202012Validator" if hasattr(jsonschema, "Draft202012Validator") else "Draft7Validator_COMPATIBILITY_MODE"


def load_and_validate_json_schema(document: Any, schema: Mapping[str, Any]) -> list[str]:
    """Validate a registry with the Python module or the system CLI."""

    try:
        import jsonschema  # type: ignore
    except ImportError:
        executable = shutil.which("jsonschema")
        if executable is None:
            return ["schema validation requested but neither the Python jsonschema module nor jsonschema CLI is available"]
        with tempfile.TemporaryDirectory(prefix="intersubmod-schema-") as directory:
            instance_path = Path(directory) / "instance.json"
            schema_path = Path(directory) / "schema.json"
            write_json(instance_path, document)
            write_json(schema_path, schema)
            result = subprocess.run(
                [executable, "-i", str(instance_path), str(schema_path)],
                capture_output=True,
                text=True,
                check=False,
            )
        if result.returncode == 0:
            return []
        message = (result.stderr or result.stdout).strip()
        return [f"schema validation failed via {executable}: {message or 'unknown validation error'}"]
    try:
        validator_class = getattr(jsonschema, "Draft202012Validator", jsonschema.Draft7Validator)
        validator_class(schema, format_checker=jsonschema.FormatChecker()).validate(document)
    except jsonschema.ValidationError as error:
        location = getattr(error, "json_path", "/" + "/".join(str(item) for item in error.path))
        return [f"schema validation failed: {location}: {error.message}"]
    return []


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authority", required=True, type=Path)
    parser.add_argument("--master-run-manifest", required=True, type=Path)
    parser.add_argument("--canonical-root", required=True, type=Path)
    parser.add_argument("--archive-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--schema-dir", type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    manifest = json.loads(args.authority.read_text(encoding="utf-8"))
    run_rows = read_tsv(args.master_run_manifest)
    datasets = build_datasets(manifest, run_rows)
    run_registry = build_runs(run_rows, args.canonical_root, args.archive_root)
    artifacts = build_artifacts(manifest, args.authority, args.canonical_root)
    machine_registry = build_machine_paths(args.authority)
    storage_manifest = build_storage_manifest(machine_registry)
    crosswalk = build_authority_crosswalk(args.authority)
    errors = validate_registries(datasets, run_registry, artifacts, machine_registry, storage_manifest, crosswalk)

    outputs = {
        "dataset_registry.json": datasets,
        "run_registry.json": run_registry,
        "artifact_registry.json": artifacts,
        "machine_path_registry.json": machine_registry,
        "storage_root_manifest.json": storage_manifest,
        "authority_superseded_crosswalk.json": crosswalk,
    }
    if args.schema_dir:
        schema_names = {
            "dataset_registry.json": "dataset-registry.schema.json",
            "run_registry.json": "run-registry.schema.json",
            "artifact_registry.json": "artifact-registry.schema.json",
            "machine_path_registry.json": "machine-path-registry.schema.json",
            "storage_root_manifest.json": "storage-root-manifest.schema.json",
            "authority_superseded_crosswalk.json": "authority-crosswalk.schema.json",
        }
        for output_name, schema_name in schema_names.items():
            schema_path = args.schema_dir / schema_name
            if not schema_path.is_file():
                errors.append(f"missing schema: {schema_path}")
                continue
            errors.extend(load_and_validate_json_schema(outputs[output_name], json.loads(schema_path.read_text(encoding="utf-8"))))
    for name, value in outputs.items():
        write_json(args.output_dir / name, value)

    reconciliation = run_registry["reconciliation"]
    tagged = [row for row in artifacts if row["artifact_type"] == "longphase_s_tagged_bam"]
    receipt = {
        "schema_version": "2.0.0",
        "receipt_type": "handoff_registry_build",
        "verified_at": VERIFY_TIME,
        "actual_host": platform.node(),
        "datasets": len(datasets),
        "biological_ids": len({row["biological_id"] for row in datasets}),
        "run_reconciliation": reconciliation,
        "registered_physical_run_rows": len(run_registry["records"]),
        "artifact_rows": len(artifacts),
        "canonical_tagged_bam_artifacts": len(tagged),
        "canonical_tagged_bam_total_bytes": sum(row["size_bytes"] for row in tagged),
        "canonical_tagged_bam_hash_policy": "NOT_COMPUTED_TIB_PAYLOAD",
        "machine_path_rows": len(machine_registry["records"]),
        "storage_root_rows": len(storage_manifest["records"]),
        "authority_crosswalk_rows": len(crosswalk["crosswalks"]),
        "schema_validation": {
            "requested": args.schema_dir is not None,
            "validator": schema_validator_name() if args.schema_dir else "NOT_REQUESTED",
            "documents_checked": len(outputs) if args.schema_dir else 0,
        },
        "errors": errors,
        "pass": not errors,
        "claim_ceiling": "Inventory and frozen provenance validation; not a TiB-scale scientific rerun.",
    }
    write_json(args.output_dir / "registry_build_receipt.json", receipt)
    print(f"[INPUT] authority={args.authority.resolve()}")
    print(f"[INPUT] master_run_manifest={args.master_run_manifest.resolve()}")
    print(f"[INPUT] canonical_root={args.canonical_root.resolve()}")
    print(f"[INPUT] archive_root={args.archive_root.resolve()}")
    print(f"[OUTPUT] registry_dir={args.output_dir.resolve()}")
    print(
        "[RESULT] "
        f"pass={str(not errors).lower()} physical_runs={len(run_registry['records'])} "
        f"logical_merged={reconciliation['logical_manifest_rows']} tagged_bams={len(tagged)} "
        f"tagged_bam_bytes={sum(row['size_bytes'] for row in tagged)} artifacts={len(artifacts)} "
        f"machine_paths={len(machine_registry['records'])} errors={len(errors)}"
    )
    print(
        "[RECONCILIATION] "
        f"18_logical={reconciliation['logical_manifest_rows']} 19_lines={reconciliation['manifest_file_lines']} "
        f"35_current={reconciliation['current_physical_directories']} "
        f"16_pending={reconciliation['pending_archive_physical_directories']} "
        f"9_current_logical={reconciliation['logical_rows_merged_current']} "
        f"9_pending_logical={reconciliation['logical_rows_merged_pending_archive']} "
        f"26_current_unregistered={reconciliation['current_physical_unregistered']} "
        f"7_archive_extras={reconciliation['pending_archive_extras']}"
    )
    for error in errors[:10]:
        print(f"[ERROR] {error}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
