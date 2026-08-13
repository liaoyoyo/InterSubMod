#!/usr/bin/env python3
"""Fail-closed validator for the public tiny InterSubMod E2E fixture."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import socket
import subprocess
import sys
from urllib.parse import urlsplit, urlunsplit
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence


EXPECTED_READ_NAMES = tuple(
    [f"tiny_hp1_{index:02d}" for index in range(1, 7)]
    + [f"tiny_hp2_{index:02d}" for index in range(1, 7)]
)
TREE_SEMANTICS = "read_dendrogram_from_methylation_distance_not_cellular_lineage"
CLAIM_CEILING = "DEMO I/O and schema reproducibility only; no biological validation"


def public_repository_locator(value: Optional[str]) -> Optional[str]:
    """Return a receipt-safe repository locator without URL credentials or tracking data."""

    if not value or "://" not in value:
        return value
    parsed = urlsplit(value)
    hostname = parsed.hostname or ""
    if parsed.port:
        hostname = f"{hostname}:{parsed.port}"
    return urlunsplit((parsed.scheme, hostname, parsed.path, "", ""))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add_check(
    checks: List[Dict[str, Any]],
    check_id: str,
    passed: bool,
    observed: Any,
    expected: Any,
) -> None:
    checks.append(
        {
            "check_id": check_id,
            "status": "PASS" if passed else "FAIL",
            "observed": observed,
            "expected": expected,
        }
    )


def read_json(path: Path) -> Dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def validate_csv(summary_path: Path, schema_path: Path, checks: List[Dict[str, Any]]) -> Dict[str, Any]:
    schema = read_json(schema_path)
    columns = schema.get("columns")
    count = schema.get("column_count")
    schema_shape_ok = (
        count == 199
        and isinstance(columns, list)
        and len(columns) == 199
        and len(set(columns)) == 199
    )
    add_check(checks, "SCHEMA_DECLARATION_199_UNIQUE", schema_shape_ok,
              {"column_count": count, "listed": len(columns or []), "unique": len(set(columns or []))},
              {"column_count": 199, "listed": 199, "unique": 199})

    with summary_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    header = rows[0] if rows else []
    data_rows = rows[1:]
    add_check(checks, "SUMMARY_HEADER_EXACT", header == columns,
              {"count": len(header), "sha256": hashlib.sha256(",".join(header).encode()).hexdigest()},
              {"count": 199, "columns": "exact schema order"})
    row_widths = [len(row) for row in data_rows]
    add_check(checks, "SUMMARY_ONE_ROW_WIDTH_199", len(data_rows) == 1 and row_widths == [199],
              {"rows": len(data_rows), "widths": row_widths}, {"rows": 1, "widths": [199]})

    indexed: Dict[str, str] = {}
    if len(data_rows) == 1 and len(header) == len(data_rows[0]):
        indexed = dict(zip(header, data_rows[0]))
    expected_values = {
        "Chr": "chrTiny",
        "Pos": "101",
        "NumReads": "12",
        "NumCpGs": "50",
        "NReadsValid": "12",
        "VerificationSchemaVersion": "2",
        "RegionStratificationSchemaVersion": "1",
        "RegionStratum_ID": "-1",
        "RegionStratum_Label": "Unassigned",
        "RegionStratum_Reason": "NOT_APPLICABLE_TUMOR_ONLY",
        "Subclone_ID": "-1",
    }
    observed_values = {key: indexed.get(key) for key in expected_values}
    add_check(checks, "SUMMARY_TINY_ROW_CONTRACT", observed_values == expected_values,
              observed_values, expected_values)
    return {
        "column_count": len(header),
        "row_count": len(data_rows),
        "required_values": observed_values,
    }


def newick_leaf_names(text: str) -> List[str]:
    return re.findall(r"(?<=[(,])([^():,;]+)(?=:)", text)


def validate_tree(output_dir: Path, checks: List[Dict[str, Any]]) -> Dict[str, Any]:
    paths = sorted(output_dir.glob("**/clustering/tree.nwk"))
    add_check(checks, "ONE_MAIN_READ_DENDROGRAM", len(paths) == 1,
              [str(path.relative_to(output_dir)) for path in paths], "exactly one clustering/tree.nwk")
    if len(paths) != 1:
        return {"path": None, "leaf_count": 0, "semantics": TREE_SEMANTICS}

    text = paths[0].read_text(encoding="utf-8").strip()
    names = newick_leaf_names(text)
    expected = sorted(EXPECTED_READ_NAMES)
    add_check(checks, "TREE_VALID_NEWICK_ENVELOPE", text.startswith("(") and text.endswith(";"),
              {"starts_paren": text.startswith("("), "ends_semicolon": text.endswith(";")}, True)
    add_check(checks, "TREE_EXACT_SYNTHETIC_READ_LEAVES", sorted(names) == expected,
              {"leaf_count": len(names), "leaf_names": sorted(names)},
              {"leaf_count": 12, "leaf_names": expected})
    return {
        "path": str(paths[0]),
        "sha256": sha256(paths[0]),
        "leaf_count": len(names),
        "semantics": TREE_SEMANTICS,
        "not_evidence_for": [
            "cellular_lineage",
            "linear_ancestry",
            "subclone_prevalence",
            "biological_accuracy",
        ],
    }


def git_metadata(
    repo_root: Path,
    *,
    source_repository: Optional[str] = None,
    requested_revision: Optional[str] = None,
    resolved_commit: Optional[str] = None,
) -> Dict[str, Any]:
    def git(*args: str) -> str:
        return subprocess.run(
            ["git", "-C", str(repo_root), *args],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout.strip()

    try:
        status = git("status", "--porcelain")
        commit = git("rev-parse", "HEAD")
        return {
            "source_repository": public_repository_locator(source_repository),
            "requested_revision": requested_revision,
            "resolved_commit": resolved_commit or commit,
            "commit": commit,
            "clean": status == "",
            "status_line_count": len(status.splitlines()),
        }
    except (OSError, subprocess.CalledProcessError):
        return {
            "source_repository": public_repository_locator(source_repository),
            "requested_revision": requested_revision,
            "resolved_commit": resolved_commit,
            "commit": "unknown",
            "clean": False,
            "status_line_count": None,
        }


def validate(args: argparse.Namespace) -> Dict[str, Any]:
    output_dir = args.output_dir.resolve()
    fixture_dir = args.fixture_dir.resolve()
    schema_path = args.schema.resolve()
    checks: List[Dict[str, Any]] = []
    add_check(checks, "CORE_PROCESS_EXIT_ZERO", args.run_exit_code == 0, args.run_exit_code, 0)

    required = (
        "significance_summary.csv",
        "significance_statistics.txt",
        "run_params.json",
        "run_summary.json",
        "region_stratification_status.tsv",
    )
    missing = [name for name in required if not (output_dir / name).is_file()]
    add_check(checks, "REQUIRED_TOP_LEVEL_ARTIFACTS", not missing, missing, [])

    fixture_required = ("reference.fa", "reference.fa.fai", "variants.vcf", "tumor.bam", "tumor.bam.bai")
    fixture_missing = [name for name in fixture_required if not (fixture_dir / name).is_file()]
    add_check(checks, "FIXTURE_RUNTIME_FILES", not fixture_missing, fixture_missing, [])

    csv_summary: Dict[str, Any] = {}
    if not missing and schema_path.is_file():
        csv_summary = validate_csv(output_dir / "significance_summary.csv", schema_path, checks)
    else:
        add_check(checks, "SUMMARY_VALIDATION_PREREQUISITES", False,
                  {"missing_outputs": missing, "schema_exists": schema_path.is_file()},
                  {"missing_outputs": [], "schema_exists": True})

    run_summary: Dict[str, Any] = {}
    if (output_dir / "run_summary.json").is_file():
        run_summary = read_json(output_dir / "run_summary.json")
        observed = {
            "regions": run_summary.get("regions"),
            "reads": run_summary.get("reads", {}).get("total"),
            "cpg_sites": run_summary.get("cpg_sites", {}).get("total"),
            "valid_pairs": run_summary.get("distance_matrix", {}).get("valid_pairs"),
            "invalid_pairs": run_summary.get("distance_matrix", {}).get("invalid_pairs_insufficient_overlap"),
        }
        expected = {
            "regions": {"total": 1, "succeeded": 1, "failed": 0},
            "reads": 12,
            "cpg_sites": 50,
            "valid_pairs": 66,
            "invalid_pairs": 0,
        }
        add_check(checks, "RUN_SUMMARY_EXACT_TINY_COUNTS", observed == expected, observed, expected)

    tree = validate_tree(output_dir, checks)
    source = git_metadata(
        args.repo_root.resolve(),
        source_repository=args.source_repository,
        requested_revision=args.requested_revision,
        resolved_commit=args.resolved_commit,
    )
    if args.source_repository:
        revision_bound = bool(
            args.requested_revision
            and args.resolved_commit
            and source.get("commit") == args.resolved_commit
        )
        add_check(
            checks,
            "SOURCE_REVISION_RESOLVED_COMMIT_BOUND",
            revision_bound,
            {
                "requested_revision": args.requested_revision,
                "declared_resolved_commit": args.resolved_commit,
                "checkout_commit": source.get("commit"),
            },
            "requested revision present and declared resolved commit equals checkout HEAD",
        )

    hashed_files: Dict[str, Dict[str, Any]] = {}
    source_fixture_dir = args.repo_root.resolve() / "tests" / "fixtures" / "tiny_public"
    source_fixture_files = (
        source_fixture_dir / "reference.fa",
        source_fixture_dir / "variants.vcf",
        source_fixture_dir / "tumor.sam",
        source_fixture_dir / "expected_significance_schema.json",
    )
    for path in [schema_path, *source_fixture_files, *(fixture_dir / name for name in fixture_required), *(output_dir / name for name in required)]:
        if path.is_file():
            hashed_files[str(path)] = {"size_bytes": path.stat().st_size, "sha256": sha256(path)}
    if args.binary and args.binary.is_file():
        hashed_files[str(args.binary.resolve())] = {
            "size_bytes": args.binary.stat().st_size,
            "sha256": sha256(args.binary),
        }
    if args.run_log and args.run_log.is_file():
        hashed_files[str(args.run_log.resolve())] = {
            "size_bytes": args.run_log.stat().st_size,
            "sha256": sha256(args.run_log),
        }
    if args.command_receipt and args.command_receipt.is_file():
        hashed_files[str(args.command_receipt.resolve())] = {
            "size_bytes": args.command_receipt.stat().st_size,
            "sha256": sha256(args.command_receipt),
        }

    binary_identity: Dict[str, Any] = {"path": None, "size_bytes": None, "sha256": None}
    if args.binary and args.binary.is_file():
        binary_identity = {
            "path": str(args.binary.resolve()),
            "size_bytes": args.binary.stat().st_size,
            "sha256": sha256(args.binary),
        }
    add_check(
        checks,
        "BINARY_IDENTITY_BOUND",
        bool(binary_identity["sha256"] and binary_identity["size_bytes"]),
        binary_identity,
        "non-empty executable path with size and SHA-256",
    )
    all_pass = all(check["status"] == "PASS" for check in checks)

    return {
        "schema_name": "intersubmod.tiny_public_e2e_receipt",
        "schema_version": "1.0.0",
        "verified_at": datetime.now(timezone.utc).isoformat(),
        "hostname": socket.gethostname(),
        "all_pass": all_pass,
        "clean_clone_e2e_pass": bool(
            all_pass and args.source_repository and source.get("clean") is True
        ),
        "overall_release_gate": "NOT_EVALUATED",
        "scope": "DEMO",
        "evidence_status": "VALIDATED_DERIVED" if all_pass else "PARTIAL",
        "claim_ceiling": CLAIM_CEILING,
        "source_checkout": source,
        "binary_identity": binary_identity,
        "paths": {
            "repo_root": str(args.repo_root.resolve()),
            "fixture_dir": str(fixture_dir),
            "output_dir": str(output_dir),
            "schema": str(schema_path),
        },
        "summary": csv_summary,
        "tree": tree,
        "checks": checks,
        "hashed_files": hashed_files,
    }


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--fixture-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--schema", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--run-log", type=Path)
    parser.add_argument("--command-receipt", type=Path)
    parser.add_argument("--run-exit-code", type=int, required=True)
    parser.add_argument("--source-repository")
    parser.add_argument("--requested-revision")
    parser.add_argument("--resolved-commit")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    try:
        receipt = validate(args)
    except Exception as exc:  # Fail closed while preserving a diagnostic receipt.
        receipt = {
            "schema_name": "intersubmod.tiny_public_e2e_receipt",
            "schema_version": "1.0.0",
            "verified_at": datetime.now(timezone.utc).isoformat(),
            "all_pass": False,
            "clean_clone_e2e_pass": False,
            "overall_release_gate": "NOT_EVALUATED",
            "scope": "DEMO",
            "evidence_status": "PARTIAL",
            "claim_ceiling": CLAIM_CEILING,
            "fatal_error": f"{type(exc).__name__}: {exc}",
        }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        f"TINY_E2E {'PASS' if receipt.get('all_pass') else 'FAIL'} "
        f"clean_clone_e2e_pass={str(receipt.get('clean_clone_e2e_pass')).lower()} receipt={args.receipt}"
    )
    return 0 if receipt.get("all_pass") else 1


if __name__ == "__main__":
    raise SystemExit(main())
