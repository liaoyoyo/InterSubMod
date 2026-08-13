#!/usr/bin/env python3
"""Run fail-closed public-claim validation and write one aggregate receipt."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from lxml import etree, html

from claim_registry_contract import (
    BUILD_RECEIPT_SCHEMA_NAME,
    BUILD_RECEIPT_SCHEMA_VERSION,
    REGISTRY_SCHEMA_NAME,
    REGISTRY_SCHEMA_VERSION,
    atomic_write_json,
    build_public_source_manifest,
    canonical_json_bytes,
    canonical_object_sha256,
    sha256_path,
)
from validate_claim_remediation_registry import validate_build_receipt


ROOT = Path(__file__).resolve().parents[3]
PROJECT = ROOT / "research/20260813_public_docs_p0_correction"
RECEIPT = PROJECT / "validation_receipt.json"
REGISTRY = PROJECT / "claim_remediation_registry.json"
BUILD_RECEIPT = PROJECT / "claim_remediation_build_receipt.json"
P0_REGISTRY = PROJECT / "scripts/p0_claim_registry.json"
BROWSER_RECEIPT = PROJECT / "html_browser_qa_receipt.json"
BROWSER_RUNNER = PROJECT / "scripts/run_html_browser_qa.py"
MINIMUM_PYTHON = (3, 10)

EXPECTED_HTML_RELATIVE_PATHS = (
    Path("docs/explain/01_background-glossary.standalone.html"),
    Path("docs/explain/02_ism-core.standalone.html"),
    Path("docs/explain/03_methylation-read-filter.standalone.html"),
    Path("docs/explain/04_subclone-reconstruction-chr2-18M.standalone.html"),
    Path("docs/explain/05_subclone-correction-audit-chr2-18M.standalone.html"),
    Path("docs/explain/06_ism-subclone-pipeline-concept.standalone.html"),
    Path("docs/explain/07_subclone-judgment-workstation-chr2-18M.standalone.html"),
    Path("docs/explain/08_subclone-logic-chain-chr2-18M.standalone.html"),
    Path("docs/explain/09_three-stats-division-of-labor.standalone.html"),
    Path("docs/explain/10_ism-cpp-vs-chr2-subclone-capability.standalone.html"),
    Path("docs/explain/11_system-map-overview.standalone.html"),
    Path("docs/explain/12_intersubmod-io.standalone.html"),
    Path("docs/explain/13_longlineage-io.standalone.html"),
    Path("docs/explain/14_upstream-data.standalone.html"),
    Path("docs/explain/15_python-html-layer.standalone.html"),
    Path("docs/explain/16_how-to-run.standalone.html"),
    Path("docs/explain/index.standalone.html"),
    Path(
        "docs/handoff/20260813_完整研究資料與軟體交接_01/"
        "20260813_完整研究交接總覽_01.standalone.html"
    ),
    Path(
        "docs/handoff/20260805_系統交接與驗收_01/"
        "00_交接總覽與驗收.standalone.html"
    ),
    Path(
        "docs/handoff/20260805_系統交接與驗收_01/"
        "05_repo整理與可觀察性盤點.standalone.html"
    ),
    Path("docs/handoff/20260806_兩repo端到端串接可行性盤點_01.standalone.html"),
)
EXPECTED_BROWSER_MODES = {
    "desktop": {"javascript_enabled": True, "print_media": False},
    "mobile": {"javascript_enabled": True, "print_media": False},
    "no_js": {"javascript_enabled": False, "print_media": False},
    "print": {"javascript_enabled": True, "print_media": True},
}
EXPECTED_BROWSER_CASES = len(EXPECTED_HTML_RELATIVE_PATHS) * len(EXPECTED_BROWSER_MODES)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rebuild and validate the 158-claim remediation registry, or run a "
            "side-effect-free validation of the currently materialized artifacts."
        )
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Do not rebuild registries, run Chromium, or write validation_receipt.json.",
    )
    parser.add_argument(
        "--run-browser-qa",
        action="store_true",
        help="Regenerate the exact 84-case Chromium receipt before validating it.",
    )
    args = parser.parse_args()
    if args.check_only and args.run_browser_qa:
        parser.error("--check-only and --run-browser-qa are mutually exclusive")
    return args


def require_supported_python(version_info: Any = None) -> str:
    """Return the interpreter version or raise when it is below Python 3.10."""
    version_info = sys.version_info if version_info is None else version_info
    version = tuple(version_info[:3])
    if version < (*MINIMUM_PYTHON, 0):
        raise RuntimeError(
            f"Python >= {MINIMUM_PYTHON[0]}.{MINIMUM_PYTHON[1]} is required; "
            f"running {version[0]}.{version[1]}.{version[2]}"
        )
    return ".".join(str(item) for item in version)


def run(command: list[str], root: Path = ROOT) -> dict[str, object]:
    result = subprocess.run(command, cwd=root, text=True, capture_output=True, check=False)
    return {
        "command": command,
        "exit_code": result.returncode,
        "stdout_tail": result.stdout[-4000:],
        "stderr_tail": result.stderr[-4000:],
    }


def resolve_html_allowlist(root: Path) -> list[Path]:
    """Resolve the exact 21-file contract and reject missing or extra pages."""
    expected = {root / relative for relative in EXPECTED_HTML_RELATIVE_PATHS}
    discovered = set((root / "docs/explain").glob("*.standalone.html"))
    discovered.update((root / "docs/handoff").rglob("*.standalone.html"))
    missing = sorted(expected - discovered)
    unexpected = sorted(discovered - expected)
    if missing or unexpected:
        details = []
        if missing:
            details.append("missing=" + ",".join(str(path.relative_to(root)) for path in missing))
        if unexpected:
            details.append(
                "unexpected=" + ",".join(str(path.relative_to(root)) for path in unexpected)
            )
        raise RuntimeError("standalone HTML allowlist mismatch: " + "; ".join(details))
    return [root / relative for relative in EXPECTED_HTML_RELATIVE_PATHS]


def html_checks(root: Path = ROOT) -> dict[str, object]:
    errors: list[str] = []
    svg_count = 0
    external_resources: list[str] = []
    resource_tags = {"script": "src", "img": "src", "iframe": "src", "source": "src", "link": "href"}
    try:
        files = resolve_html_allowlist(root)
    except RuntimeError as exc:
        files = []
        errors.append(str(exc))
    for path in files:
        relative = path.relative_to(root)
        try:
            document = html.fromstring(path.read_bytes())
        except (OSError, etree.ParserError, ValueError) as exc:
            errors.append(f"{relative}: HTML parse failed: {exc}")
            continue
        svgs = document.xpath("//*[local-name()='svg']")
        svg_count += len(svgs)
        for index, node in enumerate(svgs):
            try:
                etree.fromstring(etree.tostring(node))
            except etree.XMLSyntaxError as exc:
                errors.append(f"{relative}: svg[{index}] XML parse failed: {exc}")
        for tag, attr in resource_tags.items():
            for node in document.xpath(f"//{tag}[@{attr}]"):
                value = node.get(attr, "")
                if value.startswith(("http://", "https://", "//")):
                    external_resources.append(f"{relative}:{tag}:{value}")
    return {
        "expected_files": len(EXPECTED_HTML_RELATIVE_PATHS),
        "files": len(files),
        "paths": [str(path.relative_to(root)) for path in files],
        "inline_svg_elements": svg_count,
        "xml_svg_parse_errors": errors,
        "external_runtime_resources_declared": external_resources,
        "verdict": (
            "PASS"
            if len(files) == len(EXPECTED_HTML_RELATIVE_PATHS) and not errors and not external_resources
            else "FAIL"
        ),
    }


def _read_json_object(
    path: Path, *, require_canonical: bool = True
) -> tuple[dict[str, Any], bytes, list[str]]:
    errors: list[str] = []
    try:
        raw = path.read_bytes()
        value = json.loads(raw.decode("utf-8"))
        if not isinstance(value, dict):
            raise TypeError("root must be an object")
        if require_canonical and raw != canonical_json_bytes(value):
            errors.append(f"{path}: bytes are not canonical UTF-8 JSON")
        return value, raw, errors
    except (OSError, UnicodeError, json.JSONDecodeError, TypeError) as exc:
        return {}, b"", [f"{path}: cannot read canonical JSON: {exc}"]


def _resolve_repo_file(root: Path, relative: object, label: str, errors: list[str]) -> Path | None:
    if not isinstance(relative, str) or not relative:
        errors.append(f"{label}: path must be a non-empty string")
        return None
    path = (root / relative).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError:
        errors.append(f"{label}: path escapes repository: {relative}")
        return None
    if not path.is_file():
        errors.append(f"{label}: file is missing: {relative}")
        return None
    return path


def validate_claim_registry_bundle(
    root: Path = ROOT,
    registry_path: Path = REGISTRY,
    receipt_path: Path = BUILD_RECEIPT,
    p0_path: Path = P0_REGISTRY,
) -> dict[str, Any]:
    """Validate the registry/build-receipt contract against current source bytes."""
    registry, registry_bytes, errors = _read_json_object(registry_path)
    receipt, _receipt_bytes, receipt_errors = _read_json_object(receipt_path)
    errors.extend(receipt_errors)
    p0, _p0_bytes, p0_errors = _read_json_object(p0_path, require_canonical=False)
    errors.extend(p0_errors)
    if not errors:
        if registry.get("schema_name") != REGISTRY_SCHEMA_NAME:
            errors.append("claim registry schema_name mismatch")
        if registry.get("schema_version") != REGISTRY_SCHEMA_VERSION:
            errors.append("claim registry schema_version mismatch")
        if receipt.get("schema_name") != BUILD_RECEIPT_SCHEMA_NAME:
            errors.append("claim build receipt schema_name mismatch")
        if receipt.get("schema_version") != BUILD_RECEIPT_SCHEMA_VERSION:
            errors.append("claim build receipt schema_version mismatch")
        errors.extend(validate_build_receipt(registry, registry_bytes, receipt))

        for path_key, hash_key in (
            ("source_inventory", "source_inventory_sha256"),
            ("source_scope", "source_scope_sha256"),
            ("p0_guard_registry", "p0_guard_registry_sha256"),
            ("github_about_receipt", "github_about_receipt_sha256"),
        ):
            source = _resolve_repo_file(root, registry.get(path_key), path_key, errors)
            if source is not None and registry.get(hash_key) != sha256_path(source):
                errors.append(f"{hash_key} does not bind current {path_key} bytes")

        try:
            current_manifest = build_public_source_manifest(root, p0)
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            current_manifest = {}
            errors.append(f"public source manifest cannot be rebuilt: {exc}")
        if registry.get("public_source_manifest") != current_manifest:
            errors.append("claim registry public source manifest is stale")
        if registry.get("public_source_manifest_sha256") != canonical_object_sha256(current_manifest):
            errors.append("claim registry public source manifest hash is stale")

    return {
        "registry": str(registry_path.relative_to(root)),
        "build_receipt": str(receipt_path.relative_to(root)),
        "registry_sha256": sha256_path(registry_path) if registry_path.is_file() else None,
        "public_source_files": registry.get("public_source_manifest", {}).get("counts", {}).get("files"),
        "errors": errors,
        "verdict": "PASS" if not errors else "FAIL",
    }


def _version_tuple(value: object) -> tuple[int, int, int] | None:
    if not isinstance(value, str):
        return None
    match = re.fullmatch(r"(\d+)\.(\d+)(?:\.(\d+))?.*", value)
    if not match:
        return None
    return tuple(int(item or 0) for item in match.groups())


def _check_source_records(
    records: object,
    expected_paths: Iterable[Path],
    root: Path,
    label: str,
    errors: list[str],
) -> None:
    expected = list(expected_paths)
    if not isinstance(records, list):
        errors.append(f"browser receipt {label} source provenance must be a list")
        return
    actual_paths = [item.get("path") if isinstance(item, dict) else None for item in records]
    wanted_paths = [str(path.relative_to(root)) for path in expected]
    if actual_paths != wanted_paths:
        errors.append(f"browser receipt {label} source path set/order mismatch")
        return
    for item, path in zip(records, expected):
        if item.get("sha256") != sha256_path(path):
            errors.append(f"browser receipt {label} source hash is stale: {path.relative_to(root)}")


def validate_browser_qa_receipt(
    root: Path = ROOT,
    receipt_path: Path = BROWSER_RECEIPT,
    runner_path: Path = BROWSER_RUNNER,
) -> dict[str, Any]:
    """Bind the browser receipt to the current 21 HTML, all SVG, and runner bytes."""
    receipt, _raw, errors = _read_json_object(receipt_path)
    try:
        html_files = resolve_html_allowlist(root)
    except RuntimeError as exc:
        html_files = []
        errors.append(str(exc))
    svg_files = sorted((root / "docs").rglob("*.svg"))
    runner_relative = str(runner_path.relative_to(root))

    if receipt.get("schema_version") != 3:
        errors.append("browser receipt schema_version must be 3")
    if receipt.get("verdict") != "PASS_FOR_LOCAL_SOURCE_QA":
        errors.append("browser receipt verdict is not PASS_FOR_LOCAL_SOURCE_QA")
    runner = receipt.get("runner")
    if not isinstance(runner, dict):
        errors.append("browser receipt runner must be an object")
        runner = {}
    if runner.get("path") != runner_relative:
        errors.append("browser receipt runner path mismatch")
    if not runner_path.is_file() or runner.get("sha256") != sha256_path(runner_path):
        errors.append("browser receipt runner hash is stale")
    runner_python = _version_tuple(runner.get("python"))
    if runner_python is None or runner_python < (3, 10, 0):
        errors.append("browser QA must be executed with Python >= 3.10")

    provenance = receipt.get("source_provenance")
    if not isinstance(provenance, dict):
        errors.append("browser receipt source_provenance must be an object")
        provenance = {}
    _check_source_records(provenance.get("html"), html_files, root, "HTML", errors)
    _check_source_records(provenance.get("standalone_svg"), svg_files, root, "SVG", errors)
    runner_record = provenance.get("runner")
    if not isinstance(runner_record, dict):
        errors.append("browser receipt runner source provenance is missing")
    elif (
        runner_record.get("path") != runner_relative
        or runner_record.get("sha256") != sha256_path(runner_path)
    ):
        errors.append("browser receipt runner source provenance is stale")
    if provenance.get("stable_during_run") is not True:
        errors.append("browser receipt sources were not stable during the run")

    tested_files = receipt.get("tested_files")
    _check_source_records(tested_files, html_files, root, "tested HTML", errors)
    if isinstance(tested_files, list):
        if any(
            not isinstance(item, dict) or item.get("passed") is not True or item.get("errors") != []
            for item in tested_files
        ):
            errors.append("browser receipt HTML parsing contains failures")

    expected_combinations = {
        (str(path.relative_to(root)), mode)
        for path in html_files
        for mode in EXPECTED_BROWSER_MODES
    }
    cases = receipt.get("browser_cases")
    actual_combinations: list[tuple[object, object]] = []
    if not isinstance(cases, list):
        errors.append("browser receipt browser_cases must be a list")
        cases = []
    for case in cases:
        if not isinstance(case, dict):
            errors.append("browser receipt contains a non-object browser case")
            continue
        combination = (case.get("path"), case.get("mode"))
        actual_combinations.append(combination)
        mode_contract = EXPECTED_BROWSER_MODES.get(str(case.get("mode")))
        source_path = root / str(case.get("path"))
        if combination not in expected_combinations:
            errors.append(f"unexpected browser case: {combination}")
            continue
        if case.get("source_sha256") != sha256_path(source_path):
            errors.append(f"browser case source hash is stale: {combination}")
        if not mode_contract or any(case.get(key) != value for key, value in mode_contract.items()):
            errors.append(f"browser case mode contract mismatch: {combination}")
        if (
            case.get("passed") is not True
            or case.get("document_overflow_pass") is not True
            or case.get("page_errors") != []
            or case.get("console_errors") != []
            or case.get("external_runtime_requests") != []
            or case.get("navigation_error") is not None
        ):
            errors.append(f"browser case did not pass all runtime checks: {combination}")
        if case.get("mode") == "print" and (
            not isinstance(case.get("print_pdf_size_bytes"), int)
            or case["print_pdf_size_bytes"] <= 1000
            or not isinstance(case.get("print_pdf_sha256"), str)
        ):
            errors.append(f"browser print case has no in-memory PDF evidence: {combination}")
    if len(actual_combinations) != EXPECTED_BROWSER_CASES:
        errors.append(
            f"browser receipt must contain exactly {EXPECTED_BROWSER_CASES} cases, "
            f"got {len(actual_combinations)}"
        )
    if len(set(actual_combinations)) != len(actual_combinations):
        errors.append("browser receipt contains duplicate cases")
    if set(actual_combinations) != expected_combinations:
        errors.append("browser receipt does not cover the exact 21 x 4 matrix")

    checks = receipt.get("checks")
    browser_counts = checks.get("browser_cases") if isinstance(checks, dict) else None
    expected_per_mode = {
        mode: len(EXPECTED_HTML_RELATIVE_PATHS) for mode in sorted(EXPECTED_BROWSER_MODES)
    }
    if not isinstance(browser_counts, dict) or {
        "expected": browser_counts.get("expected"),
        "total": browser_counts.get("total"),
        "per_mode": browser_counts.get("per_mode"),
        "passed": browser_counts.get("passed"),
        "failed": browser_counts.get("failed"),
    } != {
        "expected": EXPECTED_BROWSER_CASES,
        "total": EXPECTED_BROWSER_CASES,
        "per_mode": expected_per_mode,
        "passed": EXPECTED_BROWSER_CASES,
        "failed": 0,
    }:
        errors.append(f"browser receipt aggregate {EXPECTED_BROWSER_CASES}-case counters mismatch")
    if isinstance(checks, dict):
        if checks.get("page_errors") != 0 or checks.get("console_errors") != 0:
            errors.append("browser receipt aggregate runtime errors are non-zero")
        external = checks.get("external_runtime_requests")
        if not isinstance(external, dict) or external.get("attempted") != 0:
            errors.append("browser receipt aggregate external request count is non-zero")
        if checks.get("source_provenance_stable_during_run") is not True:
            errors.append("browser receipt aggregate source stability check did not pass")

    return {
        "receipt": str(receipt_path.relative_to(root)),
        "receipt_sha256": sha256_path(receipt_path) if receipt_path.is_file() else None,
        "html_files": len(html_files),
        "standalone_svg_files": len(svg_files),
        "browser_cases": len(cases),
        "browser_case_modes": dict(
            sorted(
                Counter(
                    str(item.get("mode")) for item in cases if isinstance(item, dict)
                ).items()
            )
        ),
        "errors": errors,
        "verdict": "PASS" if not errors else "FAIL",
    }


def main() -> int:
    try:
        python_version = require_supported_python()
    except RuntimeError as exc:
        print(json.dumps({"verdict": "FAIL", "error": str(exc)}, ensure_ascii=False))
        return 2
    args = parse_args()
    commands = []
    if not args.check_only:
        commands.append(
            run([sys.executable, str(PROJECT / "scripts/build_claim_remediation_registry.py")])
        )
    if args.run_browser_qa:
        commands.append(
            run(
                [
                    sys.executable,
                    str(BROWSER_RUNNER.relative_to(ROOT)),
                    "--root",
                    str(ROOT),
                    "--output",
                    str(BROWSER_RECEIPT),
                ]
            )
        )
    commands.extend(
        [
            run([sys.executable, str(PROJECT / "scripts/validate_public_p0_claims.py")]),
            run([sys.executable, str(PROJECT / "scripts/validate_claim_remediation_registry.py")]),
            run(
                [
                    sys.executable,
                    "-m",
                    "pytest",
                    "-q",
                    "research/20260813_public_docs_p0_correction/tests",
                ]
            ),
        ]
    )
    html_result = html_checks()
    registry_result = validate_claim_registry_bundle()
    browser_result = validate_browser_qa_receipt()
    passed = bool(
        all(item["exit_code"] == 0 for item in commands)
        and html_result["verdict"] == "PASS"
        and registry_result["verdict"] == "PASS"
        and browser_result["verdict"] == "PASS"
    )
    receipt = {
        "schema_name": "intersubmod.public_claim_validation_receipt",
        "schema_version": "2.0.0",
        "created_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "task_type": ["B_COMPREHENSIVE_VALIDATION", "D_EXTERNAL_HANDOFF"],
        "scope": "LOCAL_SOURCE_PLUS_C108_LIVE_RECEIPT_PLUS_CHROMIUM_QA",
        "verdict": "PASS" if passed else "FAIL",
        "publication_status": "BLOCKED_DEFAULT_BRANCH_WIKI_PAGES_NOT_PUBLISHED_OR_REFETCHED__ABOUT_C108_CONFIRMED",
        "release_status": "BLOCKED",
        "runner": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256_path(Path(__file__).resolve()),
            "python": python_version,
        },
        "commands": commands,
        "claim_registry_contract": registry_result,
        "html_source_checks": html_result,
        "browser_qa_contract": browser_result,
        "authority_chain": {
            "policy_input": str(P0_REGISTRY.relative_to(ROOT)),
            "builder": str((PROJECT / "scripts/build_claim_remediation_registry.py").relative_to(ROOT)),
            "generated_registry": str(REGISTRY.relative_to(ROOT)),
            "build_receipt": str(BUILD_RECEIPT.relative_to(ROOT)),
            "browser_receipt": str(BROWSER_RECEIPT.relative_to(ROOT)),
            "rebuild_executed_first": not args.check_only,
            "browser_qa_executed_first": args.run_browser_qa,
        },
        "claim_ceiling": {
            "confirmed_cellular_subclone": 0,
            "confirmed_linear_ancestry": 0,
            "methylation": "ASSOCIATION_ONLY",
            "cn_loh": "NOT_INTEGRATED",
            "graph_shape_88_2579_percent": "MODEL_CONDITIONAL_NOT_BIOLOGICAL_PREVALENCE",
        },
    }
    if not args.check_only:
        atomic_write_json(RECEIPT, receipt)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
