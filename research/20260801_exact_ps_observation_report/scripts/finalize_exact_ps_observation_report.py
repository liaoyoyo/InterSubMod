#!/usr/bin/env python3
"""Build, browser-QA and publish one immutable Exact-PS observation report release."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import secrets
import shlex
import stat
import subprocess
import sys
from typing import Any


TOPIC_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = TOPIC_ROOT.parents[1]
DEFAULT_BUILDER = REPO_ROOT / "scripts/analysis/build_exact_ps_observation_report.py"
DEFAULT_QA = TOPIC_ROOT / "scripts/qa_exact_ps_observation_report.py"
SCHEMA_NAME = "intersubmod.exact_ps_observation_report.release_receipt"
SCHEMA_VERSION = "1.0.0"
AUTHORITY_SCHEMA = "intersubmod.exact_ps_tree_research.ai_handoff_authority"
BUILD_SCHEMA = "intersubmod.exact_ps_observation_report.build_receipt"
BUILD_VERSION = "1.0.0"
QA_SCHEMA = "intersubmod.exact_ps_observation_report.browser_qa"
QA_VERSION = "1.0.0"
REPORT_DATA_SCHEMA = "intersubmod.exact_ps_observation_report.data"
REPORT_DATA_VERSION = "1.0.0"
EXPECTED_ARTIFACT_IDS = {
    "strict_linkage_ready",
    "all7_cohort_receipt",
    "all7_summary",
    "all7_cpp_reader_report",
    "topology_receipt",
    "topology_summary",
    "topology_reader_report",
    "candidate_factorization_receipt",
    "read_af_decision_report",
    "methyl_manifest",
    "methyl_sidecar_receipt",
    "methyl_report_data",
    "methyl_reader_report",
}
EXPECTED_BUILD_FILES = {
    "report_data.json",
    "report_data.json.sha256",
    "report.standalone.html",
    "report.standalone.html.sha256",
    "report_build_receipt.json",
    "report_build_receipt.json.sha256",
}
EXPECTED_QA_CHECKS = {
    "report_data_schema_pass",
    "report_data_contract_pass",
    "html_data_binding_pass",
    "sample_identity_and_binding_pass",
    "four_viewports_pass",
    "no_js_core_content_pass",
    "print_pdf_pass",
    "print_page_and_overflow_pass",
    "console_errors_zero",
    "page_errors_zero",
    "external_requests_zero",
    "sample_interaction_pass",
}
RELEASE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


class FinalizeError(RuntimeError):
    """A release-blocking orchestration or evidence-contract failure."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FinalizeError(message)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def timestamp_slug() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    require(not path.is_symlink(), f"symlink is not allowed: {path}")
    resolved = path.resolve(strict=True)
    require(stat.S_ISREG(resolved.stat().st_mode), f"not a regular file: {path}")
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_path(resolved),
    }


def read_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing JSON: {path}")
    try:
        payload = json.loads(
            path.read_text(encoding="utf-8"),
            parse_constant=lambda value: (_ for _ in ()).throw(
                ValueError(f"non-finite JSON constant: {value}")
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FinalizeError(f"cannot read JSON {path}: {exc}") from exc
    require(isinstance(payload, dict), f"JSON root is not an object: {path}")
    return payload


def write_json_exclusive(path: Path, payload: dict[str, Any]) -> None:
    text = json.dumps(
        payload,
        ensure_ascii=False,
        allow_nan=False,
        indent=2,
        sort_keys=True,
    ) + "\n"
    try:
        with path.open("x", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
    except FileExistsError as exc:
        raise FinalizeError(f"refusing to overwrite release receipt: {path}") from exc


def write_sidecar_exclusive(path: Path) -> Path:
    sidecar = path.with_name(f"{path.name}.sha256")
    try:
        with sidecar.open("x", encoding="ascii") as handle:
            handle.write(f"{sha256_path(path)}  {path.name}\n")
            handle.flush()
            os.fsync(handle.fileno())
    except FileExistsError as exc:
        raise FinalizeError(f"refusing to overwrite checksum: {sidecar}") from exc
    return sidecar


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--authority-manifest", type=Path, required=True)
    parser.add_argument("--denominator-registry", type=Path, required=True)
    parser.add_argument("--release-root", type=Path, required=True)
    parser.add_argument("--release-id", required=True)
    parser.add_argument("--builder", type=Path, default=DEFAULT_BUILDER)
    parser.add_argument("--qa-script", type=Path, default=DEFAULT_QA)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--chromium", type=Path)
    parser.add_argument("--command-timeout-seconds", type=int, default=600)
    parser.add_argument(
        "--allow-nonproduction-tools",
        action="store_true",
        help="Permit custom builder/QA paths for contract tests only.",
    )
    return parser.parse_args()


def run_command(
    command: list[str], label: str, timeout_seconds: int
) -> dict[str, Any]:
    started = utc_now()
    try:
        completed = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=timeout_seconds,
        )
        returncode = completed.returncode
        stdout = completed.stdout
        stderr = completed.stderr
    except subprocess.TimeoutExpired as exc:
        returncode = -9
        stdout = exc.stdout if isinstance(exc.stdout, str) else ""
        stderr = exc.stderr if isinstance(exc.stderr, str) else ""
        stderr += f"\ncommand timed out after {timeout_seconds} seconds"
    record = {
        "label": label,
        "command": shlex.join(command),
        "started_at_utc": started,
        "finished_at_utc": utc_now(),
        "returncode": returncode,
        "stdout_tail": stdout[-4000:],
        "stderr_tail": stderr[-4000:],
    }
    return record


def require_success(record: dict[str, Any]) -> None:
    require(
        record["returncode"] == 0,
        f"{record['label']} failed with exit {record['returncode']}: "
        f"{record['stderr_tail'][-1000:] or record['stdout_tail'][-1000:]}",
    )


def compare_identity(
    expected: Any,
    path: Path,
    label: str,
    *,
    require_path: bool = True,
) -> dict[str, Any]:
    require(isinstance(expected, dict), f"{label} identity missing")
    observed = identity(path)
    if require_path:
        require(expected.get("path") == observed["path"], f"{label} path mismatch")
    require(expected.get("size_bytes") == observed["size_bytes"], f"{label} size mismatch")
    require(expected.get("sha256") == observed["sha256"], f"{label} SHA mismatch")
    return observed


def require_within(path: Path, root: Path, label: str) -> None:
    resolved = path.resolve(strict=True)
    root_resolved = root.resolve(strict=True)
    try:
        resolved.relative_to(root_resolved)
    except ValueError as exc:
        raise FinalizeError(f"{label} escapes release staging: {path}") from exc


def verify_sidecar(path: Path) -> None:
    sidecar = path.with_name(f"{path.name}.sha256")
    require(sidecar.is_file() and not sidecar.is_symlink(), f"missing sidecar: {sidecar}")
    parts = sidecar.read_text(encoding="ascii").strip().split()
    require(len(parts) == 2, f"malformed sidecar: {sidecar}")
    require(parts[0] == sha256_path(path), f"sidecar SHA mismatch: {path}")
    require(parts[1] == path.name, f"sidecar filename mismatch: {path}")


def validate_authority_manifest(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    require(manifest.get("schema_name") == AUTHORITY_SCHEMA, "authority schema mismatch")
    require(manifest.get("schema_version") == "1.0.0", "authority version mismatch")
    artifacts = manifest.get("artifacts")
    require(isinstance(artifacts, list), "authority artifacts missing")
    indexed: dict[str, dict[str, Any]] = {}
    for artifact in artifacts:
        require(isinstance(artifact, dict), "authority artifact is not an object")
        artifact_id = artifact.get("artifact_id")
        require(isinstance(artifact_id, str) and artifact_id, "artifact ID missing")
        require(artifact_id not in indexed, f"duplicate authority artifact: {artifact_id}")
        indexed[artifact_id] = artifact
    require(set(indexed) == EXPECTED_ARTIFACT_IDS, "authority artifact set mismatch")
    return indexed


def validate_report_data(payload: dict[str, Any]) -> None:
    require(payload.get("schema_name") == REPORT_DATA_SCHEMA, "report-data schema mismatch")
    require(payload.get("schema_version") == REPORT_DATA_VERSION, "report-data version mismatch")
    require(
        payload.get("report_status") == "validated-derived-observation",
        "report-data status mismatch",
    )
    require(
        payload.get("cn_loh_readiness", {}).get("status") == "NOT_INTEGRATED",
        "CN/LOH status mismatch",
    )
    require(
        payload.get("methyl_auxiliary", {}).get("status") == "association-only",
        "methyl status mismatch",
    )
    require(
        payload.get("tree_decision", {}).get("materialization_status")
        == "POLICY_ONLY_NOT_MATERIALIZED",
        "tree-decision materialization mismatch",
    )
    checks = payload.get("conservation_checks")
    require(
        isinstance(checks, list)
        and checks
        and all(item.get("status") == "PASS" for item in checks),
        "report-data conservation checks failed",
    )
    sample_ids = [row.get("sample") for row in payload.get("per_sample", [])]
    require(
        sample_ids
        == [
            "HCC1395",
            "HCC1395_DORADO",
            "COLO829",
            "H1437",
            "H2009",
            "HCC1937",
            "HCC1954",
        ],
        "report-data sample scope mismatch",
    )


def validate_build_receipt(
    payload: dict[str, Any],
    manifest: Path,
    registry: Path,
    report_data: Path,
    report_html: Path,
    build_receipt: Path,
    staging: Path,
    authority_artifacts: dict[str, dict[str, Any]],
) -> None:
    require(payload.get("schema_name") == BUILD_SCHEMA, "build receipt schema mismatch")
    require(payload.get("schema_version") == BUILD_VERSION, "build receipt version mismatch")
    require(payload.get("all_pass") is True, "build receipt is not all-pass")
    inputs = payload.get("inputs")
    outputs = payload.get("outputs")
    checks = payload.get("checks")
    require(isinstance(inputs, dict), "build receipt inputs missing")
    require(isinstance(outputs, dict), "build receipt outputs missing")
    require(isinstance(checks, dict), "build receipt checks missing")
    compare_identity(inputs.get("authority_manifest"), manifest, "build manifest")
    compare_identity(inputs.get("denominator_registry"), registry, "build registry")
    verified = inputs.get("verified_authority_artifacts")
    require(isinstance(verified, list) and len(verified) == 13, "authority proof count drift")
    verified_by_id = {
        item.get("artifact_id"): item for item in verified if isinstance(item, dict)
    }
    require(set(verified_by_id) == EXPECTED_ARTIFACT_IDS, "authority proofs incomplete")
    for artifact_id, record in authority_artifacts.items():
        source = Path(str(record.get("path", "")))
        observed = compare_identity(
            verified_by_id[artifact_id],
            source,
            f"authority proof {artifact_id}",
        )
        require(observed["sha256"] == record.get("sha256"), f"authority drift: {artifact_id}")
    nested = inputs.get("strict_nested_bundle_identities")
    require(isinstance(nested, list) and len(nested) == 9, "nested strict proofs drift")
    require(len({item.get("bundle_name") for item in nested}) == 9, "nested proof duplicate")
    for item in nested:
        compare_identity(item, Path(item["path"]), f"strict nested {item.get('bundle_name')}")
    compare_identity(outputs.get("report_data"), report_data, "build report-data")
    compare_identity(outputs.get("report_html"), report_html, "build HTML")
    compare_identity(
        outputs.get("report_data_sidecar"),
        report_data.with_name(f"{report_data.name}.sha256"),
        "build report-data sidecar",
    )
    compare_identity(
        outputs.get("report_html_sidecar"),
        report_html.with_name(f"{report_html.name}.sha256"),
        "build HTML sidecar",
    )
    verify_sidecar(report_data)
    verify_sidecar(report_html)
    verify_sidecar(build_receipt)
    require(checks.get("authority_hash_count") == 13, "authority check count drift")
    require(checks.get("strict_nested_bundle_hash_count") == 9, "nested check count drift")
    require(checks.get("denominator_row_count") == 19, "denominator check count drift")
    conservation = checks.get("conservation")
    render = checks.get("render")
    require(
        isinstance(conservation, list)
        and conservation
        and all(item.get("status") == "PASS" for item in conservation),
        "build conservation checks failed",
    )
    require(
        isinstance(render, list)
        and render
        and all(item.get("status") == "PASS" for item in render),
        "build render checks failed",
    )
    require(checks.get("canonical_data_unmodified") is True, "canonical status drift")
    require(checks.get("cn_loh_status") == "NOT_INTEGRATED", "build CN status drift")
    require(checks.get("methyl_status") == "association-only", "build methyl status drift")
    for output_path in (
        report_data,
        report_html,
        build_receipt,
        report_data.with_name(f"{report_data.name}.sha256"),
        report_html.with_name(f"{report_html.name}.sha256"),
        build_receipt.with_name(f"{build_receipt.name}.sha256"),
    ):
        require_within(output_path, staging, "builder output")


def validate_qa_receipt(
    payload: dict[str, Any],
    report_data: Path,
    report_html: Path,
    qa_receipt: Path,
    qa_dir: Path,
) -> None:
    require(payload.get("schema_name") == QA_SCHEMA, "QA receipt schema mismatch")
    require(payload.get("schema_version") == QA_VERSION, "QA receipt version mismatch")
    require(payload.get("all_pass") is True, "browser QA receipt is not all-pass")
    compare_identity(payload.get("inputs", {}).get("html"), report_html, "QA HTML input")
    compare_identity(
        payload.get("inputs", {}).get("report_data"),
        report_data,
        "QA report-data input",
    )
    checks = payload.get("checks")
    require(isinstance(checks, dict), "QA checks missing")
    require(set(checks) == EXPECTED_QA_CHECKS, "QA check set drift")
    require(all(value is True for value in checks.values()), "QA checks failed")
    require(
        set(payload.get("viewports", {}))
        == {"desktop", "laptop", "mobile", "narrow"},
        "QA viewport set drift",
    )
    require(payload.get("no_js", {}).get("observed"), "no-JS evidence missing")
    require(
        payload.get("print", {}).get("observed", {}).get("pdf_page_count", 0) >= 1,
        "print-page evidence missing",
    )
    require_within(qa_receipt, qa_dir, "QA receipt")


def move_failed_attempt(staging: Path, release_root: Path, release_id: str) -> Path | None:
    if not staging.exists():
        return None
    failed_root = release_root / "failed_attempts"
    failed_root.mkdir(parents=True, exist_ok=True)
    destination = failed_root / f"{release_id}.{timestamp_slug()}.{os.getpid()}"
    require(not destination.exists(), f"failed-attempt destination exists: {destination}")
    staging.rename(destination)
    return destination


def main() -> int:
    args = parse_args()
    require(
        RELEASE_ID_RE.fullmatch(args.release_id) is not None,
        "release-id contains unsafe characters",
    )
    require(args.command_timeout_seconds > 0, "command timeout must be positive")
    manifest = args.authority_manifest.expanduser().resolve()
    registry = args.denominator_registry.expanduser().resolve()
    builder = args.builder.expanduser().resolve()
    qa_script = args.qa_script.expanduser().resolve()
    release_root = args.release_root.expanduser().resolve()
    final_dir = release_root / args.release_id
    staging_root = release_root / ".staging"
    require(manifest.is_file(), f"missing authority manifest: {manifest}")
    require(registry.is_file(), f"missing denominator registry: {registry}")
    require(builder.is_file(), f"missing builder: {builder}")
    require(qa_script.is_file(), f"missing QA script: {qa_script}")
    if not args.allow_nonproduction_tools:
        require(builder == DEFAULT_BUILDER.resolve(), "custom builder is non-production")
        require(qa_script == DEFAULT_QA.resolve(), "custom QA script is non-production")
    manifest_payload = read_json(manifest)
    authority_artifacts = validate_authority_manifest(manifest_payload)
    release_root.mkdir(parents=True, exist_ok=True)
    require(not release_root.is_symlink(), f"release root is a symlink: {release_root}")
    require(not final_dir.exists(), f"release already exists: {final_dir}")
    staging_root.mkdir(parents=True, exist_ok=True)
    lock_root = release_root / ".locks"
    lock_root.mkdir(parents=True, exist_ok=True)
    token = secrets.token_hex(12)
    staging = staging_root / (
        f"{args.release_id}.{timestamp_slug()}.{os.getpid()}.{token}"
    )
    lock_path = lock_root / f"{args.release_id}.lock.json"
    lock_payload = {
        "schema_name": "intersubmod.exact_ps_observation_report.release_reservation",
        "schema_version": "1.0.0",
        "created_at_utc": utc_now(),
        "release_id": args.release_id,
        "token": token,
        "pid": os.getpid(),
        "staging_path": str(staging),
        "final_path": str(final_dir),
    }
    write_json_exclusive(lock_path, lock_payload)
    require(not final_dir.exists(), f"release appeared after reservation: {final_dir}")
    staging.mkdir()

    commands: list[dict[str, Any]] = []
    renamed_final = False
    try:
        build_command = [
            args.python,
            str(builder),
            "--authority-manifest",
            str(manifest),
            "--denominator-registry",
            str(registry),
            "--output-dir",
            str(staging),
        ]
        build_record = run_command(
            build_command,
            "build_report",
            args.command_timeout_seconds,
        )
        commands.append(build_record)
        require_success(build_record)
        report_data = staging / "report_data.json"
        report_html = staging / "report.standalone.html"
        build_receipt = staging / "report_build_receipt.json"
        for path in (
            report_data,
            report_html,
            build_receipt,
            staging / "report_data.json.sha256",
            staging / "report.standalone.html.sha256",
            staging / "report_build_receipt.json.sha256",
        ):
            require(path.is_file(), f"builder omitted required output: {path}")
            require(not path.is_symlink(), f"builder emitted symlink: {path}")
        require(
            {path.name for path in staging.iterdir()} == EXPECTED_BUILD_FILES,
            "builder output file set drift",
        )
        build_payload = read_json(build_receipt)
        report_payload = read_json(report_data)
        validate_report_data(report_payload)
        validate_build_receipt(
            build_payload,
            manifest,
            registry,
            report_data,
            report_html,
            build_receipt,
            staging,
            authority_artifacts,
        )
        post_build_identities = {
            "report_data": identity(report_data),
            "report_html": identity(report_html),
            "build_receipt": identity(build_receipt),
        }

        qa_dir = staging / "qa"
        qa_command = [
            args.python,
            str(qa_script),
            "--html",
            str(report_html),
            "--report-data",
            str(report_data),
            "--output-dir",
            str(qa_dir),
        ]
        if args.chromium is not None:
            qa_command.extend(["--chromium", str(args.chromium.expanduser().resolve())])
        qa_record = run_command(
            qa_command,
            "browser_qa",
            args.command_timeout_seconds,
        )
        commands.append(qa_record)
        require_success(qa_record)
        qa_receipt = qa_dir / "browser_qa_receipt.json"
        qa_payload = read_json(qa_receipt)
        require(not qa_dir.is_symlink(), f"QA directory is a symlink: {qa_dir}")
        require(
            {path.name for path in qa_dir.iterdir()}
            == {
                "desktop.png",
                "laptop.png",
                "mobile.png",
                "narrow.png",
                "no_js.png",
                "report_A4.pdf",
                "browser_qa_receipt.json",
            },
            "QA output file set drift",
        )
        validate_qa_receipt(
            qa_payload,
            report_data,
            report_html,
            qa_receipt,
            qa_dir,
        )
        for label, path in (
            ("report_data", report_data),
            ("report_html", report_html),
            ("build_receipt", build_receipt),
        ):
            compare_identity(
                post_build_identities[label],
                path,
                f"post-QA {label}",
            )
        for path in staging.rglob("*"):
            require(not path.is_symlink(), f"release staging contains symlink: {path}")

        # The directory becomes visible before its final marker.  Consumers must
        # require release_receipt.json(.sha256), so an interrupted rename remains
        # fail closed rather than appearing as a valid release.
        staging.rename(final_dir)
        renamed_final = True
        final_manifest = final_dir / "report_data.json"
        final_html = final_dir / "report.standalone.html"
        final_build_receipt = final_dir / "report_build_receipt.json"
        final_qa_receipt = final_dir / "qa/browser_qa_receipt.json"
        scientific_status = (
            "INCOMPLETE_WITH_ABSTAIN"
            if report_payload["headline"]["resource_abstain"] > 0
            or report_payload["cn_loh_readiness"]["status"] == "NOT_INTEGRATED"
            else "COMPLETE_WITHIN_MODEL"
        )
        release_checks = {
            "authority_manifest_contract_pass": True,
            "build_receipt_schema_and_identity_pass": True,
            "browser_qa_schema_and_identity_pass": True,
            "post_qa_hashes_unchanged": True,
            "symlinks_zero": True,
            "cn_loh_not_inferred": report_payload["cn_loh_readiness"]["status"]
            == "NOT_INTEGRATED",
            "methyl_association_only": report_payload["methyl_auxiliary"]["status"]
            == "association-only",
            "canonical_data_unmodified": build_payload["checks"][
                "canonical_data_unmodified"
            ]
            is True,
            "final_marker_written_last": True,
        }
        receipt = {
            "schema_name": SCHEMA_NAME,
            "schema_version": SCHEMA_VERSION,
            "created_at_utc": utc_now(),
            "release_id": args.release_id,
            "release_root": str(release_root),
            "release_path": str(final_dir),
            "statuses": {
                "data_technical_status": manifest_payload["project"]["project_status"],
                "scientific_completeness_status": scientific_status,
                "report_build_status": "PASS",
                "browser_qa_status": "PASS",
                "release_status": "VALIDATED_DERIVED_OBSERVATION",
            },
            "claim_ceiling": report_payload["claim_boundary"]["canonical_sentence"],
            "inputs": {
                "authority_manifest": identity(manifest),
                "denominator_registry": identity(registry),
                "builder": identity(builder),
                "qa_script": identity(qa_script),
                "release_reservation": identity(lock_path),
            },
            "outputs": {
                "report_data": identity(final_manifest),
                "html": identity(final_html),
                "build_receipt": identity(final_build_receipt),
                "browser_qa_receipt": identity(final_qa_receipt),
            },
            "commands": commands,
            "checks": release_checks,
        }
        receipt["all_pass"] = all(release_checks.values())
        require(receipt["all_pass"], "release checks are not all-pass")
        release_receipt = final_dir / "release_receipt.json"
        write_json_exclusive(release_receipt, receipt)
        sidecar = write_sidecar_exclusive(release_receipt)
        print(f"[INPUT] authority_manifest={manifest}")
        print(f"[INPUT] denominator_registry={registry}")
        print(f"[OUTPUT] release={final_dir}")
        print(f"[OUTPUT] receipt={release_receipt}")
        print(f"[OUTPUT] receipt_sha256={identity(sidecar)['sha256']}")
        print(
            json.dumps(
                {"all_pass": True, "statuses": receipt["statuses"]},
                ensure_ascii=False,
                sort_keys=True,
            )
        )
        return 0
    except Exception as exc:  # noqa: BLE001 - preserve failed staging for audit
        failed: Path | None = None
        try:
            owned_path = final_dir if renamed_final else staging
            if owned_path.exists():
                failure_receipt = owned_path / "failure_receipt.json"
                if not failure_receipt.exists():
                    write_json_exclusive(
                        failure_receipt,
                        {
                            "schema_name": (
                                "intersubmod.exact_ps_observation_report."
                                "failed_attempt"
                            ),
                            "schema_version": "1.0.0",
                            "created_at_utc": utc_now(),
                            "release_id": args.release_id,
                            "reservation_token": token,
                            "failure": f"{type(exc).__name__}: {exc}",
                            "commands": commands,
                        },
                    )
                failed = move_failed_attempt(
                    owned_path,
                    release_root,
                    args.release_id,
                )
            if lock_path.exists():
                if failed is None:
                    failed_root = release_root / "failed_attempts"
                    failed_root.mkdir(parents=True, exist_ok=True)
                    failed = failed_root / (
                        f"{args.release_id}.{timestamp_slug()}.{os.getpid()}"
                    )
                    failed.mkdir()
                lock_path.rename(failed / "release_reservation.json")
        except Exception as move_exc:  # noqa: BLE001
            print(f"[ERROR] failed-attempt preservation failed: {move_exc}", file=sys.stderr)
        print(f"[ERROR] {type(exc).__name__}: {exc}", file=sys.stderr)
        if failed is not None:
            print(f"[OUTPUT] failed_attempt={failed}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
