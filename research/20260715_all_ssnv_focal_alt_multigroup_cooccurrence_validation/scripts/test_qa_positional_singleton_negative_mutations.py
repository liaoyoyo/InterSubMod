#!/usr/bin/env python3
"""Exercise fail-closed standalone-report QA with retained negative mutations."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import stat
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable
from urllib.parse import quote


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qa-script", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--positive-receipt", type=Path, required=True)
    parser.add_argument("--executable-path", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def artifact(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(resolved, flags)
    try:
        status_before = os.fstat(descriptor)
        if not stat.S_ISREG(status_before.st_mode):
            raise RuntimeError(f"Expected regular file: {resolved}")
        digest = hashlib.sha256()
        size_bytes = 0
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
            size_bytes += len(chunk)
        status_after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    stable_fields = (
        "st_dev",
        "st_ino",
        "st_mode",
        "st_size",
        "st_mtime_ns",
        "st_ctime_ns",
    )
    if any(
        getattr(status_before, field) != getattr(status_after, field)
        for field in stable_fields
    ) or size_bytes != status_before.st_size:
        raise RuntimeError(f"File changed while its identity was read: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": size_bytes,
        "sha256": digest.hexdigest(),
        "mode": oct(stat.S_IMODE(status_before.st_mode)),
    }


def write_exclusive(path: Path, content: str) -> None:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    path.chmod(0o444)


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"Expected JSON object: {path}")
    return value


def mutate_external_resource(source: str) -> str:
    return source.replace(
        "</style>",
        "body{background-image:url(file:///etc/passwd)!important}</style>",
        1,
    )


def mutate_fixed_overlay(source: str) -> str:
    overlay = (
        '<div id="qa-invalid-overlay" '
        'style="position:fixed;inset:0;z-index:999999;background:#fff"></div>'
    )
    return source.replace("</body>", overlay + "</body>", 1)


def mutate_pseudo_overlay(source: str) -> str:
    return source.replace(
        "</style>",
        "body::after{content:'';position:relative;display:block;width:100%;"
        "height:10000vh;margin-top:-10000vh;z-index:999999;background:#fff;"
        "pointer-events:none}</style>",
        1,
    )


def mutate_nearblank_image(source: str) -> str:
    mutated = source.replace("fill='%23147d75'", "fill='%23ffffff'", 1)
    return mutated.replace("fill='%23d08a17'", "fill='%23ffffff'", 1)


def mutate_sparse_nearblank_image(source: str) -> str:
    stripes = "".join(
        f"<rect x='{x}' width='1' height='120' fill='#000000'/>"
        for x in (10, 45, 80, 115, 150, 185)
    )
    svg = (
        "<svg xmlns='http://www.w3.org/2000/svg' width='200' height='120'>"
        "<rect width='200' height='120' fill='#ffffff'/>"
        f"{stripes}</svg>"
    )
    replacement = "data:image/svg+xml," + quote(svg, safe="'")
    start = source.index("data:image/svg+xml,")
    end = source.index('"', start)
    return source[:start] + replacement + source[end:]


def mutate_masked_image(source: str) -> str:
    return source.replace(
        "</style>",
        "img{mask-image:linear-gradient(transparent,transparent)!important;"
        "-webkit-mask-image:linear-gradient(transparent,transparent)!important}"
        "</style>",
        1,
    )


def mutate_delayed_inline_handler(source: str) -> str:
    return source.replace(
        "<body>",
        '<body onload="setTimeout(()=>new WebSocket(\'wss://example.invalid/qa\'),600000)">',
        1,
    )


def mutate_duplicate_anchor(source: str) -> str:
    return source.replace('href="#s8"', 'href="#s7"', 1)


def mutate_print_overflow(source: str) -> str:
    return source.replace(
        "</style>",
        "@media print{.table-wrap{width:2000px!important;max-width:none!important}}"
        "</style>",
        1,
    )


def mutate_nearblank_pdf_page(source: str) -> str:
    style = (
        "#qa-nearblank-print-page{display:none}"
        "@media print{#qa-nearblank-print-page{display:block!important;"
        "break-before:page;page-break-before:always;width:70px;height:70px;"
        "margin:20px;background:#000;color:#000;font-size:1px;line-height:1px}}"
    )
    element = (
        '<div id="qa-nearblank-print-page" aria-hidden="true">'
        "ABCDEFGHIJKLMNOPQRST</div>"
    )
    return source.replace("</style>", style + "</style>", 1).replace(
        "</body>", element + "</body>", 1
    )


def layouts(payload: dict[str, Any]) -> list[dict[str, Any]]:
    observed: list[dict[str, Any]] = []
    for profile in payload.get("profiles", []):
        for key in ("diagnostics", "open_diagnostics", "post_interaction_diagnostics"):
            value = profile.get(key)
            if isinstance(value, dict):
                observed.append(value)
    return observed


def signal_external_resource(payload: dict[str, Any]) -> bool:
    return any(layout.get("resourceViolations") for layout in layouts(payload)) or any(
        profile.get("network_attempts") for profile in payload.get("profiles", [])
    )


def signal_fixed_overlay(payload: dict[str, Any]) -> bool:
    return any(layout.get("unexpectedPositionedElements") for layout in layouts(payload))


def signal_pseudo_overlay(payload: dict[str, Any]) -> bool:
    return any(
        layout.get("unexpectedPositionedPseudoElements") for layout in layouts(payload)
    )


def signal_nearblank_image(payload: dict[str, Any]) -> bool:
    return any(
        image.get("pass") is False
        for layout in layouts(payload)
        for image in layout.get("images", [])
    )


def signal_masked_image(payload: dict[str, Any]) -> bool:
    return signal_nearblank_image(payload)


def signal_inline_handler(payload: dict[str, Any]) -> bool:
    security = payload.get("static_html_security", {})
    return any(
        violation.get("reason") == "inline_event_handler"
        for violation in security.get("violations", [])
    )


def signal_duplicate_anchor(payload: dict[str, Any]) -> bool:
    return any(layout.get("counts") != layout.get("expected") for layout in layouts(payload))


def signal_print_overflow(payload: dict[str, Any]) -> bool:
    for profile in payload.get("profiles", []):
        if profile.get("media") != "print":
            continue
        diagnostic = profile.get("diagnostics", {})
        viewport = diagnostic.get("viewport", {})
        if (
            viewport.get("documentScrollWidth", 0)
            > viewport.get("documentClientWidth", 0) + 1
            or diagnostic.get("outOfBounds")
            or any(table.get("pass") is False for table in diagnostic.get("tables", []))
        ):
            return True
    return False


def signal_nearblank_pdf_page(payload: dict[str, Any]) -> bool:
    for profile in payload.get("profiles", []):
        verification = profile.get("pdf_verification")
        if not isinstance(verification, dict):
            continue
        if any(
            page.get("spatial_coverage_pass") is False
            or page.get("content_pass") is False
            for page in verification.get("raster_pages", [])
        ):
            return True
    return False


def run_mutation(
    *,
    name: str,
    mutation: Callable[[str], str],
    signal: Callable[[dict[str, Any]], bool],
    fixture_text: str,
    input_dir: Path,
    result_dir: Path,
    qa_script: Path,
    executable: Path,
) -> dict[str, Any]:
    input_path = input_dir / f"{name}.html"
    write_exclusive(input_path, mutation(fixture_text))
    expected_output = result_dir / name
    command = [
        sys.executable,
        str(qa_script),
        "--html",
        str(input_path),
        "--executable-path",
        str(executable),
        "--output-dir",
        str(expected_output),
        "--evidence-kind",
        "self_test",
        "--wait-ms",
        "50",
    ]
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    failure_dirs = sorted(result_dir.glob(f"{name}.failed.*"))
    if completed.returncode == 0:
        raise RuntimeError(f"Mutation unexpectedly passed: {name}")
    if expected_output.exists():
        raise RuntimeError(f"Mutation retained an authoritative output: {name}")
    if len(failure_dirs) != 1:
        raise RuntimeError(f"Mutation did not publish exactly one failure directory: {name}")
    failure_dir = failure_dirs[0]
    failure_receipt = failure_dir / "failure_receipt.json"
    if not failure_receipt.is_file():
        raise RuntimeError(f"Mutation lacks failure receipt: {name}")
    if (failure_dir / "qa_receipt.json").exists() or (failure_dir / "_SUCCESS.json").exists():
        raise RuntimeError(f"Mutation retained success artifacts: {name}")
    payload = load_json(failure_receipt)
    if payload.get("pass") is not False or not signal(payload):
        raise RuntimeError(f"Mutation lacks its expected rejection signal: {name}")
    after_identity_checks = {
        "input_html": payload.get("input_html_after") == payload.get("input_html"),
        "executed_html": payload.get("executed_html_after")
        == payload.get("executed_html"),
        "qa_script": payload.get("qa_script_after") == payload.get("qa_script"),
        "browser_executable": payload.get("browser_executable_after")
        == payload.get("browser", {}).get("executable"),
        "runtime_bundle": payload.get("runtime_identity_after")
        == payload.get("runtime_identity"),
        "formal_release_binding": payload.get("formal_release_binding_after")
        == payload.get("formal_release_binding"),
        "capture_errors_absent": payload.get("after_identity_capture_errors") == {},
    }
    if not all(after_identity_checks.values()):
        raise RuntimeError(
            f"Mutation failure receipt lacks complete stable after identity: "
            f"{name}: {after_identity_checks}"
        )
    return {
        "name": name,
        "input": artifact(input_path),
        "command": command,
        "returncode": completed.returncode,
        "stdout_sha256": hashlib.sha256(completed.stdout.encode()).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr.encode()).hexdigest(),
        "failure_directory": str(failure_dir.resolve()),
        "failure_receipt": artifact(failure_receipt),
        "expected_rejection_signal_observed": True,
        "after_identity_checks": after_identity_checks,
        "all_after_identities_bound_and_stable": True,
        "authoritative_output_absent": True,
        "success_artifacts_absent": True,
        "pass": True,
    }


def run_post_publication_failure_quarantine_test(
    *, qa_script: Path, result_dir: Path
) -> dict[str, Any]:
    spec = importlib.util.spec_from_file_location(
        "qa_positional_singleton_standalone_under_test", qa_script
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load QA module for publication-failure test")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    published_output = result_dir / "synthetic_post_publication_failure"
    published_output.mkdir(mode=0o700)
    write_exclusive(published_output / "qa_receipt.json", '{"pass":true}\n')
    write_exclusive(published_output / "_SUCCESS.json", '{"pass":true}\n')
    failure_dir = module.publish_failure(
        published_output,
        published_output,
        {
            "schema_name": "intersubmod.synthetic_post_publication_failure",
            "profiles": [],
            "pass": True,
        },
    )
    quarantine_dirs = sorted(
        result_dir.glob(
            f"{failure_dir.name}.quarantined_success_artifacts"
        )
    )
    if published_output.exists() or len(quarantine_dirs) != 1:
        raise RuntimeError("Post-publication failure did not quarantine atomically")
    quarantine_dir = quarantine_dirs[0]
    failure_receipt = failure_dir / "failure_receipt.json"
    if (
        not failure_receipt.is_file()
        or (failure_dir / "qa_receipt.json").exists()
        or (failure_dir / "_SUCCESS.json").exists()
        or not (quarantine_dir / "qa_receipt.json").is_file()
        or not (quarantine_dir / "_SUCCESS.json").is_file()
    ):
        raise RuntimeError("Post-publication failure quarantine contents are invalid")
    failure_payload = load_json(failure_receipt)
    publication = failure_payload.get("failure_publication", {})
    if (
        failure_payload.get("pass") is not False
        or publication.get("success_markers_present_in_failure_directory") is not False
        or len(publication.get("quarantined_success_artifacts", [])) != 2
    ):
        raise RuntimeError("Post-publication failure receipt contract drift")
    return {
        "failure_directory": str(failure_dir.resolve()),
        "failure_receipt": artifact(failure_receipt),
        "quarantine_directory": str(quarantine_dir.resolve()),
        "quarantined_qa_receipt": artifact(quarantine_dir / "qa_receipt.json"),
        "quarantined_success_marker": artifact(quarantine_dir / "_SUCCESS.json"),
        "authoritative_output_absent": True,
        "success_artifacts_absent_from_failure_directory": True,
        "success_artifacts_preserved_without_overwrite": True,
        "pass": True,
    }


def run_stdout_failure_quarantine_integration(
    *,
    qa_script: Path,
    fixture: Path,
    executable: Path,
    result_dir: Path,
) -> dict[str, Any]:
    name = "stdout_write_failure_integration"
    expected_output = result_dir / name
    command = [
        sys.executable,
        "-u",
        str(qa_script),
        "--html",
        str(fixture),
        "--executable-path",
        str(executable),
        "--output-dir",
        str(expected_output),
        "--evidence-kind",
        "self_test",
        "--wait-ms",
        "50",
    ]
    environment = os.environ.copy()
    environment["PYTHONUNBUFFERED"] = "1"
    with Path("/dev/full").open("wb", buffering=0) as stdout_target:
        completed = subprocess.run(
            command,
            check=False,
            stdout=stdout_target,
            stderr=subprocess.PIPE,
            text=True,
            env=environment,
        )
    failure_dirs = [
        path
        for path in sorted(result_dir.glob(f"{name}.failed.*"))
        if ".quarantined_success_artifacts" not in path.name
    ]
    quarantine_dirs = sorted(
        result_dir.glob(f"{name}.failed.*.quarantined_success_artifacts")
    )
    if completed.returncode == 0 or expected_output.exists():
        raise RuntimeError("Real stdout failure retained an authoritative QA output")
    if len(failure_dirs) != 1 or len(quarantine_dirs) != 1:
        raise RuntimeError("Real stdout failure did not publish one failure and quarantine")
    failure_dir = failure_dirs[0]
    quarantine_dir = quarantine_dirs[0]
    failure_receipt = failure_dir / "failure_receipt.json"
    if (
        not failure_receipt.is_file()
        or (failure_dir / "qa_receipt.json").exists()
        or (failure_dir / "_SUCCESS.json").exists()
        or not (quarantine_dir / "qa_receipt.json").is_file()
        or not (quarantine_dir / "_SUCCESS.json").is_file()
    ):
        raise RuntimeError("Real stdout failure quarantine contents are invalid")
    failure_payload = load_json(failure_receipt)
    publication = failure_payload.get("failure_publication", {})
    if (
        failure_payload.get("pass") is not False
        or failure_payload.get("exception", {}).get("type") not in {"OSError", "BlockingIOError"}
        or len(publication.get("quarantined_success_artifacts", [])) != 2
    ):
        raise RuntimeError("Real stdout failure receipt lacks the integrated failure signal")
    return {
        "command": command,
        "stdout_target": "/dev/full",
        "returncode": completed.returncode,
        "stderr_sha256": hashlib.sha256(completed.stderr.encode()).hexdigest(),
        "failure_directory": str(failure_dir.resolve()),
        "failure_receipt": artifact(failure_receipt),
        "quarantine_directory": str(quarantine_dir.resolve()),
        "quarantined_qa_receipt": artifact(quarantine_dir / "qa_receipt.json"),
        "quarantined_success_marker": artifact(quarantine_dir / "_SUCCESS.json"),
        "authoritative_output_absent": True,
        "success_artifacts_absent_from_failure_directory": True,
        "success_artifacts_preserved_without_overwrite": True,
        "pass": True,
    }


def main() -> None:
    args = parse_args()
    harness_before = artifact(Path(__file__).resolve())
    qa_script = args.qa_script.resolve(strict=True)
    fixture = args.fixture.resolve(strict=True)
    positive_receipt_path = args.positive_receipt.resolve(strict=True)
    executable = args.executable_path.resolve(strict=True)
    output_root = args.output_root.resolve(strict=False)
    qa_script_identity = artifact(qa_script)
    fixture_identity = artifact(fixture)
    executable_identity = artifact(executable)
    positive_receipt_identity = artifact(positive_receipt_path)
    if output_root.exists():
        raise SystemExit(f"Refusing to overwrite output root: {output_root}")

    positive_receipt = load_json(positive_receipt_path)
    if (
        positive_receipt.get("schema_name")
        != "intersubmod.positional_singleton_standalone_qa"
        or positive_receipt.get("schema_version") != "2.4.0"
        or positive_receipt.get("evidence_kind") != "self_test"
        or positive_receipt.get("validation_eligibility")
        != "self_test_only_not_formal_report_evidence"
        or positive_receipt.get("pass") is not True
        or positive_receipt.get("runtime_identity_stable") is not True
        or positive_receipt.get("runtime_identity")
        != positive_receipt.get("runtime_identity_after")
        or positive_receipt.get("static_html_security", {}).get("pass") is not True
        or positive_receipt.get("qa_script") != qa_script_identity
        or positive_receipt.get("qa_script_after") != qa_script_identity
        or positive_receipt.get("input_html") != fixture_identity
        or positive_receipt.get("input_html_after") != fixture_identity
        or positive_receipt.get("browser", {}).get("executable") != executable_identity
        or positive_receipt.get("browser_executable_after") != executable_identity
        or positive_receipt.get("executed_html", {}).get("size_bytes")
        != fixture_identity["size_bytes"]
        or positive_receipt.get("executed_html", {}).get("sha256")
        != fixture_identity["sha256"]
    ):
        raise RuntimeError("Positive self-test receipt contract drift")
    positive_images = [
        image
        for layout in layouts(positive_receipt)
        for image in layout.get("images", [])
    ]
    if not any(
        image.get("pass") is True
        and image.get("pixelContentMode") == "sparse_spatially_distributed"
        and image.get("contentRatioRelativeError", 1) <= 0.005
        and image.get("boxRatioRelativeError", 0) > 0.01
        for image in positive_images
    ):
        raise RuntimeError(
            "Positive self-test lacks a passing sparse wide-image border canary"
        )

    output_root.mkdir(mode=0o700, parents=True)
    input_dir = output_root / "mutated_inputs"
    result_dir = output_root / "results"
    input_dir.mkdir(mode=0o700)
    result_dir.mkdir(mode=0o700)

    fixture_text = fixture.read_text(encoding="utf-8")
    cases = (
        ("external_file_css_url", mutate_external_resource, signal_external_resource),
        ("fixed_full_page_overlay", mutate_fixed_overlay, signal_fixed_overlay),
        ("pseudo_element_overlay", mutate_pseudo_overlay, signal_pseudo_overlay),
        ("nearblank_image", mutate_nearblank_image, signal_nearblank_image),
        (
            "sparse_nearblank_image",
            mutate_sparse_nearblank_image,
            signal_nearblank_image,
        ),
        ("masked_image", mutate_masked_image, signal_masked_image),
        ("delayed_inline_handler", mutate_delayed_inline_handler, signal_inline_handler),
        ("duplicate_anchor_target", mutate_duplicate_anchor, signal_duplicate_anchor),
        ("print_horizontal_overflow", mutate_print_overflow, signal_print_overflow),
        ("nearblank_pdf_page", mutate_nearblank_pdf_page, signal_nearblank_pdf_page),
    )
    outcomes = [
        run_mutation(
            name=name,
            mutation=mutation,
            signal=signal,
            fixture_text=fixture_text,
            input_dir=input_dir,
            result_dir=result_dir,
            qa_script=qa_script,
            executable=executable,
        )
        for name, mutation, signal in cases
    ]
    post_publication_failure = run_post_publication_failure_quarantine_test(
        qa_script=qa_script,
        result_dir=result_dir,
    )
    stdout_failure_integration = run_stdout_failure_quarantine_integration(
        qa_script=qa_script,
        fixture=fixture,
        executable=executable,
        result_dir=result_dir,
    )

    formal_output = result_dir / "formal_mode_without_signed_release"
    formal_command = [
        sys.executable,
        str(qa_script),
        "--html",
        str(fixture),
        "--executable-path",
        str(executable),
        "--output-dir",
        str(formal_output),
        "--evidence-kind",
        "formal_report",
        "--wait-ms",
        "50",
    ]
    formal_result = subprocess.run(
        formal_command, check=False, capture_output=True, text=True
    )
    if (
        formal_result.returncode == 0
        or formal_output.exists()
        or "requires signed release arguments" not in formal_result.stderr
    ):
        raise RuntimeError("Formal mode accepted an unsigned release context")

    harness_after = artifact(Path(__file__).resolve())
    if harness_after != harness_before:
        raise RuntimeError("Negative mutation harness identity changed during execution")
    receipt = {
        "schema_name": "intersubmod.positional_singleton_qa_negative_mutations",
        "schema_version": "1.6.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_type": "F_demo_illustration",
        "validation_eligibility": "QA_harness_self_test_only",
        "command": [str(Path(sys.executable).resolve()), *sys.argv],
        "inputs": {
            "negative_mutation_harness": harness_before,
            "negative_mutation_harness_after": harness_after,
            "qa_script": qa_script_identity,
            "fixture": fixture_identity,
            "positive_self_test_receipt": positive_receipt_identity,
            "chromium": executable_identity,
            "python": artifact(Path(sys.executable)),
        },
        "positive_receipt_transitive_binding": {
            "qa_script_bound": True,
            "fixture_bound": True,
            "chromium_bound": True,
            "runtime_before_after_equal": True,
            "executed_html_bytes_bound": True,
            "sparse_wide_image_border_canary_passed": True,
            "pass": True,
        },
        "verification_trust_model": {
            "environment": "trusted_local_scientific_QA_without_malicious_same_uid_concurrent_mutation",
            "source_identity_scope": (
                "single_open_file_descriptor_identity_before_after_stability_check;_"
                "not_self_attested_kernel_executed_byte_provenance"
            ),
            "independent_source_review_required": True,
        },
        "mutations": outcomes,
        "post_publication_failure_quarantine": post_publication_failure,
        "stdout_failure_quarantine_integration": stdout_failure_integration,
        "formal_unsigned_rejection": {
            "command": formal_command,
            "returncode": formal_result.returncode,
            "stderr_sha256": hashlib.sha256(formal_result.stderr.encode()).hexdigest(),
            "required_error_observed": True,
            "output_absent": True,
            "pass": True,
        },
        "counts": {
            "negative_mutations": len(outcomes),
            "negative_mutations_rejected": sum(outcome["pass"] for outcome in outcomes),
            "negative_mutations_with_complete_after_identity": sum(
                outcome["all_after_identities_bound_and_stable"]
                for outcome in outcomes
            ),
            "unsigned_formal_contexts_rejected": 1,
            "post_publication_failures_quarantined": 2,
        },
        "pass": all(outcome["pass"] for outcome in outcomes)
        and post_publication_failure["pass"]
        and stdout_failure_integration["pass"],
    }
    receipt_path = output_root / "negative_mutation_receipt.json"
    write_exclusive(
        receipt_path,
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )
    input_dir.chmod(0o555)
    result_dir.chmod(0o555)
    output_root.chmod(0o555)
    print(
        json.dumps(
            {
                "output_root": str(output_root),
                "receipt": str(receipt_path),
                "negative_mutations_rejected": len(outcomes),
                "unsigned_formal_contexts_rejected": 1,
                "post_publication_failures_quarantined": 2,
                "pass": True,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
