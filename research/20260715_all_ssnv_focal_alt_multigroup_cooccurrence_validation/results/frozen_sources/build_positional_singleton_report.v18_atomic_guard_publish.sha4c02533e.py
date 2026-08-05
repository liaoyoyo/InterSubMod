#!/usr/bin/env python3
"""Build the positional-singleton methyl multigroup Markdown and HTML report."""

from __future__ import annotations

import argparse
import base64
from collections import Counter
import ctypes
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import html
import json
import math
import os
from pathlib import Path
import stat
import subprocess
import sys
from typing import Any, Iterable, Mapping
import uuid
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
TOPIC_ROOT = SCRIPT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import finalize_task_b_result_release as FINAL_RELEASE  # noqa: E402


SCHEMA_NAME = "intersubmod.positional_singleton_methyl_multigroup_report"
SCHEMA_VERSION = "1.0.0"
AUDIT_SCHEMA = "intersubmod.positional_singleton_methyl_multigroup_audit"
EXPECTED_DATASETS = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)
EXPECTED_POSITIVE_DATASETS = tuple(
    dataset for dataset in EXPECTED_DATASETS if dataset != "COLO829"
)
REQUIRED_VISUAL_CASES = (
    ("HCC1395_DORADO", "chr1", 20_467_811, "G", "C"),
    ("HCC1395_DORADO", "chr22", 47_518_662, "A", "G"),
    ("HCC1937", "chr5", 43_849_776, "T", "C"),
)
EXPECTED_COUNTS = {
    "all_dataset_sites": 469_849,
    "singleton_sites": 50_432,
    "m1_evaluable": 48_347,
    "m1_flagged": 5_961,
    "m2_pass": 30,
    "m2_fail": 18,
    "m2_not_evaluable": 5_913,
    "m2_not_run": 44_471,
}
EXPECTED_PHASE_SET_AUDIT = {
    "site_count": 30,
    "total_core_reads": 3_505,
    "total_missing_ps": 131,
    "sites_with_missing_ps": 10,
    "sites_with_multiple_nonmissing_ps": 1,
    "sites_with_missing_and_multiple_nonmissing_ps": 1,
}
EXPECTED_PHASE_SET_CASE = {
    "key": ("HCC1937", "chr5", 43_849_776, "T", "C"),
    "core_read_n": 109,
    "ps_missing_n": 88,
    "ps_nonmissing_n": 21,
    "ps_distinct_nonmissing_n": 2,
    "ps_values": [43_668_888, 43_913_176],
}
MINIMUM_CANONICAL_TEST_COUNT = 488
REQUIRED_SUPPLEMENTAL_TESTS = frozenset(
    {
        "test_select_cases_includes_dataset_maxima_and_required_caveat_cases",
        "test_summarize_phase_sets_counts_missing_and_multiple_ps",
        "test_summarize_phase_sets_rejects_core_count_drift",
        "test_summarize_phase_sets_rejects_missing_latest_ps",
        "test_validate_phase_set_audit_rejects_expected_count_drift",
        "test_validate_phase_set_audit_accepts_expected_hcc1937_case",
        "test_validate_phase_set_audit_rejects_hcc1937_case_drift",
        "test_validate_formal_release_chain_uses_canonical_fd_validators",
        "test_validate_formal_release_chain_rejects_noncanonical_path",
        "test_validate_formal_release_chain_rejects_report_contract_drift",
        "test_validate_pre_publish_state_detects_late_input_mutation",
        "test_validate_pre_publish_state_detects_transitive_authority_mutation",
        "test_main_orders_publish_guard_after_outputs_and_before_rename",
        "test_validate_and_publish_contains_only_guard_then_atomic_rename",
        "test_validate_canonical_pytest_release_evidence_accepts_signed_current",
        "test_validate_canonical_pytest_release_evidence_rejects_noncanonical_xml",
        "test_validate_canonical_pytest_release_evidence_rejects_report_source_drift",
        "test_validate_canonical_pytest_release_evidence_rejects_unretired_private_key",
    }
)
COLORS = {
    "ink": "#17202a",
    "muted": "#59636e",
    "grid": "#d9dee3",
    "teal": "#147d75",
    "blue": "#2f6ea5",
    "amber": "#d08a17",
    "red": "#b84a3a",
    "gray": "#9aa3ab",
    "light": "#eef2f3",
    "missing": "#e7e9ec",
}
GROUP_COLORS = (
    "#147d75",
    "#d08a17",
    "#2f6ea5",
    "#b84a3a",
    "#7a5aa6",
    "#568c3b",
    "#a2632f",
    "#3f8792",
    "#8f6679",
    "#68727c",
)
OPENSSL = Path("/usr/bin/openssl")
FINAL_DATASET_PUBLIC_KEY_SHA256 = (
    "54598225ab57a52393fbe63a29b24a19f39998bdf5d951fb61f4edab67bfeb24"
)
FINAL_REPORT_PUBLIC_KEY_SHA256 = (
    "98e7ca01a67ce73ac3eea0a18599db78233fd66104f1922db773396ef85f56fb"
)
SOURCE_AUTHORITY_PUBLIC_KEY_SHA256 = (
    "b71855f5fd9d0e97df0f6186420b5cec95f85d8b462fde0a890443846271bee4"
)
SOURCE_AUTHORITY_ID = "20260718_all_ssnv_focal_alt_task_b_release_v4"
CANONICAL_PYTEST_XML = (
    TOPIC_ROOT
    / "logs"
    / "pytest_full_pre_supplemental_v18_atomic_guard_publish.xml"
)
TEST_EVIDENCE_MANIFEST = (
    TOPIC_ROOT / "logs" / "supplemental_report_test_evidence.v2.json"
)
TEST_EVIDENCE_SIGNATURE = Path(str(TEST_EVIDENCE_MANIFEST) + ".ed25519.sig")
TEST_EVIDENCE_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/"
    "intersubmod_supplemental_test_authority/"
    "20260718_positional_singleton_report_v2/ed25519_public.pem"
)
TEST_EVIDENCE_PRIVATE_KEY = TEST_EVIDENCE_PUBLIC_KEY.with_name(
    "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)
TEST_EVIDENCE_PUBLIC_KEY_SHA256 = (
    "a7d7933deff31339f7f7af92ac1cf4db280c12b49b1c36b5db1a955e30220b3a"
)
TEST_EVIDENCE_AUTHORITY_ID = (
    "20260718_positional_singleton_supplemental_test_evidence_v2"
)
TEST_EVIDENCE_REPORT_BUILDER_PATH = Path(__file__).resolve()
TEST_EVIDENCE_TEST_SOURCE_PATH = (
    TOPIC_ROOT / "tests" / "test_build_positional_singleton_report.py"
)

SiteKey = tuple[str, str, int, str, str]


class ReportBuildError(RuntimeError):
    """Raised when a report input or output contract is violated."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path, *, public_path: Path | None = None) -> dict[str, Any]:
    resolved = path.expanduser().resolve(strict=True)
    info = resolved.stat()
    return {
        "path": str(
            public_path.expanduser().resolve()
            if public_path is not None
            else resolved
        ),
        "size_bytes": info.st_size,
        "mtime_ns": info.st_mtime_ns,
        "sha256": sha256(resolved),
        "mode": oct(stat.S_IMODE(info.st_mode)),
    }


def require_identity(path: Path, record: Any, *, label: str) -> None:
    if not isinstance(record, dict):
        raise ReportBuildError(f"{label} identity is missing")
    observed = artifact(path)
    for field in ("path", "size_bytes", "sha256"):
        if observed[field] != record.get(field):
            raise ReportBuildError(
                f"{label} identity drift for {field}: "
                f"{observed[field]!r} != {record.get(field)!r}"
            )
    if "mode" in record and observed["mode"] != record["mode"]:
        raise ReportBuildError(
            f"{label} identity drift for mode: "
            f"{observed['mode']!r} != {record['mode']!r}"
        )


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ReportBuildError(f"Unable to read {label}: {path}") from error
    if not isinstance(value, dict):
        raise ReportBuildError(f"{label} must be a JSON object")
    return value


def validate_pytest_root(root: ET.Element) -> dict[str, int]:
    """Validate canonical suite totals and the current named regression set."""

    suites = [root] if root.tag == "testsuite" else list(root.findall("testsuite"))
    if not suites:
        raise ReportBuildError("Canonical pytest XML contains no test suites")
    counts = {
        field: sum(int(suite.attrib.get(field, "0")) for suite in suites)
        for field in ("tests", "failures", "errors", "skipped")
    }
    test_names = [case.attrib.get("name", "") for case in root.iter("testcase")]
    name_counts = Counter(test_names)
    missing_required = sorted(
        name for name in REQUIRED_SUPPLEMENTAL_TESTS if name_counts[name] != 1
    )
    if (
        counts["tests"] < MINIMUM_CANONICAL_TEST_COUNT
        or counts["failures"] != 0
        or counts["errors"] != 0
        or counts["skipped"] != 0
        or len(test_names) != counts["tests"]
        or missing_required
    ):
        raise ReportBuildError(
            "Canonical pytest suite did not pass current supplemental contract: "
            f"counts={counts}, testcase_count={len(test_names)}, "
            f"missing_or_duplicated_required={missing_required}"
        )
    counts["testcases"] = len(test_names)
    counts["required_supplemental_tests"] = len(REQUIRED_SUPPLEMENTAL_TESTS)
    return counts


def validate_pytest_xml(path: Path) -> dict[str, int]:
    try:
        root = ET.parse(path).getroot()
    except (ET.ParseError, OSError) as error:
        raise ReportBuildError(f"Invalid canonical pytest XML: {path}") from error
    return validate_pytest_root(root)


def validate_canonical_pytest_release_evidence(
    path: Path,
) -> dict[str, Any]:
    """Verify signed test evidence bound to the live report source and tests."""

    if path.resolve() != CANONICAL_PYTEST_XML.resolve():
        raise ReportBuildError(
            "Canonical pytest XML path is not the signed release path"
        )
    try:
        private_mode = stat.S_IMODE(TEST_EVIDENCE_PRIVATE_KEY.stat().st_mode)
    except OSError as error:
        raise ReportBuildError(
            "Supplemental test-evidence private-key state is unavailable"
        ) from error
    if private_mode != 0:
        raise ReportBuildError(
            "Supplemental test-evidence private key is not retired"
        )
    authority = FINAL_RELEASE.SOURCE_AUTHORITY
    try:
        with authority.BoundArtifactReader() as reader:
            manifest_fd, manifest_bytes, manifest_artifact = reader.open(
                TEST_EVIDENCE_MANIFEST, include_mode=True
            )
            signature_fd, _, signature_artifact = reader.open(
                TEST_EVIDENCE_SIGNATURE, include_mode=True
            )
            public_key_fd, _, public_key_artifact = reader.open(
                TEST_EVIDENCE_PUBLIC_KEY, include_mode=True
            )
            _, xml_bytes, xml_artifact = reader.open(
                CANONICAL_PYTEST_XML, include_mode=True
            )
            _, _, builder_artifact = reader.open(
                TEST_EVIDENCE_REPORT_BUILDER_PATH, include_mode=True
            )
            _, _, test_source_artifact = reader.open(
                TEST_EVIDENCE_TEST_SOURCE_PATH, include_mode=True
            )
            if (
                public_key_artifact["sha256"]
                != TEST_EVIDENCE_PUBLIC_KEY_SHA256
                or any(
                    record["mode"] != "0o444"
                    for record in (
                        manifest_artifact,
                        signature_artifact,
                        public_key_artifact,
                        xml_artifact,
                        builder_artifact,
                        test_source_artifact,
                    )
                )
            ):
                raise ReportBuildError(
                    "Supplemental test-evidence identity or mode drift"
                )
            if not authority.verify_ed25519_signature_fds(
                data_fd=manifest_fd,
                public_key_fd=public_key_fd,
                signature_fd=signature_fd,
            ):
                raise ReportBuildError(
                    "Supplemental test-evidence Ed25519 signature failed"
                )
            try:
                manifest = json.loads(manifest_bytes.decode("utf-8"))
                root = ET.fromstring(xml_bytes)
            except (
                UnicodeDecodeError,
                json.JSONDecodeError,
                ET.ParseError,
            ) as error:
                raise ReportBuildError(
                    "Signed supplemental test evidence is malformed"
                ) from error
            if not isinstance(manifest, dict):
                raise ReportBuildError(
                    "Signed supplemental test-evidence manifest is not an object"
                )
            verification = validate_pytest_root(root)
            expected_artifacts = {
                "canonical_pytest_xml": xml_artifact,
                "report_builder": builder_artifact,
                "report_builder_tests": test_source_artifact,
            }
            expected_signature_contract = {
                "algorithm": "Ed25519",
                "public_key": public_key_artifact,
                "signed_artifact": str(TEST_EVIDENCE_MANIFEST.resolve()),
                "signature": str(TEST_EVIDENCE_SIGNATURE.resolve()),
                "private_key_lifecycle": (
                    "encrypted_one_time_key_chmod_000_after_signing"
                ),
                "private_key_path": str(TEST_EVIDENCE_PRIVATE_KEY.resolve()),
            }
            if (
                manifest.get("schema_name")
                != "intersubmod.supplemental_report_test_evidence"
                or manifest.get("schema_version") != "1.0.0"
                or manifest.get("authority_id")
                != TEST_EVIDENCE_AUTHORITY_ID
                or manifest.get("task_type") != "B_comprehensive_validation"
                or manifest.get("scope")
                != "positional_singleton_supplemental_report_release"
                or manifest.get("artifacts") != expected_artifacts
                or manifest.get("canonical_test_summary") != verification
                or manifest.get("required_test_names")
                != sorted(REQUIRED_SUPPLEMENTAL_TESTS)
                or manifest.get("signature_contract")
                != expected_signature_contract
                or manifest.get("checks")
                != {
                    "canonical_pytest_pass": True,
                    "one_time_private_key_retired_after_signature": True,
                    "report_builder_identity_bound": True,
                    "report_builder_tests_identity_bound": True,
                }
                or manifest.get("pass") is not True
            ):
                raise ReportBuildError(
                    "Signed supplemental test-evidence contract drift"
                )
            reader.require_paths_still_bound()
    except (
        FINAL_RELEASE.SOURCE_AUTHORITY.SourceAuthorityError,
        OSError,
    ) as error:
        raise ReportBuildError(
            f"Unable to validate signed supplemental test evidence: {error}"
        ) from error
    return {
        **verification,
        "authority_id": TEST_EVIDENCE_AUTHORITY_ID,
        "signed_evidence_verified": True,
        "evidence_manifest": manifest_artifact,
        "evidence_signature": signature_artifact,
        "evidence_public_key": public_key_artifact,
    }


def verify_signature(
    *,
    receipt_path: Path,
    signature_path: Path,
    public_key_path: Path,
    expected_public_key_sha256: str,
    label: str,
) -> None:
    for path in (receipt_path, signature_path, public_key_path):
        if stat.S_IMODE(path.stat().st_mode) != 0o444:
            raise ReportBuildError(f"{label} protected mode drift: {path}")
    if sha256(public_key_path) != expected_public_key_sha256:
        raise ReportBuildError(f"{label} public-key SHA-256 drift")
    completed = subprocess.run(
        [
            str(OPENSSL),
            "pkeyutl",
            "-verify",
            "-rawin",
            "-pubin",
            "-inkey",
            str(public_key_path),
            "-sigfile",
            str(signature_path),
            "-in",
            str(receipt_path),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise ReportBuildError(
            f"{label} detached signature failed: {completed.stderr.strip()}"
        )


def fsync_file(path: Path) -> None:
    with path.open("rb") as handle:
        os.fsync(handle.fileno())


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def rename_no_replace(source: Path, target: Path) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        raise ReportBuildError("renameat2 is required for no-replace publication")
    renameat2.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    result = renameat2(
        -100,
        os.fsencode(source),
        -100,
        os.fsencode(target),
        1,
    )
    if result != 0:
        error_number = ctypes.get_errno()
        raise FileExistsError(
            error_number,
            f"No-replace publication failed: {os.strerror(error_number)}",
            str(target),
        )


def pct(numerator: int, denominator: int, digits: int = 2) -> str:
    if denominator <= 0:
        raise ReportBuildError(
            f"Percentage denominator must be positive: {denominator}"
        )
    return f"{100.0 * numerator / denominator:.{digits}f}%"


def ratio(numerator: int, denominator: int) -> float:
    if denominator <= 0:
        raise ReportBuildError(
            f"Ratio denominator must be positive: {denominator}"
        )
    return numerator / denominator


def wilson_interval(
    numerator: int, denominator: int, *, z: float = 1.959963984540054
) -> tuple[float, float]:
    proportion = numerator / denominator
    correction = z * z / denominator
    center = (proportion + correction / 2) / (1 + correction)
    half_width = (
        z
        * math.sqrt(
            proportion * (1 - proportion) / denominator
            + z * z / (4 * denominator * denominator)
        )
        / (1 + correction)
    )
    return center - half_width, center + half_width


def rate_value(summary: Mapping[str, Any], name: str) -> float:
    rate = summary["rates"][name]
    if not isinstance(rate, dict):
        raise ReportBuildError(f"Rate is not structured: {name}")
    return float(rate["value"])


def validate_summary_contract(summary: Mapping[str, Any]) -> None:
    counts = summary.get("counts")
    if not isinstance(counts, dict):
        raise ReportBuildError("Singleton summary counts are missing")
    required_count_fields = {
        "singleton_sites",
        "m1_evaluable",
        "m1_not_evaluable",
        "m1_flagged",
        "m2_pass",
        "m2_fail",
        "m2_not_evaluable",
        "m2_not_run",
    }
    if set(counts) != required_count_fields or any(
        not isinstance(counts[field], int) or counts[field] < 0
        for field in required_count_fields
    ):
        raise ReportBuildError("Singleton summary count schema drift")
    if (
        counts["m1_evaluable"] + counts["m1_not_evaluable"]
        != counts["singleton_sites"]
    ):
        raise ReportBuildError("M1 evaluability conservation failed")
    if counts["m1_flagged"] > counts["m1_evaluable"]:
        raise ReportBuildError("M1 flags exceed M1-evaluable sites")
    if (
        counts["m2_pass"]
        + counts["m2_fail"]
        + counts["m2_not_evaluable"]
        != counts["m1_flagged"]
    ):
        raise ReportBuildError("M2 within-M1 conservation failed")
    if (
        counts["m1_flagged"] + counts["m2_not_run"]
        != counts["singleton_sites"]
    ):
        raise ReportBuildError("M2 full-singleton conservation failed")

    expected_status = {
        "m1": {
            "FLAGGED": counts["m1_flagged"],
            "NOT_FLAGGED": counts["m2_not_run"],
        },
        "m2": {
            "PASS": counts["m2_pass"],
            "FAIL": counts["m2_fail"],
            "NOT_EVALUABLE": counts["m2_not_evaluable"],
            "NOT_RUN": counts["m2_not_run"],
        },
        "g1": {"NOT_RUN": counts["singleton_sites"]},
        "g2": {"NOT_RUN": counts["singleton_sites"]},
        "r1": {"NOT_RUN": counts["singleton_sites"]},
    }
    if summary.get("status_census") != expected_status:
        raise ReportBuildError("Singleton status census drift")

    distribution = summary.get("methyl_group_count_distribution")
    if (
        not isinstance(distribution, dict)
        or any(
            not str(group_count).isdigit()
            or not 2 <= int(group_count) <= 10
            or not isinstance(site_count, int)
            or site_count < 0
            for group_count, site_count in distribution.items()
        )
        or sum(distribution.values()) != counts["m1_flagged"]
    ):
        raise ReportBuildError("M1 methyl-group distribution drift")

    screen = summary.get("screen_metadata")
    if not isinstance(screen, dict):
        raise ReportBuildError("Screen metadata are missing")
    all_rows = screen.get("all_rows")
    if (
        screen.get("datasets") != list(EXPECTED_DATASETS)
        or screen.get("unique_keys") != all_rows
        or screen.get("duplicate_keys") != 0
        or not isinstance(all_rows, int)
        or all_rows <= 0
    ):
        raise ReportBuildError("Screen universe contract drift")
    branch_counts = screen.get("branch_counts")
    if (
        not isinstance(branch_counts, dict)
        or set(branch_counts)
        != {"max_snv_excluded", "positional_singleton", "retained"}
        or any(
            not isinstance(value, int) or value < 0
            for value in branch_counts.values()
        )
        or branch_counts.get("positional_singleton")
        != counts["singleton_sites"]
        or sum(branch_counts.values()) != all_rows
    ):
        raise ReportBuildError("Screen branch conservation failed")
    biological_ids = screen.get("dataset_biological_ids")
    if (
        not isinstance(biological_ids, dict)
        or set(biological_ids) != set(EXPECTED_DATASETS)
        or len(set(biological_ids.values())) != 6
        or biological_ids.get("HCC1395")
        != biological_ids.get("HCC1395_DORADO")
    ):
        raise ReportBuildError("Biological-sample identity contract drift")

    positional = summary.get("positional_recomputation")
    if (
        not isinstance(positional, dict)
        or positional.get("component_rows") != all_rows
        or positional.get("recomputed_singletons")
        != counts["singleton_sites"]
        or positional.get("component_identity_mismatch") != 0
        or positional.get("singleton_branch_mismatch") != 0
        or positional.get("singleton_positional_contract_failure") != 0
        or not isinstance(
            positional.get("minimum_finite_singleton_nearest_gap_bp"), int
        )
        or positional["minimum_finite_singleton_nearest_gap_bp"] <= 50_000
    ):
        raise ReportBuildError("Positional singleton recomputation drift")

    total_fields = {
        "sites": counts["singleton_sites"],
        "m1_evaluable": counts["m1_evaluable"],
        "m1_flagged": counts["m1_flagged"],
        "m2_pass": counts["m2_pass"],
        "m2_fail": counts["m2_fail"],
        "m2_not_evaluable": counts["m2_not_evaluable"],
        "m2_not_run": counts["m2_not_run"],
    }
    breakdowns = summary.get("breakdowns")
    if not isinstance(breakdowns, dict):
        raise ReportBuildError("Singleton breakdowns are missing")

    def validate_partition(
        partition: Any,
        expected_keys: set[str],
        *,
        label: str,
    ) -> None:
        if not isinstance(partition, dict) or set(partition) != expected_keys:
            raise ReportBuildError(f"{label} key drift")
        for key, values in partition.items():
            if (
                not isinstance(values, dict)
                or set(values) != set(total_fields)
                or any(
                    not isinstance(values[field], int) or values[field] < 0
                    for field in total_fields
                )
                or values["m1_evaluable"] > values["sites"]
                or values["m1_flagged"]
                != values["m2_pass"]
                + values["m2_fail"]
                + values["m2_not_evaluable"]
                or values["sites"]
                != values["m1_flagged"] + values["m2_not_run"]
            ):
                raise ReportBuildError(f"{label} row drift at {key}")
        for field, expected_total in total_fields.items():
            if sum(values[field] for values in partition.values()) != expected_total:
                raise ReportBuildError(
                    f"{label} conservation failed for {field}"
                )

    dataset_partition = breakdowns.get("dataset")
    truth_partition = breakdowns.get("truth")
    validate_partition(
        dataset_partition, set(EXPECTED_DATASETS), label="Dataset breakdown"
    )
    validate_partition(
        truth_partition, {"TP", "FP", "UNASSESSED"}, label="Truth breakdown"
    )
    dataset_truth = breakdowns.get("dataset_truth")
    expected_dataset_truth = {
        f"{dataset}|{truth}"
        for dataset in EXPECTED_DATASETS
        for truth in ("TP", "FP", "UNASSESSED")
    }
    if not isinstance(dataset_truth, dict) or set(dataset_truth) != expected_dataset_truth:
        raise ReportBuildError("Dataset-by-truth breakdown key drift")
    for key, values in dataset_truth.items():
        if (
            not isinstance(values, dict)
            or set(values) != set(total_fields)
            or any(
                not isinstance(values[field], int) or values[field] < 0
                for field in total_fields
            )
            or values["m1_evaluable"] > values["sites"]
            or values["m1_flagged"]
            != values["m2_pass"]
            + values["m2_fail"]
            + values["m2_not_evaluable"]
            or values["sites"]
            != values["m1_flagged"] + values["m2_not_run"]
        ):
            raise ReportBuildError(
                f"Dataset-by-truth breakdown row drift at {key}"
            )
    for dataset in EXPECTED_DATASETS:
        for field in total_fields:
            if (
                sum(
                    dataset_truth[f"{dataset}|{truth}"][field]
                    for truth in ("TP", "FP", "UNASSESSED")
                )
                != dataset_partition[dataset][field]
            ):
                raise ReportBuildError(
                    f"Dataset-by-truth to dataset drift: {dataset} {field}"
                )
    for truth in ("TP", "FP", "UNASSESSED"):
        for field in total_fields:
            if (
                sum(
                    dataset_truth[f"{dataset}|{truth}"][field]
                    for dataset in EXPECTED_DATASETS
                )
                != truth_partition[truth][field]
            ):
                raise ReportBuildError(
                    f"Dataset-by-truth to truth drift: {truth} {field}"
                )

    dataset_counts = screen.get("dataset_counts")
    if not isinstance(dataset_counts, dict) or set(dataset_counts) != set(
        EXPECTED_DATASETS
    ):
        raise ReportBuildError("All-site dataset census drift")
    required_dataset_count_fields = {
        "all_ssnv",
        "branch_max_snv_excluded",
        "branch_positional_singleton",
        "branch_retained",
        "truth_fp",
        "truth_tp",
        "truth_unassessed",
    }
    if any(
        not isinstance(value, dict)
        or not required_dataset_count_fields.issubset(value)
        or any(
            not isinstance(value[field], int) or value[field] < 0
            for field in required_dataset_count_fields
        )
        for value in dataset_counts.values()
    ):
        raise ReportBuildError("All-site dataset row schema drift")
    if (
        sum(value["all_ssnv"] for value in dataset_counts.values()) != all_rows
        or sum(
            value["branch_positional_singleton"]
            for value in dataset_counts.values()
        )
        != counts["singleton_sites"]
        or any(
            value["all_ssnv"]
            != value["branch_max_snv_excluded"]
            + value["branch_positional_singleton"]
            + value["branch_retained"]
            for value in dataset_counts.values()
        )
    ):
        raise ReportBuildError("All-site dataset conservation failed")
    truth_counts = screen.get("truth_counts")
    if (
        not isinstance(truth_counts, dict)
        or set(truth_counts) != {"TP", "FP", "UNASSESSED"}
        or any(
            not isinstance(value, int) or value < 0
            for value in truth_counts.values()
        )
        or sum(truth_counts.values()) != all_rows
    ):
        raise ReportBuildError("All-site truth census drift")


def resolve_source_authority_input_paths(
    summary: Mapping[str, Any],
) -> dict[str, Path]:
    """Resolve every transitive artifact in the signed v4 audit chain."""

    source_chain = summary.get("source_chain")
    if not isinstance(source_chain, dict):
        raise ReportBuildError("Source chain is missing")
    authority_chain = source_chain.get("v4_source_authority")
    if (
        not isinstance(authority_chain, dict)
        or authority_chain.get("pass") is not True
        or authority_chain.get("authority_id") != SOURCE_AUTHORITY_ID
    ):
        raise ReportBuildError("Source-authority chain is incomplete")
    records = {
        "v4_authority_manifest": authority_chain.get("authority_manifest"),
        "v4_authority_approval": authority_chain.get("detached_approval"),
        "v4_authority_signature": authority_chain.get(
            "detached_approval_signature"
        ),
        "v4_authority_public_key": authority_chain.get("public_key"),
    }
    paths = {}
    for name, record in records.items():
        if not isinstance(record, dict):
            raise ReportBuildError(f"{name} identity is missing")
        paths[name] = Path(str(record.get("path", ""))).resolve(strict=True)
    code = summary.get("code")
    if not isinstance(code, dict):
        raise ReportBuildError("Singleton audit code identities are missing")
    supplemental_record = code.get("supplemental_auditor")
    if not isinstance(supplemental_record, dict):
        raise ReportBuildError("Supplemental auditor identity is missing")
    paths["supplemental_auditor"] = Path(
        str(supplemental_record.get("path", ""))
    ).resolve(strict=True)
    return paths


def validate_source_authority_chain(summary: Mapping[str, Any]) -> None:
    source_chain = summary.get("source_chain")
    if not isinstance(source_chain, dict):
        raise ReportBuildError("Source chain is missing")
    authority_chain = source_chain.get("v4_source_authority")
    if not isinstance(authority_chain, dict):
        raise ReportBuildError("Source-authority chain is incomplete")
    paths = resolve_source_authority_input_paths(summary)
    records = {
        "v4_authority_manifest": authority_chain.get("authority_manifest"),
        "v4_authority_approval": authority_chain.get("detached_approval"),
        "v4_authority_signature": authority_chain.get(
            "detached_approval_signature"
        ),
        "v4_authority_public_key": authority_chain.get("public_key"),
        "supplemental_auditor": summary.get("code", {}).get(
            "supplemental_auditor"
        ),
    }
    for name, path in paths.items():
        require_identity(path, records[name], label=name)
    verify_signature(
        receipt_path=paths["v4_authority_approval"],
        signature_path=paths["v4_authority_signature"],
        public_key_path=paths["v4_authority_public_key"],
        expected_public_key_sha256=SOURCE_AUTHORITY_PUBLIC_KEY_SHA256,
        label="source authority",
    )
    authority = load_json(
        paths["v4_authority_manifest"], label="source authority manifest"
    )
    approval = load_json(
        paths["v4_authority_approval"], label="source authority approval"
    )
    expected_contract = (
        authority.get("schema_name")
        == "intersubmod.release_source_authority"
        and authority.get("schema_version") == "1.0.0"
        and authority.get("authority_id") == SOURCE_AUTHORITY_ID
        and authority.get("approval_status")
        == "APPROVED_FOR_FULL_TASK_B_RUN"
        and authority.get("source_set_sha256")
        == authority_chain.get("source_set_sha256")
        and authority.get("repository", {}).get("git_head_at_authorization")
        == authority_chain.get("authorized_git_head")
        and approval.get("schema_name")
        == "intersubmod.release_source_authority.approval"
        and approval.get("schema_version") == "1.0.0"
        and approval.get("authority_id") == SOURCE_AUTHORITY_ID
        and approval.get("approval_status")
        == "APPROVED_FOR_FULL_TASK_B_RUN"
        and all(
            approval.get("authority_manifest", {}).get(field)
            == authority_chain.get("authority_manifest", {}).get(field)
            for field in ("path", "size_bytes", "sha256", "mode")
        )
        and all(
            approval.get("public_key", {}).get(field)
            == authority_chain.get("public_key", {}).get(field)
            for field in ("path", "size_bytes", "sha256", "mode")
        )
        and len(approval.get("review_approvals", [])) >= 2
        and all(
            review.get("verdict") == "APPROVE"
            for review in approval.get("review_approvals", [])
        )
    )
    if not expected_contract:
        raise ReportBuildError("Source-authority signed contract drift")


def validate_formal_release_chain(
    *,
    final_dataset_receipt_path: Path,
    final_dataset_signature_path: Path,
    final_dataset_public_key_path: Path,
    final_report_receipt_path: Path,
    final_report_signature_path: Path,
    final_report_public_key_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    expected_paths = {
        "final dataset receipt": FINAL_RELEASE.RELEASE_RECEIPT,
        "final dataset signature": FINAL_RELEASE.RELEASE_SIGNATURE,
        "final dataset public key": FINAL_RELEASE.RESULT_PUBLIC_KEY,
        "final report receipt": FINAL_RELEASE.REPORT_RELEASE_RECEIPT,
        "final report signature": FINAL_RELEASE.REPORT_RELEASE_SIGNATURE,
        "final report public key": FINAL_RELEASE.REPORT_PUBLIC_KEY,
    }
    observed_paths = {
        "final dataset receipt": final_dataset_receipt_path,
        "final dataset signature": final_dataset_signature_path,
        "final dataset public key": final_dataset_public_key_path,
        "final report receipt": final_report_receipt_path,
        "final report signature": final_report_signature_path,
        "final report public key": final_report_public_key_path,
    }
    mismatches = {
        label: {
            "observed": str(observed.resolve()),
            "expected": str(expected.resolve()),
        }
        for label, observed in observed_paths.items()
        for expected in (expected_paths[label],)
        if observed.resolve() != expected.resolve()
    }
    if mismatches:
        raise ReportBuildError(
            f"Formal release paths are not canonical: {mismatches}"
        )
    try:
        dataset_verification = (
            FINAL_RELEASE.validate_release_signature_artifacts()
        )
        report_verification = (
            FINAL_RELEASE.validate_report_release_signature_artifacts()
        )
    except (
        FINAL_RELEASE.FinalReleaseError,
        FINAL_RELEASE.SOURCE_AUTHORITY.SourceAuthorityError,
        OSError,
        ValueError,
    ) as error:
        raise ReportBuildError(
            f"Formal Task-B signed release chain failed: {error}"
        ) from error
    if (
        dataset_verification.get("schema_name")
        != "intersubmod.task_b_final_dataset_release_receipt.verification"
        or dataset_verification.get("signature_verified") is not True
        or dataset_verification.get("all_final_outputs_reverified") is not True
        or dataset_verification.get("pass") is not True
        or report_verification.get("schema_name")
        != "intersubmod.task_b_final_report_release_receipt.verification"
        or report_verification.get("signed_dataset_release_reverified")
        is not True
        or report_verification.get("signature_verified") is not True
        or report_verification.get("all_report_outputs_reverified") is not True
        or report_verification.get("pass") is not True
    ):
        raise ReportBuildError("Formal Task-B release verification contract drift")
    if (
        dataset_verification.get("source_authority")
        != report_verification.get("source_authority")
    ):
        raise ReportBuildError(
            "Dataset/report source-authority verification differs"
        )
    return dataset_verification, report_verification


def validate_inputs(
    *,
    summary_path: Path,
    site_audit_path: Path,
    candidate_path: Path,
    assignments_path: Path,
    final_dataset_receipt_path: Path,
    final_dataset_signature_path: Path,
    final_dataset_public_key_path: Path,
    final_report_receipt_path: Path,
    final_report_signature_path: Path,
    final_report_public_key_path: Path,
    tumor_ref_controls_path: Path,
) -> tuple[
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
]:
    summary = load_json(summary_path, label="singleton audit summary")
    if (
        summary.get("schema_name") != AUDIT_SCHEMA
        or summary.get("schema_version") != "1.0.0"
        or summary.get("pass") is not True
        or not isinstance(summary.get("checks"), dict)
        or not summary["checks"]
        or not all(value is True for value in summary["checks"].values())
    ):
        raise ReportBuildError("Singleton audit is not a complete PASS")
    observed = summary.get("counts", {})
    expected = {
        "singleton_sites": EXPECTED_COUNTS["singleton_sites"],
        "m1_evaluable": EXPECTED_COUNTS["m1_evaluable"],
        "m1_flagged": EXPECTED_COUNTS["m1_flagged"],
        "m2_pass": EXPECTED_COUNTS["m2_pass"],
        "m2_fail": EXPECTED_COUNTS["m2_fail"],
        "m2_not_evaluable": EXPECTED_COUNTS["m2_not_evaluable"],
        "m2_not_run": EXPECTED_COUNTS["m2_not_run"],
    }
    mismatches = {
        key: {"observed": observed.get(key), "expected": value}
        for key, value in expected.items()
        if observed.get(key) != value
    }
    if mismatches:
        raise ReportBuildError(f"Singleton count drift: {mismatches}")
    if summary.get("screen_metadata", {}).get("all_rows") != EXPECTED_COUNTS[
        "all_dataset_sites"
    ]:
        raise ReportBuildError("All-site denominator drift")
    validate_summary_contract(summary)
    validate_source_authority_chain(summary)
    require_identity(
        site_audit_path,
        summary.get("outputs", {}).get("site_audit"),
        label="singleton site audit",
    )
    require_identity(
        candidate_path,
        summary.get("outputs", {}).get("m2_pass_cases"),
        label="singleton M2 cases",
    )
    require_identity(
        assignments_path,
        summary.get("inputs", {}).get("stable_assignments"),
        label="stable assignments",
    )
    dataset_verification, report_verification = validate_formal_release_chain(
        final_dataset_receipt_path=final_dataset_receipt_path,
        final_dataset_signature_path=final_dataset_signature_path,
        final_dataset_public_key_path=final_dataset_public_key_path,
        final_report_receipt_path=final_report_receipt_path,
        final_report_signature_path=final_report_signature_path,
        final_report_public_key_path=final_report_public_key_path,
    )
    dataset_receipt = load_json(
        final_dataset_receipt_path, label="final dataset release receipt"
    )
    report_receipt = load_json(
        final_report_receipt_path, label="final report release receipt"
    )
    if (
        dataset_receipt.get("schema_name")
        != "intersubmod.task_b_final_dataset_release_receipt"
        or dataset_receipt.get("schema_version") != "1.0.0"
        or dataset_receipt.get("pass") is not True
    ):
        raise ReportBuildError("Final dataset release receipt contract drift")
    if (
        report_receipt.get("schema_name")
        != "intersubmod.task_b_final_report_release_receipt"
        or report_receipt.get("schema_version") != "1.0.0"
        or report_receipt.get("pass") is not True
    ):
        raise ReportBuildError("Final report release receipt contract drift")
    builder_record = dataset_receipt.get("inputs", {}).get("builder_receipt")
    if not isinstance(builder_record, dict):
        raise ReportBuildError("Signed dataset receipt lacks builder receipt identity")
    builder_receipt_path = Path(str(builder_record.get("path", ""))).resolve(
        strict=True
    )
    require_identity(
        builder_receipt_path,
        builder_record,
        label="signed final dataset builder receipt",
    )
    builder_receipt = load_json(
        builder_receipt_path, label="signed final dataset builder receipt"
    )
    if (
        builder_receipt.get("schema_name")
        != "intersubmod.all_ssnv_final_report_dataset_run_receipt"
        or builder_receipt.get("task_type") != "B_comprehensive_validation"
        or builder_receipt.get("formal_task_b_release_eligible") is not True
        or builder_receipt.get("pass") is not True
    ):
        raise ReportBuildError("Signed final dataset builder receipt contract drift")
    require_identity(
        tumor_ref_controls_path,
        builder_receipt.get("inputs", {}).get("tumor_ref_sites"),
        label="signed tumor-REF controls",
    )
    return (
        summary,
        dataset_receipt,
        report_receipt,
        builder_receipt,
        dataset_verification,
        report_verification,
    )


def load_candidates(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != EXPECTED_COUNTS["m2_pass"]:
        raise ReportBuildError(
            f"Expected {EXPECTED_COUNTS['m2_pass']} M2 PASS cases, "
            f"observed {len(rows)}"
        )
    keys: set[SiteKey] = set()
    for row in rows:
        key = (
            row["dataset"],
            row["chrom"],
            int(row["pos"]),
            row["ref"],
            row["alt"],
        )
        if key in keys:
            raise ReportBuildError(f"Duplicate M2 PASS case: {key}")
        keys.add(key)
        row["pos"] = int(row["pos"])
        row["methyl_group_count"] = int(row["methyl_group_count"])
        row["core_read_n"] = int(row["core_read_n"])
        row["min_group_n"] = int(row["min_group_n"])
    return rows


def repository_metadata(summary: Mapping[str, Any]) -> dict[str, str]:
    repository = Path(__file__).resolve().parents[3]

    def git_value(*arguments: str) -> str:
        completed = subprocess.run(
            ["/usr/bin/git", "-C", str(repository), *arguments],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            raise ReportBuildError(
                f"Unable to read repository metadata: {completed.stderr.strip()}"
            )
        return completed.stdout.strip()

    head = git_value("rev-parse", "HEAD")
    branch = git_value("rev-parse", "--abbrev-ref", "HEAD")
    authorized_head = str(
        summary["source_chain"]["v4_source_authority"][
            "authorized_git_head"
        ]
    )
    if head != authorized_head:
        raise ReportBuildError(
            f"Repository HEAD drift: {head} != {authorized_head}"
        )
    return {
        "worktree": str(repository),
        "branch": branch,
        "head": head,
        "authorized_head": authorized_head,
    }


def load_m1_singleton_keys(path: Path) -> set[SiteKey]:
    keys: set[SiteKey] = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "dataset",
            "chrom",
            "pos",
            "ref",
            "alt",
            "m1_status",
        }
        if reader.fieldnames is None or not required.issubset(
            reader.fieldnames
        ):
            raise ReportBuildError("Singleton site audit lacks M1 key fields")
        for row in reader:
            if row["m1_status"] != "FLAGGED":
                continue
            key = (
                row["dataset"],
                row["chrom"],
                int(row["pos"]),
                row["ref"],
                row["alt"],
            )
            if key in keys:
                raise ReportBuildError(f"Duplicate singleton M1 key: {key}")
            keys.add(key)
    if len(keys) != EXPECTED_COUNTS["m1_flagged"]:
        raise ReportBuildError(
            f"Singleton M1 key count drift: {len(keys)}"
        )
    return keys


def parse_lower_bool(value: str, *, field: str, key: SiteKey) -> bool:
    if value not in {"true", "false"}:
        raise ReportBuildError(
            f"Invalid tumor-REF boolean {field}={value!r} at {key}"
        )
    return value == "true"


def load_tumor_ref_controls(
    path: Path,
    *,
    m1_keys: set[SiteKey],
    m2_keys: set[SiteKey],
) -> tuple[dict[SiteKey, dict[str, Any]], dict[str, dict[str, Any]]]:
    if not m2_keys.issubset(m1_keys):
        raise ReportBuildError("M2 PASS keys are not a subset of singleton M1 keys")
    controls: dict[SiteKey, dict[str, Any]] = {}
    required = {
        "dataset",
        "chrom",
        "pos",
        "ref",
        "alt",
        "n_tumor_ref",
        "ref_status",
        "ref_evaluable",
        "ref_stable_null_multigroup",
        "joint_stable_null_multigroup",
        "joint_allele_axis_aligned",
        "background_control_interpretation",
    }
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(
            reader.fieldnames
        ):
            raise ReportBuildError("Tumor-REF controls lack required fields")
        for row in reader:
            key = (
                row["dataset"],
                row["chrom"],
                int(row["pos"]),
                row["ref"],
                row["alt"],
            )
            if key not in m1_keys:
                continue
            if key in controls:
                raise ReportBuildError(f"Duplicate tumor-REF control: {key}")
            evaluable = parse_lower_bool(
                row["ref_evaluable"], field="ref_evaluable", key=key
            )
            stable = parse_lower_bool(
                row["ref_stable_null_multigroup"],
                field="ref_stable_null_multigroup",
                key=key,
            )
            joint_stable = parse_lower_bool(
                row["joint_stable_null_multigroup"],
                field="joint_stable_null_multigroup",
                key=key,
            )
            allele_aligned = parse_lower_bool(
                row["joint_allele_axis_aligned"],
                field="joint_allele_axis_aligned",
                key=key,
            )
            if stable and not evaluable:
                raise ReportBuildError(
                    f"Non-evaluable tumor-REF control is stable at {key}"
                )
            if evaluable:
                expected_interpretation = (
                    "REF_REPLICATION_WEAKENS_ALT_SPECIFICITY"
                    if stable
                    else "REF_NONREPLICATION_DOES_NOT_PROVE_SUBCLONE"
                )
            else:
                expected_interpretation = "REF_CONTROL_NOT_TESTABLE"
            if (
                row["background_control_interpretation"]
                != expected_interpretation
            ):
                raise ReportBuildError(
                    f"Tumor-REF interpretation drift at {key}"
                )
            controls[key] = {
                "n_tumor_ref": int(row["n_tumor_ref"]),
                "ref_status": row["ref_status"],
                "ref_evaluable": evaluable,
                "ref_stable_multigroup": stable,
                "joint_stable_multigroup": joint_stable,
                "joint_allele_axis_aligned": allele_aligned,
                "interpretation": expected_interpretation,
            }
    missing = sorted(m1_keys.difference(controls))
    if missing:
        raise ReportBuildError(
            f"Singleton M1 sites lack tumor-REF controls: {missing[:3]}"
        )

    def summarize(keys: set[SiteKey]) -> dict[str, Any]:
        rows = [controls[key] for key in keys]
        evaluable = sum(row["ref_evaluable"] for row in rows)
        stable = sum(row["ref_stable_multigroup"] for row in rows)
        nonreplication = evaluable - stable
        not_evaluable = len(rows) - evaluable
        return {
            "sites": len(rows),
            "ref_evaluable": evaluable,
            "ref_stable_multigroup": stable,
            "ref_nonreplication": nonreplication,
            "ref_not_evaluable": not_evaluable,
            "joint_stable_multigroup": sum(
                row["joint_stable_multigroup"] for row in rows
            ),
            "joint_allele_axis_aligned": sum(
                row["joint_allele_axis_aligned"] for row in rows
            ),
        }

    return controls, {
        "m1": summarize(m1_keys),
        "m2_pass": summarize(m2_keys),
    }


def technical_pair_overlap(path: Path) -> dict[str, dict[str, Any]]:
    datasets = ("HCC1395", "HCC1395_DORADO")
    sets: dict[str, dict[str, set[tuple[str, int, str, str]]]] = {
        dataset: {"all": set(), "m1": set(), "m2": set()}
        for dataset in datasets
    }
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            dataset = row["dataset"]
            if dataset not in sets:
                continue
            key = (
                row["chrom"],
                int(row["pos"]),
                row["ref"],
                row["alt"],
            )
            sets[dataset]["all"].add(key)
            if row["m1_status"] == "FLAGGED":
                sets[dataset]["m1"].add(key)
            if row["m2_status"] == "PASS":
                sets[dataset]["m2"].add(key)
    result: dict[str, dict[str, Any]] = {}
    for level in ("all", "m1", "m2"):
        left = sets[datasets[0]][level]
        right = sets[datasets[1]][level]
        union = left.union(right)
        intersection = left.intersection(right)
        result[level] = {
            "hcc1395": len(left),
            "hcc1395_dorado": len(right),
            "intersection": len(intersection),
            "union": len(union),
            "jaccard": len(intersection) / len(union) if union else None,
        }
    expected = {
        "all": (7_484, 9_116),
        "m1": (407, 1_289),
        "m2": (0, 4),
    }
    for level, (intersection, union) in expected.items():
        if (
            result[level]["intersection"] != intersection
            or result[level]["union"] != union
        ):
            raise ReportBuildError(f"Technical-pair overlap drift at {level}")
    return result


def assignment_key(record: Mapping[str, Any]) -> SiteKey:
    posthoc = record.get("posthoc")
    if not isinstance(posthoc, dict):
        raise ReportBuildError("Stable assignment lacks posthoc REF/ALT")
    return (
        str(record.get("dataset") or record.get("sample")),
        str(record["chrom"]),
        int(record["pos"]),
        str(posthoc["ref"]),
        str(posthoc["alt"]),
    )


def load_candidate_assignments(
    path: Path, candidate_keys: set[SiteKey]
) -> dict[SiteKey, dict[str, Any]]:
    result: dict[SiteKey, dict[str, Any]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            key = assignment_key(record)
            if key not in candidate_keys:
                continue
            if key in result:
                raise ReportBuildError(f"Duplicate candidate assignment: {key}")
            result[key] = record
    missing = sorted(candidate_keys.difference(result))
    if missing:
        raise ReportBuildError(f"M2 candidates lack assignments: {missing[:3]}")
    return result


def summarize_phase_sets(
    candidates: Iterable[Mapping[str, Any]],
    assignments: Mapping[SiteKey, Mapping[str, Any]],
) -> dict[str, Any]:
    sites: list[dict[str, Any]] = []
    for candidate in candidates:
        key = (
            str(candidate["dataset"]),
            str(candidate["chrom"]),
            int(candidate["pos"]),
            str(candidate["ref"]),
            str(candidate["alt"]),
        )
        assignment = assignments.get(key)
        if not isinstance(assignment, Mapping):
            raise ReportBuildError(f"Candidate phase-set assignment is missing: {key}")
        core_reads = assignment.get("core_reads")
        if not isinstance(core_reads, list) or not core_reads:
            raise ReportBuildError(f"Candidate phase-set core reads are missing: {key}")
        if len(core_reads) != int(candidate["core_read_n"]):
            raise ReportBuildError(
                f"Candidate phase-set core-read count drift: {key}"
            )
        ps_values: list[int | None] = []
        labels: list[str] = []
        for read in core_reads:
            if not isinstance(read, Mapping) or "latest_ps" not in read:
                raise ReportBuildError(
                    f"Candidate core read lacks latest_ps: {key}"
                )
            value = read["latest_ps"]
            if value is not None and (
                isinstance(value, bool) or not isinstance(value, int) or value < 0
            ):
                raise ReportBuildError(
                    f"Candidate core read has invalid latest_ps: {key}"
                )
            label = str(read.get("label", ""))
            if not label:
                raise ReportBuildError(
                    f"Candidate core read lacks methyl-group label: {key}"
                )
            ps_values.append(value)
            labels.append(label)
        nonmissing = sorted({value for value in ps_values if value is not None})
        by_group = {}
        for label in sorted(set(labels)):
            indices = [
                index for index, observed in enumerate(labels) if observed == label
            ]
            group_values = [ps_values[index] for index in indices]
            by_group[label] = {
                "core_read_n": len(indices),
                "ps_missing_n": sum(value is None for value in group_values),
                "ps_values": sorted(
                    {value for value in group_values if value is not None}
                ),
            }
        sites.append(
            {
                "dataset": key[0],
                "chrom": key[1],
                "pos": key[2],
                "ref": key[3],
                "alt": key[4],
                "core_read_n": len(core_reads),
                "ps_missing_n": sum(value is None for value in ps_values),
                "ps_nonmissing_n": sum(value is not None for value in ps_values),
                "ps_distinct_nonmissing_n": len(nonmissing),
                "ps_values": nonmissing,
                "by_methyl_group": by_group,
            }
        )
    sites.sort(
        key=lambda row: (
            EXPECTED_DATASETS.index(str(row["dataset"])),
            str(row["chrom"]),
            int(row["pos"]),
        )
    )
    return {
        "contract": (
            "latest_LongPhase_S_PS_descriptive_census_not_an_M2_axis_v1"
        ),
        "m2_phase_set_axis_evaluated": False,
        "site_count": len(sites),
        "total_core_reads": sum(row["core_read_n"] for row in sites),
        "total_missing_ps": sum(row["ps_missing_n"] for row in sites),
        "sites_with_missing_ps": sum(
            row["ps_missing_n"] > 0 for row in sites
        ),
        "sites_with_multiple_nonmissing_ps": sum(
            row["ps_distinct_nonmissing_n"] > 1 for row in sites
        ),
        "sites_with_missing_and_multiple_nonmissing_ps": sum(
            row["ps_missing_n"] > 0
            and row["ps_distinct_nonmissing_n"] > 1
            for row in sites
        ),
        "sites": sites,
    }


def validate_phase_set_audit(audit: Mapping[str, Any]) -> None:
    for field, expected in EXPECTED_PHASE_SET_AUDIT.items():
        if audit.get(field) != expected:
            raise ReportBuildError(
                f"Candidate phase-set audit drift at {field}: "
                f"{audit.get(field)} != {expected}"
            )
    target_key = EXPECTED_PHASE_SET_CASE["key"]
    targets = [
        row
        for row in audit.get("sites", [])
        if (
            row.get("dataset"),
            row.get("chrom"),
            row.get("pos"),
            row.get("ref"),
            row.get("alt"),
        )
        == target_key
    ]
    if len(targets) != 1:
        raise ReportBuildError("Required HCC1937 phase-set case is missing")
    target = targets[0]
    for field in (
        "core_read_n",
        "ps_missing_n",
        "ps_nonmissing_n",
        "ps_distinct_nonmissing_n",
        "ps_values",
    ):
        expected = EXPECTED_PHASE_SET_CASE[field]
        if target.get(field) != expected:
            raise ReportBuildError(
                f"Required HCC1937 phase-set case drift at {field}"
            )


def verify_matrix_identity(assignment: Mapping[str, Any]) -> Path:
    record = assignment.get("primary_artifacts", {}).get("methylation_matrix")
    if not isinstance(record, dict):
        raise ReportBuildError("Assignment lacks methylation matrix identity")
    raw_path = record.get("path")
    if not isinstance(raw_path, str) or not raw_path:
        raise ReportBuildError("Assignment methylation matrix path is missing")
    path = Path(raw_path)
    require_identity(path, record, label="candidate methylation matrix")
    return path


def select_cases(candidates: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    by_dataset: dict[str, list[Mapping[str, Any]]] = {}
    by_key: dict[SiteKey, Mapping[str, Any]] = {}
    for row in candidates:
        by_dataset.setdefault(str(row["dataset"]), []).append(row)
        key = (
            str(row["dataset"]),
            str(row["chrom"]),
            int(row["pos"]),
            str(row["ref"]),
            str(row["alt"]),
        )
        if key in by_key:
            raise ReportBuildError(f"Duplicate illustrated candidate: {key}")
        by_key[key] = row
    selected: list[dict[str, Any]] = []
    for dataset in EXPECTED_DATASETS:
        rows = by_dataset.get(dataset, [])
        if not rows:
            continue
        chosen = max(
            rows,
            key=lambda row: (
                int(row["core_read_n"]),
                int(row["methyl_group_count"]),
                -int(row["pos"]),
            ),
        )
        selected.append(dict(chosen))
    observed_datasets = {str(row["dataset"]) for row in selected}
    if observed_datasets != set(EXPECTED_POSITIVE_DATASETS):
        raise ReportBuildError(
            "Expected illustrated datasets "
            f"{EXPECTED_POSITIVE_DATASETS}, observed {sorted(observed_datasets)}"
        )
    selected_keys = {
        (
            str(row["dataset"]),
            str(row["chrom"]),
            int(row["pos"]),
            str(row["ref"]),
            str(row["alt"]),
        )
        for row in selected
    }
    for key in REQUIRED_VISUAL_CASES:
        required = by_key.get(key)
        if required is None:
            raise ReportBuildError(f"Required illustrated case is missing: {key}")
        if key not in selected_keys:
            selected.append(dict(required))
            selected_keys.add(key)
    selected.sort(
        key=lambda row: (
            EXPECTED_DATASETS.index(str(row["dataset"])),
            str(row["chrom"]),
            int(row["pos"]),
        )
    )
    return selected


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 9.5,
            "axes.edgecolor": COLORS["grid"],
            "axes.labelcolor": COLORS["ink"],
            "text.color": COLORS["ink"],
            "xtick.color": COLORS["muted"],
            "ytick.color": COLORS["muted"],
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def save_figure(figure: plt.Figure, path: Path) -> None:
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)
    path.chmod(0o444)
    fsync_file(path)
    if path.stat().st_size < 10_000:
        raise ReportBuildError(f"Figure is unexpectedly small: {path}")


def figure_funnel(summary: Mapping[str, Any], output: Path) -> None:
    counts = summary["counts"]
    labels = (
        "All latest LPS PASS sSNV",
        "Positional singleton",
        "M1 evaluable",
        "M1 stable multigroup",
        "M2 evaluable",
        "M2 PASS",
    )
    values = (
        int(summary["screen_metadata"]["all_rows"]),
        int(counts["singleton_sites"]),
        int(counts["m1_evaluable"]),
        int(counts["m1_flagged"]),
        int(counts["m2_pass"]) + int(counts["m2_fail"]),
        int(counts["m2_pass"]),
    )
    colors = (
        COLORS["ink"],
        COLORS["blue"],
        COLORS["teal"],
        COLORS["amber"],
        COLORS["red"],
        COLORS["teal"],
    )
    figure, axis = plt.subplots(figsize=(10.2, 5.1))
    y = np.arange(len(labels))
    axis.barh(y, values, color=colors, height=0.62)
    axis.set_xscale("log")
    axis.set_yticks(y, labels)
    axis.invert_yaxis()
    axis.set_xlabel("Dataset-sites (log scale)")
    axis.set_title(
        "The broad methylation screen narrows sharply at measured-axis M2",
        loc="left",
        fontweight="bold",
    )
    axis.grid(axis="x", color=COLORS["grid"], linewidth=0.7)
    axis.set_axisbelow(True)
    axis.spines[["top", "right", "left"]].set_visible(False)
    denominators = (
        values[0],
        values[0],
        values[1],
        values[1],
        values[3],
        values[1],
    )
    denominator_labels = (
        "all LPS sites",
        "all LPS sites",
        "singletons",
        "singletons",
        "M1 flags",
        "singletons",
    )
    for index, value in enumerate(values):
        denominator = denominators[index]
        label = f"{value:,}"
        if index:
            label += (
                f"  ({pct(value, denominator, 3)} of "
                f"{denominator_labels[index]})"
            )
        axis.text(value * 1.08, index, label, va="center", fontsize=9)
    figure.text(
        0.01,
        0.01,
        "M1 is an operational screen; M2 is a selected residual-partition gate. "
        "Neither is cellular-clone confirmation.",
        color=COLORS["muted"],
        fontsize=8.5,
    )
    figure.subplots_adjust(left=0.27, right=0.90, bottom=0.14, top=0.88)
    save_figure(figure, output)


def figure_dataset(summary: Mapping[str, Any], output: Path) -> None:
    breakdown = summary["breakdowns"]["dataset"]
    m1_rates = [
        ratio(
            breakdown[dataset]["m1_flagged"],
            breakdown[dataset]["sites"],
        )
        for dataset in EXPECTED_DATASETS
    ]
    m2_rates = [
        ratio(
            breakdown[dataset]["m2_pass"],
            breakdown[dataset]["sites"],
        )
        for dataset in EXPECTED_DATASETS
    ]
    figure, axes = plt.subplots(1, 2, figsize=(12.4, 5.2))
    y = np.arange(len(EXPECTED_DATASETS))
    axes[0].barh(y, m1_rates, color=COLORS["amber"], height=0.62)
    axes[0].set_yticks(y, EXPECTED_DATASETS)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("M1 flagged / singleton sites")
    axes[0].set_title("M1 screen yield varies by dataset", loc="left", fontweight="bold")
    axes[1].barh(y, m2_rates, color=COLORS["teal"], height=0.62)
    axes[1].set_yticks(y, EXPECTED_DATASETS)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("M2 PASS / singleton sites")
    axes[1].set_title(
        "M2 residual evidence remains rare", loc="left", fontweight="bold"
    )
    for axis, values, digits in ((axes[0], m1_rates, 2), (axes[1], m2_rates, 3)):
        axis.grid(axis="x", color=COLORS["grid"], linewidth=0.7)
        axis.set_axisbelow(True)
        axis.spines[["top", "right", "left"]].set_visible(False)
        maximum = max(values) if max(values) > 0 else 1
        axis.set_xlim(0, maximum * 1.32)
        for index, value in enumerate(values):
            axis.text(
                value + maximum * 0.025,
                index,
                f"{100 * value:.{digits}f}%",
                va="center",
                fontsize=8.5,
            )
    figure.text(
        0.01,
        0.01,
        "HCC1395 and HCC1395_DORADO are one biological sample; they are shown "
        "as a technical comparison, not independent replication.",
        color=COLORS["muted"],
        fontsize=8.5,
    )
    figure.subplots_adjust(left=0.12, right=0.98, bottom=0.15, top=0.87, wspace=0.40)
    save_figure(figure, output)


def figure_truth(summary: Mapping[str, Any], output: Path) -> None:
    breakdown = summary["breakdowns"]["truth"]
    labels = ("TP", "FP", "UNASSESSED")
    m1_all = [
        ratio(
            breakdown[label]["m1_flagged"],
            breakdown[label]["sites"],
        )
        for label in labels
    ]
    m1_evaluable = [
        ratio(
            breakdown[label]["m1_flagged"],
            breakdown[label]["m1_evaluable"],
        )
        for label in labels
    ]
    x = np.arange(len(labels))
    width = 0.34
    figure, axis = plt.subplots(figsize=(8.8, 5.0))
    first = axis.bar(
        x - width / 2,
        m1_all,
        width,
        color=COLORS["blue"],
        label="All singleton sites",
    )
    second = axis.bar(
        x + width / 2,
        m1_evaluable,
        width,
        color=COLORS["amber"],
        label="M1-evaluable singleton sites",
    )
    axis.set_xticks(x, labels)
    axis.set_ylabel("M1 stable multigroup rate")
    axis.set_title(
        "Truth strata are descriptive, not a clone-specificity benchmark",
        loc="left",
        fontweight="bold",
    )
    axis.grid(axis="y", color=COLORS["grid"], linewidth=0.7)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False)
    for bars in (first, second):
        for bar in bars:
            axis.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.003,
                f"{100 * bar.get_height():.2f}%",
                ha="center",
                fontsize=8.5,
            )
    axis.set_ylim(0, max(m1_all + m1_evaluable) * 1.28)
    figure.text(
        0.01,
        0.01,
        "Only one FP singleton is M2-evaluable in the selected subset; these "
        "data cannot estimate M2 specificity.",
        color=COLORS["red"],
        fontsize=8.5,
    )
    figure.subplots_adjust(bottom=0.16, top=0.87)
    save_figure(figure, output)


def figure_m2_status(summary: Mapping[str, Any], output: Path) -> None:
    status = summary["status_census"]["m2"]
    labels = ("NOT_RUN", "NOT_EVALUABLE", "FAIL", "PASS")
    values = [int(status[label]) for label in labels]
    colors = (
        COLORS["gray"],
        COLORS["amber"],
        COLORS["red"],
        COLORS["teal"],
    )
    figure, axis = plt.subplots(figsize=(9.5, 4.8))
    bars = axis.bar(labels, values, color=colors, width=0.62)
    axis.set_yscale("log")
    axis.set_ylabel("Positional-singleton dataset-sites (log scale)")
    axis.set_title(
        "NOT_RUN and NOT_EVALUABLE are not biological negatives",
        loc="left",
        fontweight="bold",
    )
    axis.grid(axis="y", color=COLORS["grid"], linewidth=0.7)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)
    for bar, value in zip(bars, values):
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            value * 1.15,
            f"{value:,}\n{pct(value, EXPECTED_COUNTS['singleton_sites'], 3)}",
            ha="center",
            fontsize=8.5,
        )
    figure.text(
        0.01,
        0.01,
        "M2 runs only after M1. G1, G2 and formal R1 are NOT_RUN for every "
        "positional singleton because the local multi-marker opportunity is absent.",
        color=COLORS["muted"],
        fontsize=8.5,
    )
    figure.subplots_adjust(bottom=0.17, top=0.86)
    save_figure(figure, output)


def load_case_matrix(
    assignment: Mapping[str, Any],
) -> tuple[np.ndarray, list[int], list[str]]:
    matrix_path = verify_matrix_identity(assignment)
    frame = pd.read_csv(matrix_path, na_values=["NA"])
    if "read_id" not in frame.columns:
        raise ReportBuildError(f"Methylation matrix lacks read_id: {matrix_path}")
    frame["read_id"] = frame["read_id"].astype(str)
    frame = frame.set_index("read_id")
    core_reads = assignment.get("core_reads")
    if not isinstance(core_reads, list) or not core_reads:
        raise ReportBuildError("Candidate assignment lacks core reads")
    read_ids = [str(read["read_id"]) for read in core_reads]
    labels = [str(read["label"]) for read in core_reads]
    missing = sorted(set(read_ids).difference(frame.index))
    if missing:
        raise ReportBuildError(
            f"Methylation matrix lacks candidate core reads: {missing[:3]}"
        )
    frame = frame.loc[read_ids]
    positions = [int(value) for value in frame.columns]
    matrix = frame.to_numpy(dtype=float)
    minimum_called = max(3, int(np.ceil(matrix.shape[0] * 0.20)))
    keep = np.sum(np.isfinite(matrix), axis=0) >= minimum_called
    matrix = matrix[:, keep]
    positions = [position for position, selected in zip(positions, keep) if selected]
    if matrix.shape[1] > 100:
        selection = np.unique(
            np.linspace(0, matrix.shape[1] - 1, 100).round().astype(int)
        )
        matrix = matrix[:, selection]
        positions = [positions[index] for index in selection]
    unique_labels = sorted(set(labels))
    order = sorted(
        range(len(read_ids)),
        key=lambda index: (
            unique_labels.index(labels[index]),
            np.nanmean(matrix[index])
            if np.isfinite(matrix[index]).any()
            else -1,
        ),
    )
    return matrix[order], positions, [labels[index] for index in order]


def figure_cases(
    selected: list[dict[str, Any]],
    assignments: Mapping[SiteKey, Mapping[str, Any]],
    output: Path,
) -> None:
    column_count = 2
    row_count = math.ceil(len(selected) / column_count)
    figure, axes = plt.subplots(
        row_count,
        column_count,
        figsize=(13.4, 5.1 * row_count),
        squeeze=False,
    )
    color_map = matplotlib.colormaps["cividis"].copy()
    color_map.set_bad(COLORS["missing"])
    for axis, case in zip(axes.flat, selected):
        key = (
            case["dataset"],
            case["chrom"],
            int(case["pos"]),
            case["ref"],
            case["alt"],
        )
        matrix, positions, labels = load_case_matrix(assignments[key])
        ps_values = [
            read.get("latest_ps") for read in assignments[key]["core_reads"]
        ]
        ps_missing_n = sum(value is None for value in ps_values)
        ps_distinct_n = len({value for value in ps_values if value is not None})
        unique_labels = sorted(set(labels))
        axis.imshow(
            matrix,
            aspect="auto",
            interpolation="nearest",
            cmap=color_map,
            vmin=0,
            vmax=1,
        )
        numeric = np.array([unique_labels.index(label) for label in labels])[:, None]
        group_axis = axis.inset_axes((-0.045, 0, 0.022, 1), transform=axis.transAxes)
        group_axis.imshow(
            numeric,
            aspect="auto",
            interpolation="nearest",
            cmap=ListedColormap(GROUP_COLORS[: len(unique_labels)]),
            vmin=0,
            vmax=max(1, len(unique_labels) - 1),
        )
        group_axis.set_axis_off()
        for index in range(1, len(labels)):
            if labels[index] != labels[index - 1]:
                axis.axhline(index - 0.5, color="white", linewidth=1.1)
        group_sizes = Counter(labels)
        group_text = " | ".join(
            f"G{index + 1} n={group_sizes[label]}"
            for index, label in enumerate(unique_labels)
        )
        axis.set_title(
            f"{case['dataset']}  {case['chrom']}:{case['pos']:,} "
            f"{case['ref']}>{case['alt']}\n{group_text}",
            loc="left",
            fontsize=10,
            fontweight="bold",
        )
        axis.set_ylabel(f"Focal-ALT core reads (n={matrix.shape[0]})")
        if positions:
            ticks = np.unique(
                np.linspace(0, len(positions) - 1, min(4, len(positions)))
                .round()
                .astype(int)
            )
            axis.set_xticks(
                ticks,
                [f"{positions[index] / 1e6:.3f} Mb" for index in ticks],
            )
        axis.set_xlabel("CpG genomic position")
        axis.tick_params(axis="y", left=False, labelleft=False)
        axis.spines[:].set_visible(False)
        axis.text(
            0,
            -0.17,
            f"M2 PASS; truth={case['truth_label']}; "
            f"K={case['methyl_group_count']}; min group n={case['min_group_n']}; "
            f"PS missing={ps_missing_n}/{len(ps_values)}; "
            f"nonempty PS={ps_distinct_n}",
            transform=axis.transAxes,
            fontsize=8,
            color=COLORS["muted"],
            va="top",
        )
    for axis in axes.flat[len(selected) :]:
        axis.set_visible(False)
    figure.suptitle(
        "Actual read-by-CpG matrices show residual partitions, not clone identity",
        x=0.02,
        y=0.995,
        ha="left",
        fontsize=15,
        fontweight="bold",
    )
    color_axis = figure.add_axes([0.94, 0.40, 0.009, 0.18])
    scalar = matplotlib.cm.ScalarMappable(
        norm=matplotlib.colors.Normalize(vmin=0, vmax=1),
        cmap=color_map,
    )
    color_bar = figure.colorbar(scalar, cax=color_axis)
    color_bar.set_label("Methylation probability", fontsize=8)
    figure.text(
        0.01,
        0.008,
        "Rows are grouped by frozen M1 membership. M2 excludes detected "
        "alignment with eight measured axes only; PS is not an M2 axis; "
        "G1/G2/R1 remain NOT_RUN.",
        color=COLORS["muted"],
        fontsize=8.5,
    )
    figure.subplots_adjust(
        left=0.075,
        right=0.92,
        bottom=0.055,
        top=0.94,
        hspace=0.42,
        wspace=0.25,
    )
    save_figure(figure, output)


def markdown_table(headers: list[str], rows: Iterable[Iterable[Any]]) -> str:
    result = [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join("---" for _ in headers) + "|",
    ]
    for row in rows:
        result.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(result)


def html_table(headers: list[str], rows: Iterable[Iterable[Any]]) -> str:
    heading = "".join(
        f'<th scope="col">{html.escape(str(value))}</th>' for value in headers
    )
    body_rows = []
    for row in rows:
        values = list(row)
        if not values:
            continue
        cells = (
            f'<th scope="row">{html.escape(str(values[0]))}</th>'
            + "".join(
                f"<td>{html.escape(str(value))}</td>"
                for value in values[1:]
            )
        )
        body_rows.append(f"<tr>{cells}</tr>")
    body = "".join(body_rows)
    return (
        '<div class="table-wrap"><table><thead><tr>'
        + heading
        + "</tr></thead><tbody>"
        + body
        + "</tbody></table></div>"
    )


def image_data_uri(path: Path) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def dataset_rows(summary: Mapping[str, Any]) -> list[list[str]]:
    breakdown = summary["breakdowns"]["dataset"]
    rows = []
    for dataset in EXPECTED_DATASETS:
        value = breakdown[dataset]
        rows.append(
            [
                dataset,
                f"{value['sites']:,}",
                f"{value['m1_evaluable']:,}",
                f"{value['m1_flagged']:,}",
                pct(value["m1_flagged"], value["sites"], 3),
                f"{value['m2_pass']:,}",
                pct(value["m2_pass"], value["sites"], 4),
            ]
        )
    return rows


def truth_rows(summary: Mapping[str, Any]) -> list[list[str]]:
    breakdown = summary["breakdowns"]["truth"]
    rows = []
    for truth in ("TP", "FP", "UNASSESSED"):
        value = breakdown[truth]
        rows.append(
            [
                truth,
                f"{value['sites']:,}",
                f"{value['m1_evaluable']:,}",
                f"{value['m1_flagged']:,}",
                pct(value["m1_flagged"], value["sites"], 3),
                pct(value["m1_flagged"], value["m1_evaluable"], 3),
                f"{value['m2_pass']:,}",
            ]
        )
    return rows


def candidate_rows(
    candidates: Iterable[Mapping[str, Any]],
    controls: Mapping[SiteKey, Mapping[str, Any]],
    phase_set_audit: Mapping[str, Any],
) -> list[list[str]]:
    result = []
    phase_sets = {
        (
            str(row["dataset"]),
            str(row["chrom"]),
            int(row["pos"]),
            str(row["ref"]),
            str(row["alt"]),
        ): row
        for row in phase_set_audit["sites"]
    }
    ordered = sorted(
        candidates,
        key=lambda row: (
            EXPECTED_DATASETS.index(str(row["dataset"])),
            str(row["chrom"]),
            int(row["pos"]),
        ),
    )
    for row in ordered:
        key = (
            str(row["dataset"]),
            str(row["chrom"]),
            int(row["pos"]),
            str(row["ref"]),
            str(row["alt"]),
        )
        control = controls[key]
        phase_set = phase_sets.get(key)
        if phase_set is None:
            raise ReportBuildError(f"Candidate phase-set row is missing: {key}")
        if control["ref_stable_multigroup"]:
            ref_result = "REF multigroup"
        elif control["ref_evaluable"]:
            ref_result = "REF non-replication"
        else:
            ref_result = "REF not testable"
        result.append(
            [
                str(row["dataset"]),
                f"{row['chrom']}:{int(row['pos']):,}",
                f"{row['ref']}>{row['alt']}",
                str(row["truth_label"]),
                str(row["methyl_group_count"]),
                str(row["core_read_n"]),
                str(row["min_group_n"]),
                (
                    f"{phase_set['ps_missing_n']}/"
                    f"{phase_set['core_read_n']}"
                ),
                str(phase_set["ps_distinct_nonmissing_n"]),
                str(control["n_tumor_ref"]),
                ref_result,
                "yes" if control["joint_allele_axis_aligned"] else "no",
            ]
        )
    return result


def methyl_group_rows(summary: Mapping[str, Any]) -> list[list[str]]:
    distribution = summary["methyl_group_count_distribution"]
    denominator = int(summary["counts"]["m1_flagged"])
    return [
        [str(group_count), f"{count:,}", pct(int(count), denominator, 3)]
        for group_count, count in sorted(
            ((int(key), int(value)) for key, value in distribution.items())
        )
    ]


def tumor_ref_rows(
    control_summary: Mapping[str, Mapping[str, Any]],
) -> list[list[str]]:
    labels = {
        "m1": "M1 stable multigroup",
        "m2_pass": "M2 PASS",
    }
    rows = []
    for level in ("m1", "m2_pass"):
        value = control_summary[level]
        rows.append(
            [
                labels[level],
                f"{value['sites']:,}",
                f"{value['ref_evaluable']:,}",
                f"{value['ref_nonreplication']:,}",
                f"{value['ref_stable_multigroup']:,}",
                f"{value['ref_not_evaluable']:,}",
                f"{value['joint_allele_axis_aligned']:,}",
            ]
        )
    return rows


def technical_overlap_rows(
    overlap: Mapping[str, Mapping[str, Any]],
) -> list[list[str]]:
    labels = {
        "all": "All singleton loci",
        "m1": "M1 flagged loci",
        "m2": "M2 PASS loci",
    }
    return [
        [
            labels[level],
            f"{overlap[level]['hcc1395']:,}",
            f"{overlap[level]['hcc1395_dorado']:,}",
            f"{overlap[level]['intersection']:,}",
            f"{overlap[level]['union']:,}",
            (
                f"{100 * overlap[level]['jaccard']:.2f}%"
                if overlap[level]["jaccard"] is not None
                else "NA"
            ),
        ]
        for level in ("all", "m1", "m2")
    ]


def build_markdown(
    *,
    summary: Mapping[str, Any],
    candidates: list[dict[str, Any]],
    controls: Mapping[SiteKey, Mapping[str, Any]],
    control_summary: Mapping[str, Mapping[str, Any]],
    technical_overlap: Mapping[str, Mapping[str, Any]],
    phase_set_audit: Mapping[str, Any],
    figures: Mapping[str, Path],
    inputs: Mapping[str, Mapping[str, Any]],
    repository: Mapping[str, str],
    verification: Mapping[str, int],
) -> str:
    counts = summary["counts"]
    dataset_table = markdown_table(
        [
            "dataset",
            "singleton",
            "M1 evaluable",
            "M1 multigroup",
            "M1/all",
            "M2 PASS",
            "M2/all",
        ],
        dataset_rows(summary),
    )
    truth_table = markdown_table(
        [
            "truth",
            "singleton",
            "M1 evaluable",
            "M1 multigroup",
            "M1/all",
            "M1/evaluable",
            "M2 PASS",
        ],
        truth_rows(summary),
    )
    candidate_table = markdown_table(
        [
            "dataset",
            "site",
            "change",
            "truth",
            "K",
            "core reads",
            "min group n",
            "PS missing/core",
            "nonempty PS count",
            "tumor REF reads",
            "REF control",
            "joint allele aligned",
        ],
        candidate_rows(candidates, controls, phase_set_audit),
    )
    group_table = markdown_table(
        [
            "M1 methyl groups K",
            "sites",
            f"share of {counts['m1_flagged']:,}",
        ],
        methyl_group_rows(summary),
    )
    control_table = markdown_table(
        [
            "ALT methyl level",
            "sites",
            "REF evaluable",
            "REF non-replication",
            "REF multigroup",
            "REF not testable",
            "joint allele aligned",
        ],
        tumor_ref_rows(control_summary),
    )
    technical_table = markdown_table(
        [
            "technical level",
            "HCC1395",
            "DORADO",
            "intersection",
            "union",
            "Jaccard",
        ],
        technical_overlap_rows(technical_overlap),
    )
    provenance_table = markdown_table(
        ["role", "path", "SHA-256"],
        (
            (role, record["path"], record["sha256"])
            for role, record in inputs.items()
        ),
    )
    m2_evaluable = counts["m2_pass"] + counts["m2_fail"]
    m2_wilson = wilson_interval(counts["m2_pass"], m2_evaluable)
    m2_group_distribution = Counter(
        row["methyl_group_count"] for row in candidates
    )
    m2_group_summary = "、".join(
        f"K={group_count}: {site_count}"
        for group_count, site_count in sorted(m2_group_distribution.items())
    )
    biological_sample_count = len(
        set(summary["screen_metadata"]["dataset_biological_ids"].values())
    )
    passed_audit_checks = sum(summary["checks"].values())
    total_audit_checks = len(summary["checks"])
    nearest_gap = summary["positional_recomputation"][
        "minimum_finite_singleton_nearest_gap_bp"
    ]
    co = summary["breakdowns"]["dataset"]["COLO829"]
    hcc1937 = summary["breakdowns"]["dataset"]["HCC1937"]
    m1_control = control_summary["m1"]
    m2_control = control_summary["m2_pass"]
    technical_all = technical_overlap["all"]
    technical_m1 = technical_overlap["m1"]
    ps_target = next(
        row
        for row in phase_set_audit["sites"]
        if (
            row["dataset"],
            row["chrom"],
            row["pos"],
            row["ref"],
            row["alt"],
        )
        == EXPECTED_PHASE_SET_CASE["key"]
    )
    return f"""<!--
建立時間: {datetime.now(timezone.utc).isoformat()}
目標: 回答 positional singleton sSNV 中 focal-ALT 甲基是否可分析及可觀察多群的比例
處理範圍: {len(summary['screen_metadata']['datasets'])} datasets / {biological_sample_count} biological samples / chr1-22 / {counts['singleton_sites']:,} positional-singleton dataset-sites
build_branch: {repository['branch']}
build_commit: {repository['head']}
worktree: {repository['worktree']}
data_version: {SOURCE_AUTHORITY_ID}
證據等級: L1 technical census；L2 measured-axis residual partition；L4 subclone hypothesis
驗證方式: {counts['singleton_sites']:,}-row supplemental audit {passed_audit_checks}/{total_audit_checks} checks PASS；signed final dataset/report receipts；{verification['tests']}-test canonical suite（failures={verification['failures']}，errors={verification['errors']}）
關聯檔案:
  - {summary['outputs']['site_audit']['path']}
  - {summary['outputs']['m2_pass_cases']['path']}
-->

# 單一 positional sSNV 區域的 focal-ALT 甲基多群完整驗證

> **Task Type B 全量驗證；L1 technical census，claim ceiling = L2 M2
> read-level residual epigenetic partition；cellular subclone 仍是 L4 hypothesis。**

## 開頭重點結論

1. 分母定義為同一 dataset、同一 chromosome，以相鄰 sSNV 距離 <=50,000 bp
   建立 connected component 後 `component_size=1` 的位點，共
   **{counts['singleton_sites']:,} / {summary['screen_metadata']['all_rows']:,} =
   {pct(counts['singleton_sites'], summary['screen_metadata']['all_rows'], 6)}**。
2. **{counts['m1_evaluable']:,} / {counts['singleton_sites']:,} =
   {pct(counts['m1_evaluable'], counts['singleton_sites'], 6)}**
   技術上有足夠 focal-ALT methylation matrix 可進 M1；這是「可分析率」。
3. **{counts['m1_flagged']:,} / {counts['singleton_sites']:,} =
   {pct(counts['m1_flagged'], counts['singleton_sites'], 6)}**
   在全部 singleton 中看到 M1 穩定多群；在可分析者中為
   **{pct(counts['m1_flagged'], counts['m1_evaluable'], 6)}**。這是 operational
   screen yield，不是 subclone prevalence。
4. M1 的 {counts['m1_flagged']:,} 個位點中，**{m2_evaluable:,} = {counts['m2_pass']:,}
   PASS + {counts['m2_fail']:,} FAIL** 可完成八個 measured axes 的 M2 判定；
   另有 **{counts['m2_not_evaluable']:,} NOT_EVALUABLE**。因此 M2 PASS 是
   **{counts['m2_pass']:,} / {counts['singleton_sites']:,} =
   {pct(counts['m2_pass'], counts['singleton_sites'], 6)}**。
   `{counts['m2_pass']}/{m2_evaluable}={pct(counts['m2_pass'], m2_evaluable, 1)}`
   的 Wilson 95% interval 為
   [{100 * m2_wilson[0]:.2f}%, {100 * m2_wilson[1]:.2f}%]，但只描述高度選擇
   的 conditional subset，不可外推。
5. Tumor-REF 控制把 {counts['m1_flagged']:,} 個 M1 分成
   **{m1_control['ref_nonreplication']:,} REF non-replication**、
   **{m1_control['ref_stable_multigroup']:,} REF 也有多群**與
   **{m1_control['ref_not_evaluable']:,} 不可測**。前者占全部 singleton
   **{pct(m1_control['ref_nonreplication'], counts['singleton_sites'], 6)}**，
   只增加 ALT-specific candidate 支持，不證明 subclone。
6. PS 不是 M2 的八個 measured axes 之一。30 個 M2 PASS 中有
   **{phase_set_audit['sites_with_missing_ps']} 個**至少一條 core read 缺 PS，
   共 **{phase_set_audit['total_missing_ps']}/{phase_set_audit['total_core_reads']}**
   core reads 缺 PS；因此不能宣稱已排除 phase-set confounding。
7. 所有 {counts['singleton_sites']:,} 個 positional singleton 的 G1、G2、
   formal R1 都是
   **NOT_RUN**。甲基可提出共同 ALT carrier reads 內的 latent partition，
   但目前不能確認 clone 數、cellular subclone、linear ancestry 或唯一演化樹。

![全體到M2的分母漏斗](figures/{figures['funnel'].name})

## 分母與方法

`positional singleton` 不是「read-sharing graph degree=0」，也不是保證全基因組永遠
沒有可共現的 mutation；它只表示在本次 {summary['screen_metadata']['all_rows']:,} 個最新 LongPhase-S recalibrated
`FILTER=PASS` autosomal biallelic sSNV universe 中，依 50 kb adjacency contract
沒有同 component partner；重算後最小有限 nearest gap 是 {nearest_gap:,} bp。
甲基視窗是 focal site +/-5 kb。

- M1：`coarse_ng>=2 AND not unstable AND modal_assignment_ari_min>=0.8`
- M2：M1 後，2-10 群；HP exact、HP family、strand、read start/end/length、
  MAPQ、called CpGs 八軸均可判定，且未偵測到 aligned confound。
- PS/phase-set 不在 M2 八軸內；`latest_ps` 只在本報告做描述性 census，
  不可把 HP-axis PASS 改寫成 PS confound 已排除。
- M2 PASS 的最高敘述是「測得八軸後仍存在 read-level residual epigenetic
  partition」，不是「確認 subclone」。

![M2狀態分母](figures/{figures['m2_status'].name})

### M1 多群數量

{group_table}

M1 的 {counts['m1_flagged']:,} 個 flags 中，K=2 占
{pct(summary['methyl_group_count_distribution']['2'], counts['m1_flagged'], 3)}；
M2 PASS {counts['m2_pass']} 個的群數分布為 {m2_group_summary}。兩群結構與兩個後代 clone 相容，
但 K=2 也可由兩種甲基狀態或其他二元軸產生。

## Tumor-REF 背景控制

{control_table}

在全部 M1 flags 中，{m1_control['ref_evaluable']:,} /
{m1_control['sites']:,}（{pct(m1_control['ref_evaluable'], m1_control['sites'], 3)}）
有足夠 tumor-REF reads；其中
{m1_control['ref_stable_multigroup']:,} /
{m1_control['ref_evaluable']:,}
（{pct(m1_control['ref_stable_multigroup'], m1_control['ref_evaluable'], 3)}）
在 REF 也重現多群，削弱 ALT-specific 解釋；另
{m1_control['ref_nonreplication']:,} / {m1_control['ref_evaluable']:,}
（{pct(m1_control['ref_nonreplication'], m1_control['ref_evaluable'], 3)}）
不重現，但 non-replication 不等於 subclone confirmation。

M2 PASS 的 {counts['m2_pass']} 個位點中，只有 {m2_control['ref_evaluable']} 個 REF-evaluable，
且 {m2_control['ref_nonreplication']} 個都未在 REF 重現多群；其餘 {m2_control['ref_not_evaluable']} 個無足夠
REF control。M2 + evaluable REF non-replication 的觀察下限是
**{m2_control['ref_nonreplication']} / {counts['singleton_sites']:,} =
{pct(m2_control['ref_nonreplication'], counts['singleton_sites'], 6)}**；這是優先追蹤
候選率，不是 biological prevalence。

## Dataset 結果

{dataset_table}

![各dataset結果](figures/{figures['dataset'].name})

HCC1395 與 HCC1395_DORADO 是同一 biological sample 的 technical comparison，
不可當兩個獨立樣本重現。COLO829 的 M1 multigroup 為
{co['m1_flagged']:,}/{co['sites']:,}，但 M2 PASS={co['m2_pass']}；
HCC1937 的 M2 PASS 數最多（{hcc1937['m2_pass']}），仍只屬 M2。

### HCC1395 technical-pair overlap

{technical_table}

兩個 dataset 的 singleton loci Jaccard 為
{100 * technical_all['jaccard']:.2f}%，M1 flag loci 降至
{100 * technical_m1['jaccard']:.2f}%，且各自
{technical_overlap['m2']['hcc1395']} 個 M2 PASS 沒有共同 locus。這不否定甲基
訊號，但顯示目前沒有 locus-level technical replication。

## Truth strata

{truth_table}

![truth分層](figures/{figures['truth'].name})

FP 不是本研究 universe；它只是 truth label 的一個 strata。這次分母包含 TP、FP 與
UNASSESSED 的全部 LongPhase-S PASS sSNV。FP singleton 只有 1 個 M2-evaluable
位點，因此不能由這裡估 M2 specificity。

## 真實 M2 PASS 個案

![真實read-by-CpG個案](figures/{figures['cases'].name})

熱圖列是真實 focal-ALT core reads，欄是 CpG genomic positions，顏色是 methylation
probability；左側色條只表示 frozen methyl-group membership。視覺分群與 M2 PASS
支持 read-level partition，但不提供 cellular identity。

PS 描述性重算顯示：30 個 M2 PASS 中
{phase_set_audit['sites_with_missing_ps']} 個有 missing PS，
{phase_set_audit['sites_with_multiple_nonmissing_ps']} 個跨越多個非空 PS。
最需保留 caveat 的 `HCC1937 chr5:43,849,776 T>C` 有
**{ps_target['ps_missing_n']}/{ps_target['core_read_n']}** core reads 缺 PS；
其餘 {ps_target['ps_nonmissing_n']} 條分布於
{ps_target['ps_distinct_nonmissing_n']} 個非空 PS
（{", ".join(str(value) for value in ps_target['ps_values'])}）。
這不是 M2 FAIL，但也不能宣稱 phase-set confounding 已排除。

{candidate_table}

## 對 clone 假說的解釋

使用者提出的 linear model 在邏輯上可行：若 clone 1 與 clone 2 共享較早 ancestral
ALT，單看該 ALT 的 reads 會混合兩個 clone，而後續 methyl state 可能把它們分開。
但同樣資料也可由 allele state、tumor state、局部 methylation drift、copy number、
purity、read geometry 或未測技術因素產生。M2 只排除目前八個 measured axes；
positional singleton 又缺 G1/G2 multi-marker co-segregation，所以模型目前是
**compatible hypothesis**，不是 identifiable conclusion。

升級至 cellular subclone 至少需要同一 cell population 的 single-cell、colony、
spatial 或 multi-region orthogonal evidence；升級 linear order 還需 >=3 mutation
states、CCF/order 一致並排除 branching、recurrence 與 CN 替代解釋。

## 限制

- {summary['screen_metadata']['all_rows']:,} 與 {counts['singleton_sites']:,}
  都是 dataset-sites，不是病人數或 globally unique loci。
- M1 是資料相依 operational yield，不是 population prevalence。
- {counts['m2_not_evaluable']:,} 個 M2 NOT_EVALUABLE 不能當 M2 negative；
  {counts['m2_not_run']:,} 個 NOT_RUN 更不能。
- `{counts['m2_pass']}/{m2_evaluable}` 是經 M1 與 M2 evaluability 選出的條件比例。
- matched truth strata 不足以提供 M2 specificity，尤其 FP M2-evaluable n=1。
- PS 未列入 M2 axis；missing PS 或多個非空 PS 必須保留為未排除的 phase-set
  caveat。
- Tumor-REF non-replication 只增加 ALT specificity 支持；REF 不可測不能當
  non-replication。
- bulk ONT methylation partition 不自行識別 clone number 或 ancestry topology。

## 外部方法對照

PRISM 與 EVOFLUx 類方法能由甲基推估演化，是因為使用多個 epiloci、明確混合模型、
跨時間或其他正交證據；不能把它們的結論直接移植到單一 focal ALT：

- [PRISM: methylation-based inference of subclonal evolution](https://academic.oup.com/bioinformatics/article/35/14/i520/5529252)
- [EVOFLUx: fluctuating DNA methylation traces clonal evolution](https://www.nature.com/articles/s41586-025-09374-4)

## 證據鏈

{provenance_table}
"""


def build_html(
    *,
    summary: Mapping[str, Any],
    candidates: list[dict[str, Any]],
    controls: Mapping[SiteKey, Mapping[str, Any]],
    control_summary: Mapping[str, Mapping[str, Any]],
    technical_overlap: Mapping[str, Mapping[str, Any]],
    phase_set_audit: Mapping[str, Any],
    illustrated_cases: list[dict[str, Any]],
    figures: Mapping[str, Path],
    inputs: Mapping[str, Mapping[str, Any]],
    repository: Mapping[str, str],
    verification: Mapping[str, int],
) -> str:
    counts = summary["counts"]
    dataset_table = html_table(
        [
            "Dataset",
            "Singleton",
            "M1 evaluable",
            "M1 multigroup",
            "M1 / all",
            "M2 PASS",
            "M2 / all",
        ],
        dataset_rows(summary),
    )
    truth_table = html_table(
        [
            "Truth",
            "Singleton",
            "M1 evaluable",
            "M1 multigroup",
            "M1 / all",
            "M1 / evaluable",
            "M2 PASS",
        ],
        truth_rows(summary),
    )
    candidates_table = html_table(
        [
            "Dataset",
            "Site",
            "Change",
            "Truth",
            "K",
            "Core reads",
            "Min group n",
            "PS missing / core",
            "Nonempty PS count",
            "Tumor REF reads",
            "REF control",
            "Joint allele aligned",
        ],
        candidate_rows(candidates, controls, phase_set_audit),
    )
    group_table = html_table(
        ["M1 methyl groups K", "Sites", "Share of M1 flags"],
        methyl_group_rows(summary),
    )
    control_table = html_table(
        [
            "ALT methyl level",
            "Sites",
            "REF evaluable",
            "REF non-replication",
            "REF multigroup",
            "REF not testable",
            "Joint allele aligned",
        ],
        tumor_ref_rows(control_summary),
    )
    technical_table = html_table(
        [
            "Technical level",
            "HCC1395",
            "DORADO",
            "Intersection",
            "Union",
            "Jaccard",
        ],
        technical_overlap_rows(technical_overlap),
    )
    provenance_table = html_table(
        ["Role", "Path", "SHA-256"],
        (
            (role, record["path"], record["sha256"])
            for role, record in inputs.items()
        ),
    )
    claim_table = html_table(
        ["層級", "本次狀態", "可說 / 不可說"],
        [
            (
                "M1",
                f"{counts['m1_flagged']:,} FLAGGED",
                "可說 ALT-read methyl heterogeneity；不可說 prevalence",
            ),
            (
                "M2",
                f"{counts['m2_pass']} PASS",
                "可說 measured-axis residual partition；不可說 clone",
            ),
            (
                "G1/G2/R1",
                "全部 NOT_RUN",
                "沒有 local multi-marker molecular-haplotype support",
            ),
            (
                "細胞級正交確認",
                "無資料",
                (
                    "需 single-cell/colony/spatial/multi-region "
                    "與 order evidence"
                ),
            ),
        ],
    )
    image = {name: image_data_uri(path) for name, path in figures.items()}
    m2_evaluable = counts["m2_pass"] + counts["m2_fail"]
    m2_wilson = wilson_interval(counts["m2_pass"], m2_evaluable)
    m2_group_distribution = Counter(
        row["methyl_group_count"] for row in candidates
    )
    m2_group_summary = "、".join(
        f"K={group_count}: {site_count}"
        for group_count, site_count in sorted(m2_group_distribution.items())
    )
    escaped_repository_head = html.escape(repository["head"])
    escaped_repository_branch = html.escape(repository["branch"])
    nearest_gap = summary["positional_recomputation"][
        "minimum_finite_singleton_nearest_gap_bp"
    ]
    m1_control = control_summary["m1"]
    m2_control = control_summary["m2_pass"]
    technical_all = technical_overlap["all"]
    technical_m1 = technical_overlap["m1"]
    ps_target = next(
        row
        for row in phase_set_audit["sites"]
        if (
            row["dataset"],
            row["chrom"],
            row["pos"],
            row["ref"],
            row["alt"],
        )
        == EXPECTED_PHASE_SET_CASE["key"]
    )
    return f"""<!doctype html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="color-scheme" content="light">
<title>單一 positional sSNV 區域的 focal-ALT 甲基多群完整驗證</title>
<style>
:root {{
  --ink: #17202a; --muted: #59636e; --paper: #ffffff; --soft: #eef2f3;
  --line: #d9dee3; --teal: #147d75; --blue: #2f6ea5; --amber: #d08a17;
  --red: #b84a3a; --max: 1120px;
}}
* {{ box-sizing: border-box; }}
html {{ scroll-behavior: smooth; }}
body {{
  margin: 0; color: var(--ink); background: var(--paper);
  font-family: "Noto Sans TC", "Source Han Sans TC", "PingFang TC", sans-serif;
  font-size: 16px; line-height: 1.72; letter-spacing: 0;
}}
a {{ color: var(--blue); text-decoration-thickness: 1px; }}
a:focus-visible, summary:focus-visible {{
  outline: 3px solid var(--amber); outline-offset: 3px;
}}
.skip-link {{
  position: absolute; left: 12px; top: -60px; z-index: 100;
  background: var(--paper); color: var(--ink); padding: 10px 14px;
  border: 2px solid var(--ink);
}}
.skip-link:focus {{ top: 10px; }}
.ribbon {{
  background: var(--ink); color: white; padding: 10px 24px; font-size: 13px;
  position: sticky; top: 0; z-index: 20;
}}
.ribbon strong {{ color: #78d1c8; }}
header {{ border-bottom: 1px solid var(--line); }}
.header-inner, .inner {{ width: min(var(--max), calc(100% - 40px)); margin: 0 auto; }}
.header-inner {{ padding: 52px 0 42px; }}
.eyebrow {{ color: var(--teal); font-weight: 700; font-size: 14px; margin: 0 0 8px; }}
h1 {{ font-family: "Noto Serif TC", "Source Han Serif TC", serif; font-size: clamp(31px, 5vw, 54px);
  line-height: 1.16; max-width: 900px; margin: 0 0 18px; letter-spacing: 0; }}
.lead {{ max-width: 900px; color: var(--muted); font-size: 19px; margin: 0; }}
nav {{ border-bottom: 1px solid var(--line); background: #fafbfb; }}
nav .inner {{ display: flex; gap: 22px; overflow-x: auto; padding: 12px 0; white-space: nowrap; }}
nav a {{ color: var(--ink); text-decoration: none; font-size: 14px; }}
section {{ border-bottom: 1px solid var(--line); }}
section .inner {{ padding: 42px 0 48px; }}
h2 {{ font-size: 27px; line-height: 1.25; margin: 0 0 18px; letter-spacing: 0; }}
h3 {{ font-size: 19px; margin: 30px 0 10px; letter-spacing: 0; }}
p {{ max-width: 900px; }}
.kpis {{ display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); border: 1px solid var(--line); margin: 26px 0 28px; }}
.kpi {{ padding: 20px; border-right: 1px solid var(--line); min-width: 0; }}
.kpi:last-child {{ border-right: 0; }}
.kpi b {{ display: block; font-size: 28px; color: var(--teal); line-height: 1.15; overflow-wrap: anywhere; }}
.kpi span {{ display: block; color: var(--muted); font-size: 13px; margin-top: 6px; }}
.assertion {{ border-left: 5px solid var(--teal); padding: 15px 18px; background: var(--soft); max-width: 940px; }}
.warning {{ border-left-color: var(--amber); }}
.danger {{ border-left-color: var(--red); }}
.tier {{
  display: inline-block; border: 1px solid currentColor; padding: 2px 7px;
  font-size: 12px; font-weight: 700; margin-right: 7px;
}}
.tier-l1 {{ color: var(--teal); }}
.tier-l2 {{ color: var(--blue); }}
.tier-l4 {{ color: #865900; }}
figure {{ margin: 30px 0 8px; }}
figure img {{ display: block; width: 100%; height: auto; border: 1px solid var(--line); }}
figcaption {{ color: var(--muted); font-size: 13px; margin-top: 8px; max-width: 920px; }}
.table-wrap {{ width: 100%; overflow-x: auto; border: 1px solid var(--line); margin: 20px 0; }}
table {{ border-collapse: collapse; width: 100%; min-width: 760px; font-size: 13px; }}
th, td {{ padding: 10px 12px; border-bottom: 1px solid var(--line); text-align: left; vertical-align: top; }}
th {{ background: var(--soft); font-weight: 700; position: sticky; top: 0; }}
td:nth-child(n+2) {{ font-variant-numeric: tabular-nums; }}
code {{ font-family: "IBM Plex Mono", "SFMono-Regular", monospace; font-size: .9em; background: var(--soft); padding: 1px 4px; border-radius: 3px; }}
details {{ border-top: 1px solid var(--line); padding: 14px 0; }}
summary {{ cursor: pointer; font-weight: 700; }}
footer {{ background: var(--ink); color: #dce2e5; }}
footer .inner {{ padding: 34px 0; font-size: 13px; overflow-wrap: anywhere; }}
@media (max-width: 760px) {{
  .header-inner, .inner {{ width: min(100% - 28px, var(--max)); }}
  .header-inner {{ padding: 34px 0 30px; }}
  .ribbon {{ padding: 9px 14px; position: static; }}
  .lead {{ font-size: 17px; }}
  section .inner {{ padding: 32px 0 36px; }}
  .kpis {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
  .kpi:nth-child(2) {{ border-right: 0; }}
  .kpi:nth-child(-n+2) {{ border-bottom: 1px solid var(--line); }}
}}
@media print {{
  .ribbon, nav, .skip-link {{ display: none; }}
  details > * {{ display: block; }}
  body {{ color: #000; background: #fff; }}
  section {{ break-inside: avoid; }}
}}
</style>
</head>
<body>
<a class="skip-link" href="#main-content">跳至主要內容</a>
<div class="ribbon"><strong>Task Type B</strong> · 全 {counts['singleton_sites']:,} positional singleton · claim ceiling M2 · 不是 clone confirmation</div>
<header>
  <div class="header-inner">
    <p class="eyebrow">InterSubMod · Latest LongPhase-S PASS sSNV</p>
    <h1>單一 positional sSNV 區域的 focal-ALT 甲基多群完整驗證</h1>
    <p class="lead">回答「無法做局部 sSNV 共現的位點中，有多少仍可用甲基觀察 ALT-carrier reads 內多群」，並把可分析、M1 多群、M2 residual evidence 與 subclone claim 清楚分開。</p>
  </div>
</header>
<nav><div class="inner">
  <a href="#answer">直答</a><a href="#denominator">分母</a><a href="#tumor-ref">Tumor-REF</a><a href="#dataset">Dataset</a>
  <a href="#truth">Truth</a><a href="#cases">真實個案</a><a href="#clone">Clone 解釋</a><a href="#provenance">證據鏈</a>
</div></nav>
<main id="main-content">
<section id="answer"><div class="inner">
  <h2>直答：{pct(counts['m1_evaluable'], counts['singleton_sites'], 2)} 可分析，{pct(counts['m1_flagged'], counts['singleton_sites'], 2)} 見 M1 多群，{pct(counts['m2_pass'], counts['singleton_sites'], 4)} 通過 M2</h2>
  <p><span class="tier tier-l1">L1 technical census</span><span class="tier tier-l2">L2 measured-axis residual</span><span class="tier tier-l4">L4 subclone hypothesis</span></p>
  <div class="kpis">
    <div class="kpi"><b>{pct(counts['m1_evaluable'], counts['singleton_sites'], 2)}</b><span>{counts['m1_evaluable']:,} / {counts['singleton_sites']:,} 技術上 M1-evaluable</span></div>
    <div class="kpi"><b>{pct(counts['m1_flagged'], counts['singleton_sites'], 2)}</b><span>{counts['m1_flagged']:,} / {counts['singleton_sites']:,} M1 stable multigroup</span></div>
    <div class="kpi"><b>{pct(counts['m2_pass'], counts['singleton_sites'], 4)}</b><span>{counts['m2_pass']:,} / {counts['singleton_sites']:,} M2 PASS</span></div>
    <div class="kpi"><b>0 tested</b><span>G1、G2、formal R1 structural NOT_RUN；不是陰性</span></div>
  </div>
  <p class="assertion"><strong>可確認：</strong>甲基在大部分 singleton 位點技術上可分析，且 {pct(counts['m1_flagged'], counts['singleton_sites'], 2)} 有穩定 ALT-read 多群 screen signal。<br><strong>不可確認：</strong>多群就是 cellular subclone，或兩群一定是 linear clone 1 → clone 2。</p>
  <p>M2 的 {counts['m2_pass']} PASS / {m2_evaluable} evaluable =
  {pct(counts['m2_pass'], m2_evaluable, 1)}，Wilson 95% interval
  [{100 * m2_wilson[0]:.2f}%, {100 * m2_wilson[1]:.2f}%]；這只適用高度選擇的
  conditional subset。{counts['m2_not_evaluable']:,} 個 M1 flags 因 measured-axis
  power/判定不足而 NOT_EVALUABLE，不能當作陰性。</p>
  <p class="assertion warning"><strong>PS caveat：</strong>PS 不在 M2 八個 measured axes 內。30 個 M2 PASS 中有 {phase_set_audit['sites_with_missing_ps']} 個至少一條 core read 缺 PS，合計 {phase_set_audit['total_missing_ps']} / {phase_set_audit['total_core_reads']} core reads；HP-axis PASS 不等於 PS confounding 已排除。</p>
  <figure><img src="{image['funnel']}" alt="全體到 M2 的分母漏斗"><figcaption>橫軸為 dataset-site 數量（log scale）。M2 的 {m2_evaluable} 個 evaluable 位點是經 M1 與八軸可判定條件選出的子集。</figcaption></figure>
</div></section>
<section id="denominator"><div class="inner">
  <h2>分母：50 kb positional contract，不是 read-sharing degree-zero</h2>
  <p>同一 dataset、同一 chromosome，以相鄰 sSNV 距離 ≤50,000 bp 建立 connected component，component size=1 定義為 positional singleton。共有 <strong>{counts['singleton_sites']:,} / {summary['screen_metadata']['all_rows']:,} = {pct(counts['singleton_sites'], summary['screen_metadata']['all_rows'], 6)}</strong>；最小有限 nearest gap 為 {nearest_gap:,} bp。</p>
  <p>這不表示全基因組永遠沒有任何可共現 mutation，也不等於 read-sharing graph degree=0。它只說本次最新 LongPhase-S PASS、chr1-22 biallelic sSNV universe 的 50 kb local contract 沒有 partner。甲基使用 focal site ±5 kb 的 read-by-CpG matrix。</p>
  <h3>M1 群數分布</h3>
  {group_table}
  <p>M1 flags 中 K=2 占 {pct(summary['methyl_group_count_distribution']['2'], counts['m1_flagged'], 3)}；M2 PASS 的 {counts['m2_pass']} 個位點群數分布為 {m2_group_summary}。兩群與兩個後代 clone 相容，但也可由兩種 methyl states 或其他二元因素產生。</p>
  <figure><img src="{image['m2_status']}" alt="M2 PASS FAIL NOT EVALUABLE NOT RUN 分母"><figcaption>{counts['m2_not_run']:,} 是 M1 未標記而 M2 NOT_RUN；{counts['m2_not_evaluable']:,} 是 M1 後八軸仍不可完整判定。兩者均不可寫成 biological negative。</figcaption></figure>
</div></section>
<section id="tumor-ref"><div class="inner">
  <h2>Tumor-REF 控制：區分 ALT-specific candidate 與局部背景多群</h2>
  {control_table}
  <p>在 {m1_control['sites']:,} 個 M1 flags 中，
  {m1_control['ref_evaluable']:,}（{pct(m1_control['ref_evaluable'], m1_control['sites'], 3)}）有足夠 tumor-REF reads。
  其中 {m1_control['ref_stable_multigroup']:,} / {m1_control['ref_evaluable']:,}
  （{pct(m1_control['ref_stable_multigroup'], m1_control['ref_evaluable'], 3)}）
  在 REF 也重現多群，削弱 ALT specificity；另
  {m1_control['ref_nonreplication']:,} / {m1_control['ref_evaluable']:,}
  （{pct(m1_control['ref_nonreplication'], m1_control['ref_evaluable'], 3)}）
  不重現，但這只支持 ALT-specific prioritization。</p>
  <p class="assertion warning">M2 PASS {counts['m2_pass']} 個位點只有 {m2_control['ref_evaluable']} 個 REF-evaluable，其中 {m2_control['ref_nonreplication']} 個未重現 REF 多群；{m2_control['ref_not_evaluable']} 個無足夠 REF control。M2 + REF non-replication 為 {m2_control['ref_nonreplication']} / {counts['singleton_sites']:,} = {pct(m2_control['ref_nonreplication'], counts['singleton_sites'], 6)}，仍不是 subclone confirmation。</p>
</div></section>
<section id="dataset"><div class="inner">
  <h2>Dataset 分層</h2>
  {dataset_table}
  <figure><img src="{image['dataset']}" alt="各 dataset M1 與 M2 比例"><figcaption>左圖是 M1 operational yield；右圖是全部 singleton 中的 M2 PASS。HCC1395 與 DORADO 只算一個 biological sample。</figcaption></figure>
  <h3>HCC1395 technical-pair overlap</h3>
  {technical_table}
  <p>Singleton loci Jaccard 為 {100 * technical_all['jaccard']:.2f}%，M1 flag loci 降為 {100 * technical_m1['jaccard']:.2f}%；兩邊各有 {technical_overlap['m2']['hcc1395']} 個 M2 PASS，但共同 PASS locus 為 {technical_overlap['m2']['intersection']}。這是 technical comparison，並未形成 locus-level replication。</p>
</div></section>
<section id="truth"><div class="inner">
  <h2>Truth strata：不是只有 FP，也不能由此估 specificity</h2>
  <p>研究 universe 包含 TP、FP、UNASSESSED 的全部 LongPhase-S PASS sSNV；FP 只是其中一個 truth strata。</p>
  {truth_table}
  <figure><img src="{image['truth']}" alt="TP FP UNASSESSED M1 比例"><figcaption>FP 的 M1/all 較低，但 FP M2-evaluable 僅 n=1，因此不能把 M2 結果解讀為已證明 clone-specific 或 caller-specific。</figcaption></figure>
</div></section>
<section id="cases"><div class="inner">
  <h2>真實 M2 PASS read-by-CpG 個案</h2>
  <p class="assertion warning">熱圖顯示真實 ALT-carrier read methylation partition；左側群標是 frozen M1 assignment。M2 PASS 只表示八個 measured axes 下仍有 residual partition。</p>
  <figure><img src="{image['cases']}" alt="{len(illustrated_cases)} 個實際 M2 PASS methylation heatmap"><figcaption>每個有 M2 PASS 的 dataset 保留 core-read 數最多的一例，另固定呈現 DORADO chr1/chr22 與 HCC1937 PS-caveat 個案；COLO829 沒有 M2 PASS。圖中不得把群組命名為 clone 1 / clone 2。</figcaption></figure>
  <p class="assertion warning"><strong>Phase-set 個案：</strong>30 個 M2 PASS 中 {phase_set_audit['sites_with_missing_ps']} 個有 missing PS，{phase_set_audit['sites_with_multiple_nonmissing_ps']} 個跨越多個非空 PS。HCC1937 chr5:43,849,776 T&gt;C 有 {ps_target['ps_missing_n']} / {ps_target['core_read_n']} core reads 缺 PS；其餘 {ps_target['ps_nonmissing_n']} 條分布於 {ps_target['ps_distinct_nonmissing_n']} 個非空 PS（{", ".join(str(value) for value in ps_target['ps_values'])}）。此位點仍是 M2 PASS，但不能宣稱 phase-set confounding 已排除。</p>
  <details open><summary>全部 {counts['m2_pass']} 個 M2 PASS 位點</summary>{candidates_table}</details>
</div></section>
<section id="clone"><div class="inner">
  <h2>對 ancestral ALT 與 linear subclone 假說的結論</h2>
  <p>假說在邏輯上成立：若 clone 1 與 clone 2 共享較早 ancestral ALT，單看 focal ALT reads 會混合兩個 clone，而 methyl state 可能再把它們分開。這正是甲基最有價值的候選提出功能。</p>
  <p class="assertion danger"><strong>但目前不可識別：</strong>同一 heatmap 也可能來自 allele/tumor state、局部 methylation drift、CN、purity、read geometry 或未測因素。positional singleton 缺 G1/G2 multi-marker co-segregation，因此只能稱 compatible latent substructure hypothesis。</p>
  <p><a href="https://academic.oup.com/bioinformatics/article/35/14/i520/5529252">PRISM</a> 與 <a href="https://www.nature.com/articles/s41586-025-09374-4">EVOFLUx</a> 類研究依賴多 epiloci、混合模型、跨時間或正交證據；不能把其演化推論直接移植到單一 focal ALT。</p>
  {claim_table}
</div></section>
<section id="provenance"><div class="inner">
  <h2>證據鏈與限制</h2>
  <ul>
    <li>{summary['screen_metadata']['all_rows']:,} 與 {counts['singleton_sites']:,} 是 dataset-sites，不是獨立病人或 globally unique loci。</li>
    <li>M1 是資料相依 screen yield；{counts['m2_pass']}/{m2_evaluable} 是高度選擇的條件比例。</li>
    <li>M2 僅處理八個 measured axes，不排除所有 technical/biological/CN confounding。</li>
    <li>PS 未列入 M2 axis；missing PS 或多個非空 PS 是未排除的 phase-set caveat。</li>
    <li>Tumor-REF non-replication 只增加 ALT specificity 支持；REF 不可測不能當作 non-replication。</li>
    <li>bulk ONT methylation partition 不自行識別 clone number、linear order 或唯一演化樹。</li>
  </ul>
  <details><summary>輸入路徑與 SHA-256</summary>{provenance_table}</details>
</div></section>
</main>
<footer><div class="inner">Schema {SCHEMA_NAME} v{SCHEMA_VERSION} · Source-authorized commit {escaped_repository_head} · Branch {escaped_repository_branch} · Canonical tests {verification['tests']} / failures {verification['failures']} / errors {verification['errors']} · Generated {datetime.now(timezone.utc).isoformat()} · Offline portable HTML · No external assets</div></footer>
</body>
</html>
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-summary", type=Path, required=True)
    parser.add_argument("--site-audit", type=Path, required=True)
    parser.add_argument("--candidate-cases", type=Path, required=True)
    parser.add_argument("--stable-assignments", type=Path, required=True)
    parser.add_argument("--tumor-ref-controls", type=Path, required=True)
    parser.add_argument("--canonical-pytest-xml", type=Path, required=True)
    parser.add_argument("--final-dataset-receipt", type=Path, required=True)
    parser.add_argument("--final-dataset-signature", type=Path, required=True)
    parser.add_argument("--final-dataset-public-key", type=Path, required=True)
    parser.add_argument("--final-report-receipt", type=Path, required=True)
    parser.add_argument("--final-report-signature", type=Path, required=True)
    parser.add_argument("--final-report-public-key", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def validate_pre_publish_state(
    *,
    input_paths: Mapping[str, Path],
    summary: Mapping[str, Any],
    identities_expected: Mapping[str, Mapping[str, Any]],
    source_expected: Mapping[str, Any],
    repository_expected: Mapping[str, str],
    verification_expected: Mapping[str, Any],
    dataset_verification_expected: Mapping[str, Any],
    report_verification_expected: Mapping[str, Any],
) -> None:
    identities_observed = {
        name: artifact(path) for name, path in input_paths.items()
    }
    source_observed = artifact(Path(__file__).resolve())
    repository_observed = repository_metadata(summary)
    validate_source_authority_chain(summary)
    verification_observed = validate_canonical_pytest_release_evidence(
        input_paths["canonical_pytest_xml"]
    )
    dataset_observed, report_observed = validate_formal_release_chain(
        final_dataset_receipt_path=input_paths["final_dataset_receipt"],
        final_dataset_signature_path=input_paths["final_dataset_signature"],
        final_dataset_public_key_path=input_paths["final_dataset_public_key"],
        final_report_receipt_path=input_paths["final_report_receipt"],
        final_report_signature_path=input_paths["final_report_signature"],
        final_report_public_key_path=input_paths["final_report_public_key"],
    )
    if (
        identities_observed != identities_expected
        or source_observed != source_expected
        or repository_observed != repository_expected
        or verification_observed != verification_expected
        or dataset_observed != dataset_verification_expected
        or report_observed != report_verification_expected
    ):
        raise ReportBuildError(
            "Input, source, test, or formal release identity changed "
            "before atomic publish"
        )


def validate_and_publish(
    *,
    staging: Path,
    output_dir: Path,
    input_paths: Mapping[str, Path],
    summary: Mapping[str, Any],
    identities_expected: Mapping[str, Mapping[str, Any]],
    source_expected: Mapping[str, Any],
    repository_expected: Mapping[str, str],
    verification_expected: Mapping[str, Any],
    dataset_verification_expected: Mapping[str, Any],
    report_verification_expected: Mapping[str, Any],
) -> None:
    validate_pre_publish_state(
        input_paths=input_paths,
        summary=summary,
        identities_expected=identities_expected,
        source_expected=source_expected,
        repository_expected=repository_expected,
        verification_expected=verification_expected,
        dataset_verification_expected=dataset_verification_expected,
        report_verification_expected=report_verification_expected,
    )
    rename_no_replace(staging, output_dir)


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.expanduser().resolve()
    if os.path.lexists(output_dir):
        raise FileExistsError(f"Refusing to overwrite output directory: {output_dir}")
    input_paths = {
        "audit_summary": args.audit_summary.expanduser().resolve(strict=True),
        "site_audit": args.site_audit.expanduser().resolve(strict=True),
        "candidate_cases": args.candidate_cases.expanduser().resolve(strict=True),
        "stable_assignments": args.stable_assignments.expanduser().resolve(
            strict=True
        ),
        "tumor_ref_controls": args.tumor_ref_controls.expanduser().resolve(
            strict=True
        ),
        "canonical_pytest_xml": args.canonical_pytest_xml.expanduser().resolve(
            strict=True
        ),
        "final_dataset_receipt": args.final_dataset_receipt.expanduser().resolve(
            strict=True
        ),
        "final_dataset_signature": args.final_dataset_signature.expanduser().resolve(
            strict=True
        ),
        "final_dataset_public_key": args.final_dataset_public_key.expanduser().resolve(
            strict=True
        ),
        "final_report_receipt": args.final_report_receipt.expanduser().resolve(
            strict=True
        ),
        "final_report_signature": args.final_report_signature.expanduser().resolve(
            strict=True
        ),
        "final_report_public_key": args.final_report_public_key.expanduser().resolve(
            strict=True
        ),
        "test_evidence_manifest": TEST_EVIDENCE_MANIFEST.resolve(strict=True),
        "test_evidence_signature": TEST_EVIDENCE_SIGNATURE.resolve(strict=True),
        "test_evidence_public_key": TEST_EVIDENCE_PUBLIC_KEY.resolve(strict=True),
    }
    summary_preview = load_json(
        input_paths["audit_summary"], label="singleton audit summary preview"
    )
    input_paths.update(resolve_source_authority_input_paths(summary_preview))
    identities_before = {
        name: artifact(path) for name, path in input_paths.items()
    }
    source_before = artifact(Path(__file__).resolve())
    verification = validate_canonical_pytest_release_evidence(
        input_paths["canonical_pytest_xml"]
    )
    (
        summary,
        dataset_receipt,
        report_receipt,
        builder_receipt,
        dataset_verification,
        report_verification,
    ) = validate_inputs(
        summary_path=input_paths["audit_summary"],
        site_audit_path=input_paths["site_audit"],
        candidate_path=input_paths["candidate_cases"],
        assignments_path=input_paths["stable_assignments"],
        final_dataset_receipt_path=input_paths["final_dataset_receipt"],
        final_dataset_signature_path=input_paths["final_dataset_signature"],
        final_dataset_public_key_path=input_paths["final_dataset_public_key"],
        final_report_receipt_path=input_paths["final_report_receipt"],
        final_report_signature_path=input_paths["final_report_signature"],
        final_report_public_key_path=input_paths["final_report_public_key"],
        tumor_ref_controls_path=input_paths["tumor_ref_controls"],
    )
    candidates = load_candidates(input_paths["candidate_cases"])
    m1_keys = load_m1_singleton_keys(input_paths["site_audit"])
    candidate_keys = {
        (
            row["dataset"],
            row["chrom"],
            int(row["pos"]),
            row["ref"],
            row["alt"],
        )
        for row in candidates
    }
    controls, control_summary = load_tumor_ref_controls(
        input_paths["tumor_ref_controls"],
        m1_keys=m1_keys,
        m2_keys=candidate_keys,
    )
    technical_overlap = technical_pair_overlap(input_paths["site_audit"])
    assignments = load_candidate_assignments(
        input_paths["stable_assignments"], candidate_keys
    )
    phase_set_audit = summarize_phase_sets(candidates, assignments)
    validate_phase_set_audit(phase_set_audit)
    selected = select_cases(candidates)
    repository = repository_metadata(summary)

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = output_dir.parent / (
        f".{output_dir.name}.staging.{os.getpid()}.{uuid.uuid4().hex}"
    )
    staging.mkdir(mode=0o755)
    figure_dir = staging / "figures"
    figure_dir.mkdir(mode=0o755)
    setup_style()
    figures = {
        "funnel": figure_dir / "01_singleton_evidence_funnel.png",
        "dataset": figure_dir / "02_dataset_m1_m2_yield.png",
        "truth": figure_dir / "03_truth_strata_m1_yield.png",
        "m2_status": figure_dir / "04_m2_status_denominators.png",
        "cases": figure_dir / "05_actual_m2_pass_methylation_heatmaps.png",
    }
    figure_funnel(summary, figures["funnel"])
    figure_dataset(summary, figures["dataset"])
    figure_truth(summary, figures["truth"])
    figure_m2_status(summary, figures["m2_status"])
    figure_cases(selected, assignments, figures["cases"])

    repository_after = repository_metadata(summary)
    identities_after = {
        name: artifact(path) for name, path in input_paths.items()
    }
    source_after = artifact(Path(__file__).resolve())
    if (
        identities_before != identities_after
        or source_before != source_after
        or repository_after != repository
    ):
        raise ReportBuildError("Input or report source identity changed during build")

    markdown_name = "20260718_單位點sSNV_ALT內甲基多群完整驗證_01.md"
    html_name = "20260718_單位點sSNV_ALT內甲基多群完整驗證_01.standalone.html"
    markdown_path = staging / markdown_name
    html_path = staging / html_name
    markdown_path.write_text(
        build_markdown(
            summary=summary,
            candidates=candidates,
            controls=controls,
            control_summary=control_summary,
            technical_overlap=technical_overlap,
            phase_set_audit=phase_set_audit,
            figures=figures,
            inputs=identities_after,
            repository=repository,
            verification=verification,
        ),
        encoding="utf-8",
    )
    html_path.write_text(
        build_html(
            summary=summary,
            candidates=candidates,
            controls=controls,
            control_summary=control_summary,
            technical_overlap=technical_overlap,
            phase_set_audit=phase_set_audit,
            illustrated_cases=selected,
            figures=figures,
            inputs=identities_after,
            repository=repository,
            verification=verification,
        ),
        encoding="utf-8",
    )
    for path in (markdown_path, html_path):
        path.chmod(0o444)
        fsync_file(path)
    figure_dir.chmod(0o555)

    final_markdown = output_dir / markdown_name
    final_html = output_dir / html_name
    receipt = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_type": "B_comprehensive_validation",
        "scope": "all_50432_positional_singleton_dataset_sites",
        "command": [str(Path(sys.executable).resolve()), *sys.argv],
        "inputs": identities_after,
        "code": {"report_builder": source_after},
        "formal_release_links": {
            "final_dataset_receipt_pass": dataset_receipt["pass"],
            "final_report_receipt_pass": report_receipt["pass"],
            "final_dataset_builder_receipt_pass": builder_receipt["pass"],
        },
        "formal_release_verifications": {
            "final_dataset": dataset_verification,
            "final_report": report_verification,
        },
        "selected_visual_cases": [
            {
                "dataset": row["dataset"],
                "chrom": row["chrom"],
                "pos": row["pos"],
                "ref": row["ref"],
                "alt": row["alt"],
                "truth_label": row["truth_label"],
                "core_read_n": row["core_read_n"],
                "tumor_ref_control": controls[
                    (
                        row["dataset"],
                        row["chrom"],
                        row["pos"],
                        row["ref"],
                        row["alt"],
                    )
                ],
            }
            for row in selected
        ],
        "technical_pair_overlap": technical_overlap,
        "phase_set_audit": phase_set_audit,
        "tumor_ref_control_summary": control_summary,
        "repository": repository,
        "canonical_test_summary": verification,
        "outputs": {
            "markdown": artifact(markdown_path, public_path=final_markdown),
            "portable_html": artifact(html_path, public_path=final_html),
            "figures": {
                name: artifact(
                    path,
                    public_path=output_dir / "figures" / path.name,
                )
                for name, path in figures.items()
            },
        },
        "claim_ceiling": "M2_read_level_residual_epigenetic_partition",
        "pass_semantics": (
            "complete_reader_facing_singleton_report_not_cellular_subclone_"
            "or_lineage_confirmation"
        ),
        "pass": True,
    }
    receipt_path = staging / "report_build_receipt.json"
    final_receipt = output_dir / receipt_path.name
    receipt_path.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    receipt_path.chmod(0o444)
    fsync_file(receipt_path)
    success_path = staging / "_SUCCESS.json"
    final_success = output_dir / success_path.name
    success_path.write_text(
        json.dumps(
            {
                "schema_name": "intersubmod.atomic_release_marker",
                "schema_version": "1.0.0",
                "receipt": artifact(
                    receipt_path, public_path=final_receipt
                ),
                "pass": True,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    success_path.chmod(0o444)
    fsync_file(success_path)
    fsync_directory(staging)
    staging.chmod(0o555)
    validate_and_publish(
        staging=staging,
        output_dir=output_dir,
        input_paths=input_paths,
        summary=summary,
        identities_expected=identities_before,
        source_expected=source_before,
        repository_expected=repository,
        verification_expected=verification,
        dataset_verification_expected=dataset_verification,
        report_verification_expected=report_verification,
    )
    fsync_directory(output_dir.parent)
    print(
        json.dumps(
            {
                "output_dir": str(output_dir),
                "markdown": str(final_markdown),
                "portable_html": str(final_html),
                "receipt": str(final_receipt),
                "success_marker": str(final_success),
                "m2_pass_cases": len(candidates),
                "illustrated_cases": len(selected),
                "m1_tumor_ref_evaluable": control_summary["m1"][
                    "ref_evaluable"
                ],
                "m2_pass_tumor_ref_evaluable": control_summary["m2_pass"][
                    "ref_evaluable"
                ],
                "canonical_tests": verification["tests"],
                "pass": True,
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
