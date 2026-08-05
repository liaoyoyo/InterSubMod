from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
SCRIPT = SCRIPTS / "merge_all_ssnv_screen_recovery.py"
SPEC = importlib.util.spec_from_file_location("merge_all_ssnv_screen_recovery", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_sites(path: Path, samples: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, MODULE.A.SITE_FIELDS, delimiter="\t")
        writer.writeheader()
        for index, sample in enumerate(samples, 1):
            row = {field: "" for field in MODULE.A.SITE_FIELDS}
            row.update(
                {
                    "dataset": sample,
                    "sample": sample,
                    "chrom": "chr1",
                    "pos": str(index),
                    "ref": "A",
                    "alt": "C",
                }
            )
            writer.writerow(row)


def write_assignments(path: Path, samples: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for index, sample in enumerate(samples, 1):
            handle.write(
                json.dumps(
                    {
                        "sample": sample,
                        "chrom": "chr1",
                        "pos": index,
                        "posthoc": {"ref": "A", "alt": "C"},
                    }
                )
                + "\n"
            )


def test_prefix_iterators_stop_before_replacement(tmp_path: Path) -> None:
    sites = tmp_path / "sites.tsv.gz"
    assignments = tmp_path / "assignments.jsonl.gz"
    samples = ["HCC1937", "HCC1954", "HCC1954"]
    write_sites(sites, samples)
    write_assignments(assignments, samples)

    observed_sites = list(MODULE.iter_prefix_sites(sites, "HCC1954"))
    observed_assignments = list(MODULE.iter_prefix_assignments(assignments, "HCC1954"))

    assert [row["sample"] for row in observed_sites] == ["HCC1937"]
    assert [row["sample"] for row in observed_assignments] == ["HCC1937"]


def test_csv_false_values_remain_false_for_summary() -> None:
    row = {field: "False" for field in MODULE.BOOL_FIELDS}
    normalized = MODULE.summary_row(row)
    assert all(normalized[field] is False for field in MODULE.BOOL_FIELDS)
    row["stable_null_multigroup"] = "not-a-bool"
    with pytest.raises(RuntimeError, match="Invalid boolean"):
        MODULE.summary_row(row)


def test_assignment_key_requires_complete_posthoc_identity() -> None:
    row = {
        "sample": "HCC1954",
        "chrom": "chr17",
        "pos": 100,
        "posthoc": {"ref": "a", "alt": "c"},
    }
    assert MODULE.assignment_key(row) == ("HCC1954", "chr17", 100, "A", "C")
    with pytest.raises(RuntimeError, match="missing posthoc"):
        MODULE.assignment_key({"sample": "HCC1954", "chrom": "chr17", "pos": 100})


def test_replacement_receipt_requires_recorded_seed_parallelism(tmp_path: Path) -> None:
    replacement = tmp_path / "replacement"
    replacement.mkdir()
    sites = replacement / "all_ssnv_site_results.tsv.gz"
    assignments = replacement / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
    summary_path = replacement / "all_ssnv_summary.json"
    receipt_path = replacement / "run_manifest.json"
    manifest = tmp_path / "manifest.json"
    manifest.write_text("{}\n", encoding="utf-8")
    analyzer = tmp_path / "pinned_analyzer.py"
    analyzer.write_text("# pinned\n", encoding="utf-8")
    analyzer.chmod(0o444)
    started = datetime.now(timezone.utc) + timedelta(seconds=1)
    finished = started + timedelta(seconds=1)
    write_sites(sites, ["HCC1954"])
    write_assignments(assignments, [])
    summary = {
        "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_screen",
        "pass": True,
        "scope": {
            "full_469849": False,
            "selected_samples": ["HCC1954"],
            "processed_sites": 1,
        },
    }
    summary_path.write_text(json.dumps(summary) + "\n", encoding="utf-8")

    def write_receipt(workers: int) -> None:
        receipt = {
            "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_run_manifest",
            "pass": True,
            "started_at_utc": started.isoformat(),
            "finished_at_utc": finished.isoformat(),
            "input_manifest": MODULE.A.artifact(manifest),
            "source_code": {
                role: MODULE.A.artifact(analyzer)
                for role in (
                    "analyzer",
                    "focal_alt_cluster_lib",
                    "latest_tag_join",
                    "claim_contract_v2",
                )
            },
            "execution": {
                "selected_samples": ["HCC1954"],
                "phylo_parallel_workers": workers,
            },
            "outputs": {
                "site_results": MODULE.A.artifact(sites),
                "stable_assignments": MODULE.A.artifact(assignments),
                "summary": MODULE.A.artifact(summary_path),
            },
        }
        receipt_path.write_text(json.dumps(receipt) + "\n", encoding="utf-8")

    write_receipt(1)
    with pytest.raises(RuntimeError, match="did not record seed-level parallel"):
        MODULE.validate_replacement(replacement, manifest, "HCC1954", 1)
    write_receipt(11)
    _, _, observed_summary, observed_receipt = MODULE.validate_replacement(
        replacement, manifest, "HCC1954", 1
    )
    assert observed_summary["pass"] is True
    assert observed_receipt["execution"]["phylo_parallel_workers"] == 11


def test_replacement_requires_read_only_source_older_than_start(tmp_path: Path) -> None:
    replacement = tmp_path / "replacement"
    replacement.mkdir()
    sites = replacement / "all_ssnv_site_results.tsv.gz"
    assignments = replacement / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
    summary_path = replacement / "all_ssnv_summary.json"
    receipt_path = replacement / "run_manifest.json"
    manifest = tmp_path / "manifest.json"
    analyzer = tmp_path / "analyzer.py"
    manifest.write_text("{}\n", encoding="utf-8")
    analyzer.write_text("# mutable\n", encoding="utf-8")
    write_sites(sites, ["HCC1954"])
    write_assignments(assignments, [])
    summary_path.write_text(
        json.dumps(
            {
                "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_screen",
                "pass": True,
                "scope": {
                    "full_469849": False,
                    "selected_samples": ["HCC1954"],
                    "processed_sites": 1,
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    started = datetime.now(timezone.utc) + timedelta(seconds=1)
    receipt_path.write_text(
        json.dumps(
            {
                "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_run_manifest",
                "pass": True,
                "started_at_utc": started.isoformat(),
                "finished_at_utc": (started + timedelta(seconds=1)).isoformat(),
                "input_manifest": MODULE.A.artifact(manifest),
                "source_code": {
                    role: MODULE.A.artifact(analyzer)
                    for role in (
                        "analyzer",
                        "focal_alt_cluster_lib",
                        "latest_tag_join",
                        "claim_contract_v2",
                    )
                },
                "execution": {
                    "selected_samples": ["HCC1954"],
                    "phylo_parallel_workers": 11,
                },
                "outputs": {
                    "site_results": MODULE.A.artifact(sites),
                    "stable_assignments": MODULE.A.artifact(assignments),
                    "summary": MODULE.A.artifact(summary_path),
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(RuntimeError, match="not a read-only pinned source"):
        MODULE.validate_replacement(replacement, manifest, "HCC1954", 1)


def test_prefix_requires_exact_source_lock_child_receipt(tmp_path: Path) -> None:
    prefix = tmp_path / "prefix"
    prefix.mkdir()
    sites = prefix / "all_ssnv_site_results.tsv.gz"
    assignments = prefix / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
    summary_path = prefix / "all_ssnv_summary.json"
    receipt_path = prefix / "run_manifest.json"
    lock_path = prefix / "source_lock_receipt.json"
    manifest = tmp_path / "manifest.json"
    manifest.write_text("{}\n", encoding="utf-8")
    samples = MODULE.A.DATASETS[:-1]
    write_sites(sites, samples)
    write_assignments(assignments, [])
    summary_path.write_text(
        json.dumps(
            {
                "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_screen",
                "pass": True,
                "scope": {
                    "full_469849": False,
                    "selected_samples": samples,
                    "processed_sites": len(samples),
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    receipt = {
        "schema_name": "intersubmod.all_ssnv_focal_alt_multigroup_run_manifest",
        "pass": True,
        "input_manifest": MODULE.A.artifact(manifest),
        "execution": {"selected_samples": samples, "phylo_parallel_workers": 1},
        "outputs": {
            "site_results": MODULE.A.artifact(sites),
            "stable_assignments": MODULE.A.artifact(assignments),
            "summary": MODULE.A.artifact(summary_path),
        },
    }
    receipt_path.write_text(json.dumps(receipt) + "\n", encoding="utf-8")
    identity = {"/pinned.py": {"sha256": "a" * 64}}
    lock = {
        "schema_name": MODULE.SOURCE_LOCK_SCHEMA,
        "pass": True,
        "scope": samples,
        "child_run_manifest": MODULE.A.artifact(receipt_path),
        "source_identity_before": identity,
        "source_identity_after": identity,
        "checks": {"all_locked_sources_unchanged": True},
    }
    lock_path.write_text(json.dumps(lock) + "\n", encoding="utf-8")
    observed = MODULE.validate_prefix(
        prefix, manifest, samples, len(samples)
    )
    assert observed[2]["pass"] is True
    lock["child_run_manifest"]["sha256"] = "0" * 64
    lock_path.write_text(json.dumps(lock) + "\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="child receipt mismatch"):
        MODULE.validate_prefix(prefix, manifest, samples, len(samples))


def test_equivalence_receipt_requires_all_machine_gates(tmp_path: Path) -> None:
    analyzer = tmp_path / "pinned.py"
    analyzer.write_text("# pinned\n", encoding="utf-8")
    receipt_path = tmp_path / "equivalence.json"
    checks = {
        "pinned_analyzer_sha256_exact": True,
        "source_identity_unchanged": True,
        "synthetic_full_payload_exact": True,
        "real_nested_full_payload_exact": True,
        "real_nested_parallel_triggered": True,
    }
    payload = {
        "schema_name": MODULE.EQUIVALENCE_SCHEMA,
        "pass": True,
        "scope": "algorithm_and_real_nested_high_read_fixture",
        "checks": checks,
        "real_fixture": {
            "full_row_and_assignment_payload_exact": True,
            "n_alt_after_peel": 250,
            "parallel_min_reads": 200,
        },
        "inputs": {
            "source_identity_before": {"analyzer": MODULE.A.artifact(analyzer)}
        },
    }
    receipt_path.write_text(json.dumps(payload) + "\n", encoding="utf-8")
    assert MODULE.validate_equivalence(receipt_path)["pass"] is True
    payload["checks"].pop("real_nested_full_payload_exact")
    receipt_path.write_text(json.dumps(payload) + "\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="checks are incomplete"):
        MODULE.validate_equivalence(receipt_path)
