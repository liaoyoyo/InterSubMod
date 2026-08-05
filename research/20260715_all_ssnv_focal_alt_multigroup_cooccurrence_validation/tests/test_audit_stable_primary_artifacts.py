import ast
import csv
import gzip
import importlib.util
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).parents[1]
    / "scripts"
    / "audit_stable_primary_artifacts.py"
)
SPEC = importlib.util.spec_from_file_location("audit_stable_primary_artifacts", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)

SCREEN_PRODUCER = (
    Path(__file__).parents[1]
    / "scripts"
    / "analyze_all_ssnv_focal_alt_multigroup.py"
)
CANONICAL_SCREEN_CONTRACT = "phylo-v4.1_column_null95_modal_K10_RNULL40_min_group3"


def write_site_table(path: Path, *, stable: bool = True) -> None:
    fields = [
        "sample",
        "chrom",
        "pos",
        "ref",
        "alt",
        "stable_null_multigroup",
        "screen_contract",
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "sample": "S1",
                "chrom": "chr1",
                "pos": 100,
                "ref": "A",
                "alt": "C",
                "stable_null_multigroup": str(stable).lower(),
                "screen_contract": MODULE.SCREEN_CONTRACT,
            }
        )


def make_fixture(tmp_path: Path) -> tuple[Path, Path, Path, dict]:
    region = tmp_path / "S1" / "chr1_100_A_C" / "region_1"
    paths = {
        "reads": region / "reads" / "reads.tsv",
        "distance_matrix": region / "distance" / "BERNOULLI" / "matrix.csv",
        "methylation_matrix": region / "methylation" / "methylation.csv",
    }
    for role, path in paths.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"{role}\n", encoding="utf-8")
    site_results = tmp_path / "sites.tsv.gz"
    write_site_table(site_results)
    assignment = {
        "schema_name": MODULE.ASSIGNMENT_SCHEMA,
        "schema_version": "1.0.0",
        "screen_contract": MODULE.SCREEN_CONTRACT,
        "artifact_identity_contract": "sha256_size_path_v1",
        "strict_confirm_candidate": True,
        "sample": "S1",
        "chrom": "chr1",
        "pos": 100,
        "region_dir": str(region.resolve()),
        "posthoc": {"ref": "A", "alt": "C"},
        "primary_artifacts": {
            role: {
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "sha256": MODULE.sha256(path),
            }
            for role, path in paths.items()
        },
    }
    assignments = tmp_path / "assignments.jsonl.gz"
    with gzip.open(assignments, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps(assignment) + "\n")
    return site_results, assignments, tmp_path / "audit.json", paths


def run(
    site_results: Path,
    assignments: Path,
    output: Path,
    *,
    workers: int = 1,
    consumer_receipts: tuple[Path, ...] = (),
) -> dict:
    return MODULE.run_audit(
        site_results,
        assignments,
        output,
        consumer_receipts=consumer_receipts,
        workers=workers,
        chunk_size=1,
        max_pending_chunks=1,
        progress_every=1,
        command=["synthetic-test"],
    )


def source_string_constant(path: Path, name: str) -> str:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    values = [
        node.value.value
        for node in tree.body
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == name for target in node.targets)
        and isinstance(node.value, ast.Constant)
        and isinstance(node.value.value, str)
    ]
    assert len(values) == 1
    return values[0]


def test_screen_contract_matches_canonical_producer_source() -> None:
    producer_contract = source_string_constant(SCREEN_PRODUCER, "SCREEN_CONTRACT")
    assert producer_contract == CANONICAL_SCREEN_CONTRACT
    assert MODULE.SCREEN_CONTRACT == producer_contract


def test_full_primary_artifact_identity_audit(tmp_path: Path) -> None:
    site_results, assignments, output, _ = make_fixture(tmp_path)
    receipt = run(site_results, assignments, output)
    assert receipt["counts"] == {
        "stable_sites": 1,
        "assignment_records": 1,
        "primary_artifacts_expected": 3,
        "primary_artifacts_verified": 3,
    }
    assert len(receipt["verification"]["artifact_set_sha256"]) == 64
    assert receipt["execution"] == {
        "workers": 1,
        "chunk_size": 1,
        "max_pending_chunks": 1,
        "progress_every": 1,
    }
    assert receipt["pass"] is True
    original = output.read_bytes()
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        run(site_results, assignments, output)
    assert output.read_bytes() == original


def test_tampered_primary_artifact_fails_before_receipt(tmp_path: Path) -> None:
    site_results, assignments, output, paths = make_fixture(tmp_path)
    original_size = paths["reads"].stat().st_size
    paths["reads"].write_text("xxxxx\n", encoding="utf-8")
    assert paths["reads"].stat().st_size == original_size
    with pytest.raises(MODULE.AuditContractError, match="reads sha256 mismatch"):
        run(site_results, assignments, output)
    assert not output.exists()


def test_artifact_set_digest_is_deterministic_across_input_order() -> None:
    first = {
        "key": ["S1", "chr1", 100, "A", "C"],
        "artifacts": [
            {"role": "reads", "path": "/tmp/a", "size_bytes": 1, "sha256": "a" * 64},
            {"role": "distance_matrix", "path": "/tmp/b", "size_bytes": 2, "sha256": "b" * 64},
        ],
    }
    second = {
        "key": ["S2", "chr2", 200, "G", "T"],
        "artifacts": [
            {"role": "methylation_matrix", "path": "/tmp/c", "size_bytes": 3, "sha256": "c" * 64},
        ],
    }
    expected = MODULE.artifact_set_digest([first, second])
    assert MODULE.artifact_set_digest([second, first]) == expected
    assert MODULE.artifact_set_digest(
        [
            {**first, "artifacts": list(reversed(first["artifacts"]))},
            second,
        ]
    ) == expected


def test_workers_two_propagates_worker_contract_error(tmp_path: Path) -> None:
    site_results, assignments, output, paths = make_fixture(tmp_path)
    paths["reads"].write_text("xxxxx\n", encoding="utf-8")
    with pytest.raises(MODULE.AuditContractError, match="reads sha256 mismatch"):
        run(site_results, assignments, output, workers=2)
    assert not output.exists()


def test_consumer_receipts_and_timestamp_shape(tmp_path: Path) -> None:
    site_results, assignments, output, _ = make_fixture(tmp_path)
    consumers = (tmp_path / "cooccurrence_receipt.json", tmp_path / "tumor_ref_receipt.json")
    for index, path in enumerate(consumers):
        path.write_text(json.dumps({"receipt": index}) + "\n", encoding="utf-8")

    receipt = run(
        site_results,
        assignments,
        output,
        consumer_receipts=consumers,
    )

    assert receipt["created_at_utc"] == receipt["finished_at_utc"]
    assert datetime.fromisoformat(receipt["started_at_utc"]) <= datetime.fromisoformat(
        receipt["finished_at_utc"]
    )
    records = receipt["inputs"]["consumer_receipts"]
    assert [record["path"] for record in records] == [str(path.resolve()) for path in consumers]
    assert [record["sha256"] for record in records] == [MODULE.sha256(path) for path in consumers]
    assert all({"path", "size_bytes", "mtime_ns", "sha256"} == set(record) for record in records)


def test_input_artifact_change_during_verification_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    site_results, assignments, output, _ = make_fixture(tmp_path)
    consumer = tmp_path / "consumer.json"
    consumer.write_text('{"state":1}\n', encoding="utf-8")
    original_verify = MODULE.verify_task

    def mutating_verify(task: dict) -> dict:
        result = original_verify(task)
        consumer.write_text('{"state":2}\n', encoding="utf-8")
        return result

    monkeypatch.setattr(MODULE, "verify_task", mutating_verify)
    with pytest.raises(MODULE.AuditContractError, match="changed during verification"):
        run(
            site_results,
            assignments,
            output,
            consumer_receipts=(consumer,),
        )
    assert not output.exists()


def test_cli_accepts_repeatable_consumer_receipts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    consumers = [tmp_path / "first.json", tmp_path / "second.json"]
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT),
            "--site-results",
            "sites.tsv.gz",
            "--stable-assignments",
            "assignments.jsonl.gz",
            "--consumer-receipt",
            str(consumers[0]),
            "--consumer-receipt",
            str(consumers[1]),
            "--output",
            "audit.json",
        ],
    )
    assert MODULE.parse_args().consumer_receipt == consumers


def test_canonical_command_records_isolated_python(tmp_path: Path) -> None:
    command = MODULE.canonical_command(
        site_results=tmp_path / "sites.tsv.gz",
        assignments=tmp_path / "assignments.jsonl.gz",
        consumer_receipts=(),
        output=tmp_path / "audit.json",
        workers=40,
        chunk_size=64,
        max_pending_chunks=80,
        progress_every=1000,
    )
    assert command[:4] == MODULE.canonical_python_prefix()
    assert command[4] == str(SCRIPT.resolve())


def test_dangling_symlink_is_treated_as_existing_output(tmp_path: Path) -> None:
    site_results, assignments, output, _ = make_fixture(tmp_path)
    output.symlink_to(tmp_path / "missing-target.json")
    assert os.path.lexists(output)
    assert not output.exists()
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        run(site_results, assignments, output)
    assert output.is_symlink()


def test_atomic_create_does_not_overwrite_racing_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    site_results, assignments, output, _ = make_fixture(tmp_path)
    original_loader = MODULE.load_stable_sites

    def racing_loader(path: Path) -> set[MODULE.SiteKey]:
        output.write_text("created-by-another-process\n", encoding="utf-8")
        return original_loader(path)

    monkeypatch.setattr(MODULE, "load_stable_sites", racing_loader)
    with pytest.raises(FileExistsError):
        run(site_results, assignments, output)
    assert output.read_text(encoding="utf-8") == "created-by-another-process\n"


def test_stable_site_assignment_key_set_must_be_exact(tmp_path: Path) -> None:
    site_results, assignments, output, _ = make_fixture(tmp_path)
    write_site_table(site_results, stable=False)
    with pytest.raises(MODULE.AuditContractError, match="key-set mismatch"):
        run(site_results, assignments, output)
    assert not output.exists()
