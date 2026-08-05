from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = TOPIC_ROOT.parents[1]
RUNNER_PATH = TOPIC_ROOT / "scripts" / "run_exact_ps_cpp_topology_all7.py"
PARTITION_BINARY = REPO_ROOT / "build" / "bin" / "exact_ps_partition"

SPEC = importlib.util.spec_from_file_location("all7_pipeline_runner", RUNNER_PATH)
assert SPEC is not None and SPEC.loader is not None
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


SITE_FIELDS = ("site_index", "chrom", "pos1", "ref", "alt")
MEMBERSHIP_FIELDS = (
    "dataset",
    "chrom",
    "linkage_basis",
    "phase_set",
    "phase_set_status",
    "inference_role",
    "threshold",
    "site_index",
    "pos1",
    "component_id",
)
CALL_FIELDS = (
    "dataset",
    "chrom",
    "molecule_id",
    "hp_family",
    "phase_set",
    "site_indices",
    "positions1",
    "call_codes",
)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def write_sidecar(path: Path) -> None:
    path.with_name(f"{path.name}.sha256").write_text(
        f"{sha256_path(path)}  {path.name}\n",
        encoding="ascii",
    )


def write_tsv_gz(
    path: Path,
    fields: tuple[str, ...],
    rows: list[dict[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def build_upstream(tmp_path: Path) -> tuple[Path, Path]:
    sample = "HCC1395"
    chrom = "chr1"
    extraction_root = tmp_path / "source" / "extraction_root"
    strict_root = tmp_path / "source" / "strict_root"
    extraction = extraction_root / "chromosomes" / chrom / "extraction"
    strict = strict_root / "chromosomes" / chrom
    extraction.mkdir(parents=True)
    strict.mkdir(parents=True)

    site_catalog = extraction / f"{sample}.{chrom}.site_catalog.tsv.gz"
    extraction_membership = (
        extraction / f"{sample}.{chrom}.site_component_membership.tsv.gz"
    )
    strict_membership = (
        strict / f"{sample}.{chrom}.site_component_membership.tsv.gz"
    )
    molecule_calls = extraction / f"{sample}.{chrom}.molecule_sparse_calls.tsv.gz"
    site_rows = [
        {"site_index": 0, "chrom": chrom, "pos1": 100, "ref": "C", "alt": "T"},
        {"site_index": 1, "chrom": chrom, "pos1": 200, "ref": "G", "alt": "A"},
    ]
    membership_rows = [
        {
            "dataset": sample,
            "chrom": chrom,
            "linkage_basis": "PS_HP1",
            "phase_set": "100",
            "phase_set_status": "KNOWN_PS_PRIMARY",
            "inference_role": "PRIMARY_PS_AWARE",
            "threshold": "3",
            "site_index": index,
            "pos1": pos1,
            "component_id": "CP1",
        }
        for index, pos1 in ((0, 100), (1, 200))
    ]
    molecule_rows = [
        {
            "dataset": sample,
            "chrom": chrom,
            "molecule_id": f"molecule_{index}",
            "hp_family": "1",
            "phase_set": "100",
            "site_indices": "0,1",
            "positions1": "100,200",
            "call_codes": "AA",
        }
        for index in range(3)
    ]
    write_tsv_gz(site_catalog, SITE_FIELDS, site_rows)
    write_tsv_gz(extraction_membership, MEMBERSHIP_FIELDS, membership_rows)
    write_tsv_gz(strict_membership, MEMBERSHIP_FIELDS, membership_rows)
    write_tsv_gz(molecule_calls, CALL_FIELDS, molecule_rows)

    extraction_receipt = extraction / "receipt.json"
    extraction_receipt.write_text(
        json.dumps(
            {
                "schema_name": RUNNER.EXTRACTION_SCHEMA,
                "schema_version": RUNNER.EXTRACTION_VERSION,
                "all_pass": True,
                "scope": {"dataset": sample, "chrom": chrom, "n_sSNV": 2},
                "checks": {"synthetic_upstream_valid": True},
                "outputs": {
                    path.name: file_identity(path)
                    for path in (
                        site_catalog,
                        extraction_membership,
                        molecule_calls,
                    )
                },
            }
        ),
        encoding="utf-8",
    )
    write_sidecar(extraction_receipt)

    strict_receipt = strict / "receipt.json"
    strict_receipt.write_text(
        json.dumps(
            {
                "schema_name": RUNNER.STRICT_SCHEMA,
                "schema_version": RUNNER.STRICT_VERSION,
                "all_pass": True,
                "scope": {"dataset": sample, "chrom": chrom, "candidate_loci_S": 2},
                "checks": {"synthetic_strict_valid": True},
                "outputs": {"membership": file_identity(strict_membership)},
            }
        ),
        encoding="utf-8",
    )
    write_sidecar(strict_receipt)
    return extraction_root, strict_root


def write_manifest(
    path: Path,
    extraction_root: Path,
    strict_root: Path,
) -> None:
    samples = {
        sample: {
            "extraction_root": str(extraction_root),
            "strict_root": str(strict_root),
        }
        for sample in RUNNER.CANONICAL_DATASETS
    }
    path.write_text(
        json.dumps(
            {
                "schema_name": RUNNER.MANIFEST_SCHEMA,
                "schema_version": RUNNER.MANIFEST_VERSION,
                "dataset_order": list(RUNNER.CANONICAL_DATASETS),
                "samples": samples,
                "source_report": "synthetic",
            }
        ),
        encoding="utf-8",
    )


def run_partial(
    tmp_path: Path,
    *,
    resume: bool = False,
) -> subprocess.CompletedProcess[str]:
    extraction_root, strict_root = build_upstream(tmp_path) if not resume else (
        tmp_path / "source" / "extraction_root",
        tmp_path / "source" / "strict_root",
    )
    manifest = tmp_path / "manifest.json"
    if not resume:
        write_manifest(manifest, extraction_root, strict_root)
    command = [
        sys.executable,
        str(RUNNER_PATH),
        "--manifest",
        str(manifest),
        "--output-root",
        str(tmp_path / "output"),
        "--samples",
        "HCC1395",
        "--chroms",
        "chr1",
        "--partition-binary",
        str(PARTITION_BINARY),
    ]
    if resume:
        command.append("--resume")
    return subprocess.run(
        command,
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )


@pytest.mark.skipif(
    not PARTITION_BINARY.is_file(),
    reason="build/bin/exact_ps_partition is required for the integration test",
)
def test_synthetic_one_chromosome_partition_pipeline_passes(tmp_path: Path) -> None:
    completed = run_partial(tmp_path)
    assert completed.returncode == 0, (
        f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    )
    sample_root = tmp_path / "output" / "samples" / "HCC1395"
    partial = json.loads(
        (sample_root / "partial_partition_receipt.json").read_text()
    )
    chrom = json.loads(
        (
            sample_root
            / "chromosomes"
            / "chr1"
            / "stage_receipt.json"
        ).read_text()
    )
    comparison = json.loads(
        (
            sample_root
            / "chromosomes"
            / "chr1"
            / "comparison.json"
        ).read_text()
    )
    assert partial["all_pass"] is True
    assert partial["scope"]["autosomes_complete"] is False
    assert not (sample_root / "run_receipt.json").exists()
    assert not (sample_root / "HCC1395.exact_ps_mlhp.json").exists()
    assert chrom["checks"]["partition_uses_strict_membership"] is True
    assert comparison["all_pass"] is True
    assert comparison["mismatch_count"] == 0


@pytest.mark.parametrize(
    "argv",
    (
        ("--output-root", "/tmp/unused", "--samples", "HCC1395,HCC1395"),
        ("--output-root", "/tmp/unused", "--samples", "COLO829,HCC1395"),
        ("--output-root", "/tmp/unused", "--chroms", "chr2,chr1"),
        ("--output-root", "/tmp/unused", "--chroms", "chr23"),
        ("--output-root", "/tmp/unused", "--max-family-size", "-1"),
    ),
)
def test_argument_contract_rejects_invalid_scope(
    argv: tuple[str, ...],
) -> None:
    with pytest.raises(SystemExit) as exc:
        RUNNER.parse_args(argv)
    assert exc.value.code == 2


@pytest.mark.skipif(
    not PARTITION_BINARY.is_file(),
    reason="build/bin/exact_ps_partition is required for the integration test",
)
def test_resume_rejects_tampered_bound_artifact(tmp_path: Path) -> None:
    first = run_partial(tmp_path)
    assert first.returncode == 0, first.stderr
    comparison = (
        tmp_path
        / "output"
        / "samples"
        / "HCC1395"
        / "chromosomes"
        / "chr1"
        / "comparison.json"
    )
    with comparison.open("a", encoding="utf-8") as handle:
        handle.write("\n")

    resumed = run_partial(tmp_path, resume=True)
    assert resumed.returncode == 2
    assert "identity mismatch" in resumed.stderr


def test_topology_receipt_requires_separate_completeness_flags(
    tmp_path: Path,
) -> None:
    mlhp = tmp_path / "input.json"
    output = tmp_path / "topology.jsonl"
    receipt = tmp_path / "topology.receipt.json"
    mlhp.write_text('{"groups":[]}\n', encoding="utf-8")
    output.write_text('{"region_id":"R1","status":"ABSTAIN"}\n', encoding="utf-8")
    receipt.write_text(
        json.dumps(
            {
                "all_pass": True,
                "input": file_identity(mlhp),
                "output": file_identity(output),
                "all_mutation_bearing_families_complete": False,
                "all_units_objective_certified": False,
                "parameters": {
                    "max_family_size": 100_000,
                    "max_search_nodes": 100_000_000,
                },
                "counts": {"total_units": 1},
                "status_census": {
                    "unit_status": {"ABSTAIN": 1},
                    "objective_status": {"ABSTAIN": 1},
                    "family_status": {"ABSTAIN": 1},
                    "read_af_status": {"not_evaluable": 1},
                },
            }
        ),
        encoding="utf-8",
    )
    observed = RUNNER.validate_topology_receipt(
        receipt,
        input_path=mlhp,
        output_path=output,
        max_family_size=100_000,
        max_search_nodes=100_000_000,
    )
    assert observed["all_pass"] is True
    assert observed["all_mutation_bearing_families_complete"] is False
    assert observed["all_units_objective_certified"] is False
    assert observed["_verified_jsonl_rows"] == 1


def test_stage_retry_preserves_failed_attempt_logs(tmp_path: Path) -> None:
    logs = tmp_path / "logs"
    with pytest.raises(RUNNER.RunnerError):
        RUNNER.run_command(
            stage="adapter",
            command=[sys.executable, "-c", "raise SystemExit(3)"],
            logs_dir=logs,
        )
    first_stderr = logs / "adapter.stderr.log"
    first_receipt = logs / "adapter.stage.json"
    first_stderr_sha = sha256_path(first_stderr)
    assert json.loads(first_receipt.read_text())["attempt"] == 1

    retry = RUNNER.run_command(
        stage="adapter",
        command=[sys.executable, "-c", "print('retry-pass')"],
        logs_dir=logs,
    )
    assert retry["attempt"] == 2
    assert (logs / "adapter.attempt2.stdout.log").read_text().strip() == "retry-pass"
    assert sha256_path(first_stderr) == first_stderr_sha
