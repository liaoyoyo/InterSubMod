import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    REPO_ROOT
    / "research/20260718_k_gt8_read_supported_segmentation"
    / "scripts/run_hcc1395_cached_span_sensitivity.py"
)
DATASET = "HCC1395"


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _identity(path):
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }


def _write_json(path, value, *, sidecar=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if sidecar:
        path.with_name(f"{path.name}.sha256").write_text(
            f"{_sha256(path)}  {path.name}\n",
            encoding="ascii",
        )


def _canonical_sha256(value):
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def _write_tsv_gz(path, fieldnames, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _build_source_test_run(tmp_path, chrom, positions, calls):
    source = tmp_path / f"source_{chrom}"
    chrom_root = source / "chromosomes" / chrom
    extraction_dir = chrom_root / "extract"
    extraction_dir.mkdir(parents=True)
    extractor = tmp_path / f"dummy_extractor_{chrom}.py"
    extractor.write_text("# synthetic source authority\n", encoding="utf-8")
    manifest = tmp_path / f"synthetic_manifest_{chrom}.json"
    manifest.write_text('{"synthetic":true}\n', encoding="utf-8")

    site_catalog = extraction_dir / f"{DATASET}.{chrom}.site_catalog.tsv.gz"
    molecule_calls = (
        extraction_dir / f"{DATASET}.{chrom}.molecule_sparse_calls.tsv.gz"
    )
    _write_tsv_gz(
        site_catalog,
        ("site_index", "chrom", "pos1", "ref", "alt"),
        (
            {
                "site_index": index,
                "chrom": chrom,
                "pos1": position,
                "ref": "C",
                "alt": "T",
            }
            for index, position in enumerate(positions)
        ),
    )
    _write_tsv_gz(
        molecule_calls,
        (
            "dataset",
            "chrom",
            "molecule_id",
            "hp_family",
            "phase_set",
            "site_indices",
            "call_codes",
        ),
        calls,
    )

    components = []
    for position in positions:
        if not components or position - components[-1][-1] > 50_000:
            components.append([position])
        else:
            components[-1].append(position)
    targets = [component for component in components if len(component) > 8]
    inventory = {
        "chrom": chrom,
        "ssnv_sites": len(positions),
        "positional_components_all": len(components),
        "k_gt8_components": len(targets),
        "k_gt8_sites": sum(len(component) for component in targets),
        "k_gt8_max_k": max((len(component) for component in targets), default=0),
        "positions_sha256": _canonical_sha256(list(positions)),
    }
    contract = {
        "schema_name": (
            "intersubmod.hcc1395_full_k_gt8_segmentation.run_contract"
        ),
        "schema_version": "1.0.0",
        "sample": DATASET,
        "scope": {
            "chromosomes": [chrom],
            "comprehensive_chr1_22": False,
            "test_mode": True,
        },
        "parameters": {
            "legacy_gap_bp": 50_000,
            "max_block_size": 8,
            "mapq_min": 20,
            "baseq_min": 20,
            "bridge_thresholds": "1,2,3,5",
            "samtools_threads": 1,
        },
        "frozen_manifest": _identity(manifest),
        "tools": {
            "extractor": _identity(extractor),
            "python": _identity(Path(sys.executable)),
            "time": _identity(Path("/usr/bin/time")),
        },
        "vcf_derived_inventory": [inventory],
    }
    contract_path = source / "run_contract.json"
    _write_json(contract_path, contract)
    contract_sha = _sha256(contract_path)

    child = {
        "schema_name": (
            "intersubmod.lossless_read_linkage_chromosome_receipt"
        ),
        "schema_version": "1.3.0",
        "scope": {
            "dataset": DATASET,
            "chrom": chrom,
            "n_sSNV": len(positions),
        },
        "provenance": {
            "extractor": {
                "path": str(extractor.resolve()),
                "sha256": _sha256(extractor),
            },
            "manifest": {
                "path": str(manifest.resolve()),
                "sha256": _sha256(manifest),
            },
        },
        "outputs": {
            site_catalog.name: _identity(site_catalog),
            molecule_calls.name: _identity(molecule_calls),
        },
        "checks": {"synthetic_contract": True},
        "all_pass": True,
    }
    child_path = extraction_dir / "receipt.json"
    _write_json(child_path, child)

    logs = {}
    for name, filename, payload in (
        ("stdout", "extract.stdout.log", b"synthetic\n"),
        ("stderr", "extract.stderr.log", b""),
        ("resource", "extract.time_v.txt", b"synthetic time\n"),
    ):
        path = chrom_root / filename
        path.write_bytes(payload)
        logs[name] = _identity(path)
    stage = {
        "schema_name": (
            "intersubmod.hcc1395_full_k_gt8_segmentation.stage_receipt"
        ),
        "schema_version": "1.0.0",
        "all_pass": True,
        "sample": DATASET,
        "chrom": chrom,
        "stage": "extract",
        "status": "COMPLETED",
        "exit_code": 0,
        "contract_sha256": contract_sha,
        "tool": _identity(extractor),
        "logs": logs,
        "child_receipt": _identity(child_path),
    }
    extraction_command = [
        str(Path(sys.executable).resolve()),
        str(extractor.resolve()),
        "--manifest",
        str(manifest.resolve()),
        "--sample",
        DATASET,
        "--chrom",
        chrom,
        "--output-dir",
        str(extraction_dir.resolve()),
        "--mapq-min",
        "20",
        "--baseq-min",
        "20",
        "--bridge-thresholds",
        "1,2,3,5",
        "--samtools-threads",
        "1",
    ]
    stage["command"] = extraction_command
    stage["timed_command"] = [
        str(Path("/usr/bin/time").resolve()),
        "-v",
        "-o",
        str((chrom_root / "extract.time_v.txt").resolve()),
        "--",
        *extraction_command,
    ]
    stage_path = chrom_root / "extract.stage_receipt.json"
    _write_json(stage_path, stage)

    source_summary = source / "chromosome_summary.tsv"
    source_summary.write_text(
        "chrom\tstatus\n"
        f"{chrom}\tCOMPLETED\n",
        encoding="utf-8",
    )
    receipt = {
        "schema_name": (
            "intersubmod.hcc1395_full_k_gt8_segmentation.run_receipt"
        ),
        "schema_version": "1.0.0",
        "all_pass": True,
        "comprehensive_all_pass": False,
        "sample": DATASET,
        "scope": {
            "chromosomes": [chrom],
            "test_mode": True,
        },
        "contract": _identity(contract_path),
        "totals": {
            "chromosomes": 1,
            "ssnv_sites": len(positions),
            "k_gt8_components": inventory["k_gt8_components"],
            "k_gt8_sites": inventory["k_gt8_sites"],
            "k_gt8_max_k": inventory["k_gt8_max_k"],
            "partitioned_chromosomes": int(bool(targets)),
            "zero_target_skipped_chromosomes": int(not targets),
        },
        "chromosomes": [
            {
                "chrom": chrom,
                "inventory": inventory,
                "extraction": _identity(stage_path),
                "partition": _identity(stage_path),
                "partition_status": (
                    "COMPLETED" if targets else "SKIP_NO_K_GT8_TARGET"
                ),
            }
        ],
        "outputs": {
            "chromosome_summary": _identity(source_summary),
        },
    }
    receipt_path = source / "receipt.json"
    _write_json(receipt_path, receipt)
    marker = {
        "schema_name": (
            "intersubmod.hcc1395_full_k_gt8_segmentation.success_marker"
        ),
        "schema_version": "1.0.0",
        "all_pass": True,
        "comprehensive": False,
        "sample": DATASET,
        "scope": [chrom],
        "run_receipt": {
            "path": str(receipt_path.resolve()),
            "sha256": _sha256(receipt_path),
        },
    }
    _write_json(source / "_TEST_SUCCESS", marker, sidecar=False)
    return source, molecule_calls


def _run(source, output, chrom, cap, *, resume=False):
    command = [
        sys.executable,
        str(SCRIPT),
        "--source-full-root",
        str(source),
        "--output-root",
        str(output),
        "--span-caps",
        str(cap),
        "--chromosomes",
        chrom,
        "--test-mode",
    ]
    if resume:
        command.append("--resume")
    return subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )


def test_chr22_single_cap_cached_run_and_resume(tmp_path):
    positions = (
        100,
        200,
        300,
        400,
        10_000,
        10_100,
        10_200,
        10_300,
        10_400,
    )
    calls = [
        {
            "dataset": DATASET,
            "chrom": "chr22",
            "molecule_id": molecule_id,
            "hp_family": "1",
            "phase_set": "100",
            "site_indices": ",".join(map(str, indices)),
            "call_codes": "A" * len(indices),
        }
        for molecule_id, indices in (
            ("left", (0, 1, 2, 3)),
            ("bridge", (3, 4)),
            ("right", (4, 5, 6, 7, 8)),
        )
    ]
    source, _ = _build_source_test_run(
        tmp_path, "chr22", positions, calls
    )
    output = tmp_path / "span_output"

    completed = _run(source, output, "chr22", 500)
    assert completed.returncode == 0, (
        f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    )
    receipt = json.loads((output / "receipt.json").read_text())
    assert receipt["all_pass"] is True
    assert receipt["comprehensive_all_pass"] is False
    assert receipt["scope"]["span_caps_bp"] == [500]
    assert receipt["totals"]["completed_partition_tasks"] == 1
    assert receipt["totals"]["zero_target_skipped_tasks"] == 0
    assert receipt["timing_scope"].startswith(
        "cached_partition_wall_seconds excludes"
    )
    cap_root = output / "caps/span_500bp"
    stage_path = (
        cap_root / "chromosomes/chr22/partition.stage_receipt.json"
    )
    time_path = cap_root / "chromosomes/chr22/partition.time_v.txt"
    child_path = cap_root / "chromosomes/chr22/partition/receipt.json"
    assert stage_path.is_file()
    assert stage_path.with_name(f"{stage_path.name}.sha256").is_file()
    assert time_path.is_file()
    child = json.loads(child_path.read_text())
    assert child["parameters"]["max_block_span_bp"] == 500
    assert child["counts"]["unavoidable_span_cap_patterns"] == 1
    assert (output / "summary.tsv.sha256").is_file()
    assert (output / "run_contract.json.sha256").is_file()
    assert (output / "receipt.json.sha256").is_file()
    assert (output / "_TEST_SUCCESS.sha256").is_file()

    before = {
        "stage": _sha256(stage_path),
        "time": _sha256(time_path),
        "child": _sha256(child_path),
    }
    resumed = _run(source, output, "chr22", 500, resume=True)
    assert resumed.returncode == 0, (
        f"stdout:\n{resumed.stdout}\nstderr:\n{resumed.stderr}"
    )
    assert before == {
        "stage": _sha256(stage_path),
        "time": _sha256(time_path),
        "child": _sha256(child_path),
    }

    no_resume = _run(source, output, "chr22", 500)
    assert no_resume.returncode == 2
    assert "output root exists" in no_resume.stderr


def test_chr21_zero_target_writes_skip_and_time_file(tmp_path):
    source, _ = _build_source_test_run(
        tmp_path,
        "chr21",
        (100,),
        (),
    )
    output = tmp_path / "skip_output"
    completed = _run(source, output, "chr21", 50_000)
    assert completed.returncode == 0, (
        f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    )
    chrom_root = (
        output / "caps/span_50000bp/chromosomes/chr21"
    )
    skip_path = chrom_root / "partition.skip_receipt.json"
    skip = json.loads(skip_path.read_text())
    assert skip["status"] == "SKIP_NO_K_GT8_TARGET"
    assert skip["wall_seconds"] == 0.0
    assert (chrom_root / "partition.time_v.txt").is_file()
    assert not (chrom_root / "partition").exists()
    receipt = json.loads((output / "receipt.json").read_text())
    assert receipt["totals"]["completed_partition_tasks"] == 0
    assert receipt["totals"]["zero_target_skipped_tasks"] == 1
    skip_sha = _sha256(skip_path)
    resumed = _run(
        source,
        output,
        "chr21",
        50_000,
        resume=True,
    )
    assert resumed.returncode == 0, resumed.stderr
    assert _sha256(skip_path) == skip_sha


def test_source_output_identity_drift_fails_before_output_creation(tmp_path):
    positions = tuple(range(100, 1_000, 100))
    calls = [
        {
            "dataset": DATASET,
            "chrom": "chr22",
            "molecule_id": "all",
            "hp_family": "1",
            "phase_set": "100",
            "site_indices": ",".join(map(str, range(9))),
            "call_codes": "A" * 9,
        }
    ]
    source, molecule_calls = _build_source_test_run(
        tmp_path, "chr22", positions, calls
    )
    with molecule_calls.open("ab") as handle:
        handle.write(b"drift")
    output = tmp_path / "must_not_exist"

    completed = _run(source, output, "chr22", 500)
    assert completed.returncode == 2
    assert "identity mismatch" in completed.stderr
    assert not output.exists()


def test_production_rejects_noncanonical_cap_grid_before_output(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    output = tmp_path / "must_not_exist"
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--source-full-root",
            str(source),
            "--output-root",
            str(output),
            "--span-caps",
            "50000",
        ],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 2
    assert "caps must be exact 50000,100000,200000" in completed.stderr
    assert not output.exists()
