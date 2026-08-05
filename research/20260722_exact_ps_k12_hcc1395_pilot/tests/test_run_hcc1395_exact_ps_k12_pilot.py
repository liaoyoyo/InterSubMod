#!/usr/bin/env python3
"""Synthetic orchestration tests for the HCC1395 exact-PS k<=12 runner."""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import textwrap

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
RUNNER = TOPIC_ROOT / "scripts" / "run_hcc1395_exact_ps_k12_pilot.py"

SPEC = importlib.util.spec_from_file_location("hcc1395_exact_ps_runner", RUNNER)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_script(path: Path, source: str, *, executable: bool = False) -> Path:
    path.write_text(textwrap.dedent(source).lstrip(), encoding="utf-8")
    if executable:
        path.chmod(0o755)
    return path


def fake_verifier(path: Path) -> Path:
    return write_script(
        path,
        r"""
        import argparse
        import hashlib
        import json
        from pathlib import Path

        def canonical_sha256(value):
            payload = json.dumps(
                value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False
            ).encode("utf-8")
            return hashlib.sha256(payload).hexdigest()

        parser = argparse.ArgumentParser()
        parser.add_argument("--manifest", required=True)
        parser.add_argument("--output", required=True)
        parser.add_argument("--samtools", required=True)
        args = parser.parse_args()
        manifest = Path(args.manifest).resolve()
        manifest_payload = manifest.read_bytes()
        samtools = str(Path(args.samtools).resolve())
        receipt = {
            "schema_name": "intersubmod.big7_hcc1395_input_binding_receipt",
            "schema_version": "1.0.0",
            "task_type": "exploratory_pilot",
            "claim_status": "PARTIAL",
            "validation_evidence_eligible": False,
            "sample": "HCC1395",
            "all_pass": True,
            "failure": None,
            "manifest": {
                "path": str(manifest),
                "size_bytes": len(manifest_payload),
                "sha256": hashlib.sha256(manifest_payload).hexdigest(),
                "semantic_sha256": canonical_sha256(json.loads(manifest_payload)),
            },
            "source_authority": {"fixture": "authority"},
            "binding_contract": {"fixture": "binding"},
            "observed": {
                "samtools_quickcheck": {
                    "executable": samtools,
                    "returncode": 0,
                    "command": [samtools, "quickcheck", "fixture.bam"],
                }
            },
            "checks": [{"name": "fixture", "pass": True}],
        }
        receipt["receipt_integrity"] = {
            "policy": "semantic_json_sha256_without_receipt_integrity",
            "sha256": canonical_sha256(receipt),
        }
        Path(args.output).write_text(json.dumps(receipt) + "\n", encoding="utf-8")
        print(json.dumps(receipt))
        """,
    )


def fake_extractor(path: Path) -> Path:
    return write_script(
        path,
        r'''
        import argparse
        import csv
        import gzip
        import hashlib
        import json
        from pathlib import Path

        def write_gz(path, fields, rows):
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)

        def identity(path):
            payload = path.read_bytes()
            return {
                "path": str(path.resolve()),
                "size_bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }

        parser = argparse.ArgumentParser()
        parser.add_argument("--manifest", required=True)
        parser.add_argument("--sample", required=True)
        parser.add_argument("--chrom", required=True)
        parser.add_argument("--output-dir", required=True)
        parser.add_argument("--mapq-min")
        parser.add_argument("--baseq-min")
        parser.add_argument("--bridge-thresholds")
        parser.add_argument("--samtools-threads")
        parser.add_argument("--samtools", required=True)
        args = parser.parse_args()
        out = Path(args.output_dir)
        out.mkdir()
        prefix = f"{args.sample}.{args.chrom}"

        artifacts = []
        bed = out / f"{prefix}.targets.bed"
        bed.write_text(
            f"{args.chrom}\t99\t100\n{args.chrom}\t199\t200\n{args.chrom}\t299\t300\n",
            encoding="utf-8",
        )
        artifacts.append(bed)

        calls = out / f"{prefix}.molecule_sparse_calls.tsv.gz"
        write_gz(
            calls,
            ("dataset", "chrom", "molecule_id", "hp_family", "phase_set", "site_indices", "positions1", "call_codes"),
            [{
                "dataset": args.sample, "chrom": args.chrom, "molecule_id": "M1",
                "hp_family": "1", "phase_set": "100", "site_indices": "0,1,2",
                "positions1": "100,200,300", "call_codes": "AAA",
            }],
        )
        artifacts.append(calls)

        catalog = out / f"{prefix}.site_catalog.tsv.gz"
        write_gz(
            catalog,
            ("site_index", "chrom", "pos1", "ref", "alt"),
            [
                {"site_index": str(i), "chrom": args.chrom, "pos1": str(100 * (i + 1)), "ref": "C", "alt": "T"}
                for i in range(3)
            ],
        )
        artifacts.append(catalog)

        cuts = out / f"{prefix}.cut_support.tsv.gz"
        write_gz(cuts, ("cut_index", "left_pos1", "right_pos1"), [])
        artifacts.append(cuts)

        components = out / f"{prefix}.components.tsv.gz"
        write_gz(
            components,
            ("dataset", "chrom", "component_id"),
            [{"dataset": args.sample, "chrom": args.chrom, "component_id": "COMP1"}],
        )
        artifacts.append(components)

        membership = out / f"{prefix}.site_component_membership.tsv.gz"
        write_gz(
            membership,
            ("dataset", "chrom", "linkage_basis", "phase_set", "phase_set_status", "inference_role", "threshold", "site_index", "pos1", "component_id"),
            [
                {
                    "dataset": args.sample, "chrom": args.chrom, "linkage_basis": "PS_HP1",
                    "phase_set": "100", "phase_set_status": "KNOWN_PS_PRIMARY",
                    "inference_role": "PRIMARY_PS_AWARE", "threshold": "3",
                    "site_index": str(i), "pos1": str(100 * (i + 1)), "component_id": "COMP1",
                }
                for i in range(3)
            ],
        )
        artifacts.append(membership)

        receipt = {
            "schema_name": "intersubmod.lossless_read_linkage_chromosome_receipt",
            "schema_version": "1.3.0",
            "all_pass": True,
            "scope": {"dataset": args.sample, "chrom": args.chrom},
            "provenance": {
                "manifest": {
                    "path": str(Path(args.manifest).resolve()),
                    "sha256": hashlib.sha256(Path(args.manifest).read_bytes()).hexdigest(),
                }
            },
            "parameters": {"samtools_executable": str(Path(args.samtools).resolve())},
            "command": [str(Path(args.samtools).resolve()), "view", "fixture.bam"],
            "outputs": {item.name: identity(item) for item in artifacts},
        }
        receipt_path = out / "receipt.json"
        receipt_path.write_text(json.dumps(receipt, sort_keys=True) + "\n", encoding="utf-8")
        checksum = hashlib.sha256(receipt_path.read_bytes()).hexdigest()
        (out / "receipt.json.sha256").write_text(f"{checksum}  receipt.json\n", encoding="ascii")
        print(json.dumps(receipt))
        ''',
    )


def fake_partitioner(path: Path) -> Path:
    return write_script(
        path,
        r'''
        import argparse
        import csv
        import gzip
        import hashlib
        import json
        from pathlib import Path

        def write_gz(path, fields, rows):
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)

        def identity(path):
            payload = path.read_bytes()
            return {"path": str(path.resolve()), "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}

        parser = argparse.ArgumentParser()
        parser.add_argument("--dataset", required=True)
        parser.add_argument("--chrom", required=True)
        parser.add_argument("--site-catalog", required=True)
        parser.add_argument("--site-component-membership", required=True)
        parser.add_argument("--molecule-calls", required=True)
        parser.add_argument("--output-dir", required=True)
        parser.add_argument("--threshold", required=True)
        parser.add_argument("--max-block-size", required=True)
        args = parser.parse_args()
        out = Path(args.output_dir)
        out.mkdir()

        unit = {
            "dataset": args.dataset, "chrom": args.chrom, "unit_id": "U1",
            "component_id": "COMP1", "linkage_basis": "PS_HP1", "hp_family": "1",
            "phase_set": "100", "threshold": "3", "k": "3", "positions1": "100,200,300",
        }
        constraint = {
            "dataset": args.dataset, "chrom": args.chrom, "unit_id": "U1",
            "constraint_id": "C1", "hp_family": "1", "phase_set": "100",
            "positions1": "100,200,300", "call_codes": "AAA", "molecule_weight": "5",
            "pattern_count": "1",
        }
        block = {
            "dataset": args.dataset, "chrom": args.chrom, "unit_id": "U1",
            "component_id": "COMP1", "linkage_basis": "PS_HP1", "hp_family": "1",
            "phase_set": "100", "threshold": "3", "block_id": "B1", "block_index": "1",
            "start1": "100", "end1": "300", "k": "3", "positions1": "100,200,300",
            "retained_molecule_weight": "5", "retained_pattern_count": "1",
        }
        memberships = [
            {
                "dataset": args.dataset, "chrom": args.chrom, "unit_id": "U1",
                "component_id": "COMP1", "linkage_basis": "PS_HP1", "hp_family": "1",
                "phase_set": "100", "threshold": "3", "site_index": str(i),
                "pos1": str(100 * (i + 1)), "unit_local_index": str(i + 1),
                "block_id": "B1", "block_index": "1",
            }
            for i in range(3)
        ]
        disposition = {
            "dataset": args.dataset, "chrom": args.chrom, "unit_id": "U1",
            "constraint_id": "C1", "hp_family": "1", "phase_set": "100",
            "positions1": "100,200,300", "call_codes": "AAA", "molecule_weight": "5",
            "pattern_count": "1", "disposition": "retained", "span_sites": "3",
            "crossed_cut_count": "0", "retained_block_index": "1",
        }

        outputs = {
            "units.tsv.gz": (
                ("dataset", "chrom", "unit_id", "component_id", "linkage_basis", "hp_family", "phase_set", "threshold", "k", "positions1"),
                [unit],
            ),
            "constraints.tsv.gz": (
                ("dataset", "chrom", "unit_id", "constraint_id", "hp_family", "phase_set", "positions1", "call_codes", "molecule_weight", "pattern_count"),
                [constraint],
            ),
            "blocks.tsv.gz": (
                ("dataset", "chrom", "unit_id", "component_id", "linkage_basis", "hp_family", "phase_set", "threshold", "block_id", "block_index", "start1", "end1", "k", "positions1", "retained_molecule_weight", "retained_pattern_count"),
                [block],
            ),
            "membership.tsv.gz": (
                ("dataset", "chrom", "unit_id", "component_id", "linkage_basis", "hp_family", "phase_set", "threshold", "site_index", "pos1", "unit_local_index", "block_id", "block_index"),
                memberships,
            ),
            "dispositions.tsv.gz": (
                ("dataset", "chrom", "unit_id", "constraint_id", "hp_family", "phase_set", "positions1", "call_codes", "molecule_weight", "pattern_count", "disposition", "span_sites", "crossed_cut_count", "retained_block_index"),
                [disposition],
            ),
        }
        for name, (fields, rows) in outputs.items():
            write_gz(out / name, fields, rows)
        receipt = {
            "schema_name": "fake.python.partition",
            "schema_version": "1",
            "all_pass": True,
            "outputs": {name: identity(out / name) for name in outputs},
        }
        (out / "receipt.json").write_text(json.dumps(receipt) + "\n", encoding="utf-8")
        print(json.dumps(receipt))
        ''',
    )


def fake_compiler(path: Path) -> Path:
    return write_script(
        path,
        r"""
        #!/usr/bin/env python3
        from pathlib import Path
        import sys

        output = Path(sys.argv[sys.argv.index("-o") + 1])
        program = r'''#!/usr/bin/env python3
        import argparse
        import json
        from pathlib import Path

        parser = argparse.ArgumentParser()
        parser.add_argument("--units")
        parser.add_argument("--constraints")
        parser.add_argument("--output-dir")
        parser.add_argument("--max-block-size")
        args = parser.parse_args()
        out = Path(args.output_dir)
        out.mkdir()
        (out / "blocks.tsv").write_text("unit_id\\thp_family\\tphase_set\\nU1\\t1\\t100\\n", encoding="utf-8")
        (out / "membership.tsv").write_text("unit_id\\thp_family\\tphase_set\\nU1\\t1\\t100\\n", encoding="utf-8")
        (out / "dispositions.tsv").write_text("unit_id\\thp_family\\tphase_set\\nU1\\t1\\t100\\n", encoding="utf-8")
        (out / "summary.json").write_text(json.dumps({"schema_name": "fake.cpp", "all_pass": True}) + "\n", encoding="utf-8")
        print(json.dumps({"all_pass": True}))
        '''
        output.write_text(program, encoding="utf-8")
        output.chmod(0o755)
        """,
        executable=True,
    )


def fake_comparator(path: Path, *, all_pass: bool) -> Path:
    exit_code = 0 if all_pass else 1
    mismatch = 0 if all_pass else 1
    return write_script(
        path,
        f'''
        import argparse
        import json
        from pathlib import Path

        parser = argparse.ArgumentParser()
        parser.add_argument("--python-dir")
        parser.add_argument("--cpp-input-units")
        parser.add_argument("--cpp-input-constraints")
        parser.add_argument("--cpp-dir")
        parser.add_argument("--output")
        args = parser.parse_args()
        receipt = {{
            "schema_name": "fake.comparator",
            "schema_version": "1",
            "all_pass": {all_pass!r},
            "mismatch_count": {mismatch},
        }}
        Path(args.output).write_text(json.dumps(receipt) + "\\n", encoding="utf-8")
        raise SystemExit({exit_code})
        ''',
    )


def build_fake_tools(root: Path, *, comparator_pass: bool) -> dict[str, Path]:
    root.mkdir()
    return {
        "verifier": fake_verifier(root / "fake_verifier.py"),
        "extractor": fake_extractor(root / "fake_extractor.py"),
        "partitioner": fake_partitioner(root / "fake_partitioner.py"),
        "compiler": fake_compiler(root / "fake_cxx"),
        "comparator": fake_comparator(
            root / "fake_comparator.py", all_pass=comparator_pass
        ),
    }


def runner_command(
    tmp_path: Path,
    tools: dict[str, Path],
    output_root: Path,
    *,
    chromosomes: str,
    workers: int = 1,
) -> tuple[list[str], Path]:
    manifest = tmp_path / "big7_manifest.json"
    if not manifest.exists():
        manifest.write_text('{"fixture":"big7"}\n', encoding="utf-8")
    command = [
        sys.executable,
        str(RUNNER),
        "--manifest",
        str(manifest),
        "--output-root",
        str(output_root),
        "--chromosomes",
        chromosomes,
        "--workers",
        str(workers),
        "--input-verifier",
        str(tools["verifier"]),
        "--extractor",
        str(tools["extractor"]),
        "--partitioner",
        str(tools["partitioner"]),
        "--comparator",
        str(tools["comparator"]),
        "--cxx",
        str(tools["compiler"]),
        "--samtools",
        sys.executable,
    ]
    return command, manifest.resolve()


def run_fixture(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=TOPIC_ROOT.parents[1],
        text=True,
        capture_output=True,
        check=False,
        timeout=60,
    )


def read_summary(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_chromosome_list_validation_and_default() -> None:
    assert MODULE.parse_chromosomes("chr1,chr22") == ("chr1", "chr22")
    with pytest.raises(argparse.ArgumentTypeError, match="duplicates"):
        MODULE.parse_chromosomes("chr1,chr1")
    with pytest.raises(argparse.ArgumentTypeError, match="chr1-chr22"):
        MODULE.parse_chromosomes("chr1,chrX")
    parsed = MODULE.parse_args(
        ["--manifest", "/tmp/manifest", "--output-root", "/tmp/output"]
    )
    assert parsed.chromosomes == MODULE.AUTOSOMES
    assert parsed.workers == 1


def test_success_aggregate_and_big7_manifest_is_passed_to_children(
    tmp_path: Path,
) -> None:
    tools = build_fake_tools(tmp_path / "tools", comparator_pass=True)
    output_root = tmp_path / "run"
    command, manifest = runner_command(
        tmp_path,
        tools,
        output_root,
        chromosomes="chr1,chr2",
        workers=2,
    )

    completed = run_fixture(command)
    assert completed.returncode == 0, (
        f"command: {' '.join(command)}\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    )
    assert (output_root / "_SUCCESS").is_file()
    receipt = json.loads((output_root / "run_receipt.json").read_text())
    assert receipt["all_pass"] is True
    assert receipt["manifest"]["path"] == str(manifest)

    verifier_command = receipt["global_stages"][0]["command"]
    assert verifier_command[verifier_command.index("--manifest") + 1] == str(manifest)
    for chrom in ("chr1", "chr2"):
        child = json.loads(
            (output_root / "chromosomes" / chrom / "stage_receipt.json").read_text()
        )
        extract_command = child["stages"][0]["command"]
        assert extract_command[extract_command.index("--manifest") + 1] == str(manifest)

    aggregate = receipt["aggregate"]
    assert aggregate == {
        "S": 6,
        "blocks": 2,
        "chromosomes_failed": 0,
        "chromosomes_with_metrics": 2,
        "chromosomes_included_in_totals": 2,
        "failed_chromosomes_with_metrics_excluded": 0,
        "aggregation_policy": "PASS_chromosomes_with_complete_metrics_only",
        "chromosomes_passed": 2,
        "chromosomes_requested": 2,
        "constraint_dispositions": {"cut": 0, "retained": 2, "unavoidable": 0},
        "constraints": 2,
        "cross_hp_violations": 0,
        "cross_ps_violations": 0,
        "historical_comparison": {
            "artifacts_compared": 0,
            "physical_only_differences": 0,
            "semantic_mismatches": 0,
        },
        "k_bins": {"k=1": 0, "k=2-8": 2, "k=9-12": 0, "k>12": 0},
        "molecule_weights": {"cut": 0, "retained": 10, "total": 10, "unavoidable": 0},
        "pattern_count_total": 2,
        "python_cpp_mismatch_count": 0,
        "unique_sites": 6,
        "unit_memberships": 6,
        "units": 2,
    }
    rows = read_summary(output_root / "chromosome_summary.tsv")
    assert [row["chrom"] for row in rows] == ["chr1", "chr2"]
    assert {row["status"] for row in rows} == {"PASS"}
    assert {row["molecule_weight_total"] for row in rows} == {"5"}


def test_failed_child_receipt_preserved_and_no_success_marker(tmp_path: Path) -> None:
    tools = build_fake_tools(tmp_path / "tools", comparator_pass=False)
    output_root = tmp_path / "failed_run"
    command, _ = runner_command(
        tmp_path, tools, output_root, chromosomes="chr1", workers=1
    )

    completed = run_fixture(command)
    assert completed.returncode == 1
    assert not (output_root / "_SUCCESS").exists()
    comparison_path = output_root / "chromosomes" / "chr1" / "comparison.json"
    assert json.loads(comparison_path.read_text())["all_pass"] is False
    chromosome_receipt = json.loads(
        (output_root / "chromosomes" / "chr1" / "stage_receipt.json").read_text()
    )
    assert chromosome_receipt["all_pass"] is False
    assert chromosome_receipt["stages"][-1]["child_receipt"]["all_pass"] is False
    assert "all_pass=false" in chromosome_receipt["stages"][-1]["failure"]
    run_receipt = json.loads((output_root / "run_receipt.json").read_text())
    assert run_receipt["all_pass"] is False
    assert run_receipt["aggregate"]["chromosomes_failed"] == 1
    assert (output_root / "run_receipt.json.sha256").is_file()
    assert read_summary(output_root / "chromosome_summary.tsv")[0]["status"] == "FAIL"


def test_input_binding_schema_is_validated_beyond_child_all_pass(tmp_path: Path) -> None:
    tools = build_fake_tools(tmp_path / "tools", comparator_pass=True)
    output_root = tmp_path / "run"
    command, manifest = runner_command(
        tmp_path, tools, output_root, chromosomes="chr1", workers=1
    )
    completed = run_fixture(command)
    assert completed.returncode == 0

    receipt_path = output_root / "input_binding_receipt.pre.json"
    document = json.loads(receipt_path.read_text())
    document["schema_name"] = "untrusted.fake.schema"
    document.pop("receipt_integrity")
    document["receipt_integrity"] = {
        "policy": "semantic_json_sha256_without_receipt_integrity",
        "sha256": MODULE.canonical_json_sha256(document),
    }
    receipt_path.write_text(json.dumps(document) + "\n", encoding="utf-8")

    with pytest.raises(MODULE.RunnerError, match="schema/claim mismatch"):
        MODULE.validate_input_binding_receipt(
            receipt_path,
            expected_manifest=MODULE.file_identity(manifest),
            expected_samtools=Path(sys.executable).resolve(),
        )


def test_failed_chromosome_metrics_are_excluded_from_aggregate_totals() -> None:
    metrics = {field: 1 for field in MODULE.SUMMARY_FIELDS[2:-2]}
    aggregate = MODULE.aggregate_results(
        ("chr1", "chr2"),
        {
            "chr1": {"all_pass": True, "metrics": dict(metrics)},
            "chr2": {"all_pass": False, "metrics": dict(metrics)},
        },
    )
    assert aggregate["chromosomes_with_metrics"] == 2
    assert aggregate["chromosomes_included_in_totals"] == 1
    assert aggregate["failed_chromosomes_with_metrics_excluded"] == 1
    assert aggregate["aggregation_policy"] == "PASS_chromosomes_with_complete_metrics_only"
    assert aggregate["S"] == 1
    assert aggregate["molecule_weights"]["total"] == 1


def test_catalog_unique_sites_use_chrom_pos_and_reject_duplicate_locus(
    tmp_path: Path,
) -> None:
    tools = build_fake_tools(tmp_path / "tools", comparator_pass=True)
    output_root = tmp_path / "run"
    command, _ = runner_command(
        tmp_path, tools, output_root, chromosomes="chr1", workers=1
    )
    completed = run_fixture(command)
    assert completed.returncode == 0

    chrom_root = output_root / "chromosomes" / "chr1"
    catalog = chrom_root / "extraction" / "HCC1395.chr1.site_catalog.tsv.gz"
    with gzip.open(catalog, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        rows = list(reader)
    assert fields is not None
    rows[1]["pos1"] = rows[0]["pos1"]
    with gzip.open(catalog, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    with pytest.raises(MODULE.RunnerError, match=r"duplicate catalog \(chrom,pos1\)"):
        MODULE.chromosome_metrics(
            chrom="chr1",
            extraction_dir=chrom_root / "extraction",
            python_dir=chrom_root / "python_partition",
            comparison_path=chrom_root / "comparison.json",
            historical=None,
        )


def test_gzip_comparison_uses_decompressed_semantic_bytes(tmp_path: Path) -> None:
    first = tmp_path / "first.tsv.gz"
    second = tmp_path / "second.tsv.gz"
    payload = b"a\tb\n1\t2\n"
    with first.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            handle.write(payload)
    with second.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=123) as handle:
            handle.write(payload)

    comparison = MODULE.compare_semantic_artifacts(
        first, second, gzip_encoded=True
    )
    assert comparison["semantic_equal"] is True
    assert comparison["physical_equal"] is False
    assert comparison["physical_only_difference"] is True

    first_bed = tmp_path / "first.bed"
    second_bed = tmp_path / "second.bed"
    first_bed.write_bytes(b"chr1\t0\t1\n")
    second_bed.write_bytes(b"chr1\t0\t2\n")
    bed_comparison = MODULE.compare_semantic_artifacts(
        first_bed, second_bed, gzip_encoded=False
    )
    assert bed_comparison["semantic_equal"] is False
    assert bed_comparison["physical_only_difference"] is False


def test_historical_semantic_mismatch_fails_closed(tmp_path: Path) -> None:
    extraction_dir = tmp_path / "current"
    historical_root = tmp_path / "historical"
    logs_dir = tmp_path / "logs"
    extraction_dir.mkdir()
    historical_extract = historical_root / "chr1" / "extract"
    historical_extract.mkdir(parents=True)

    for suffix in MODULE.EXTRACTION_SUFFIXES:
        filename = f"HCC1395.chr1.{suffix}"
        current = extraction_dir / filename
        historical = historical_extract / filename
        payload = f"{suffix}\n".encode("utf-8")
        if suffix in MODULE.GZIP_EXTRACTION_SUFFIXES:
            with gzip.open(current, "wb") as handle:
                handle.write(payload)
            with gzip.open(historical, "wb") as handle:
                handle.write(payload)
        else:
            current.write_bytes(payload)
            historical.write_bytes(payload)

    mismatched = historical_extract / "HCC1395.chr1.targets.bed"
    mismatched.write_bytes(b"different\n")
    result = MODULE.historical_comparison_stage(
        chrom="chr1",
        extraction_dir=extraction_dir,
        historical_root=historical_root,
        logs_dir=logs_dir,
    )

    assert result["all_pass"] is False
    assert result["exit_code"] == 1
    assert result["artifact_count"] == len(MODULE.EXTRACTION_SUFFIXES)
    assert result["semantic_mismatch_count"] == 1
    assert result["scientific_equivalence_all"] is False
    assert "semantic_mismatches=1" in result["failure"]
