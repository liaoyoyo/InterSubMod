from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys


TOPIC_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = TOPIC_ROOT / "scripts" / "summarize_exact_ps_cpp_topology_cohort.py"


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )


def write_sidecar(path: Path) -> None:
    path.with_name(f"{path.name}.sha256").write_text(
        f"{sha256_path(path)}  {path.name}\n",
        encoding="ascii",
    )


def build_completed_sample(root: Path, sample: str = "SYNTH") -> Path:
    sample_root = root / "samples" / sample
    sample_root.mkdir(parents=True)
    mlhp = sample_root / f"{sample}.exact_ps_mlhp.json"
    mlhp.write_text('{"groups":[{},{},{},{},{}],"sample":"SYNTH"}\n', encoding="utf-8")
    mlhp_receipt = sample_root / f"{sample}.exact_ps_mlhp.json.receipt.json"
    write_json(
        mlhp_receipt,
        {
            "schema_name": "intersubmod.exact_ps_layered_tree_input.receipt",
            "schema_version": "1.0.0",
            "all_pass": True,
            "counts": {"groups_tree_input": 5},
            "output": identity(mlhp),
        },
    )

    common = {
        "schema_name": "intersubmod.exact_ps_cpp_topology_af.unit",
        "schema_version": "1.0.0",
        "sample": sample,
    }
    rows = [
        {
            **common,
            "active_bit_count": 0,
            "unit_status": "no_active_alt",
            "family_status": "FAMILY_COMPLETE",
            "objective_status": "OBJECTIVE_CERTIFIED",
            "read_af_status": "not_applicable_no_active_alt",
            "best_tree_unique": None,
            "search_nodes": 1,
            "solver_elapsed_microseconds": 5,
        },
        {
            **common,
            "active_bit_count": 2,
            "unit_status": "ranked",
            "family_status": "FAMILY_COMPLETE",
            "objective_status": "OBJECTIVE_CERTIFIED",
            "read_af_status": "ranked_complete",
            "best_tree_unique": True,
            "best_tree_tie_count": "1",
            "representative_best_morphology": {
                "has_sister": False,
                "has_direct": True,
            },
            "search_nodes": 3,
            "solver_elapsed_microseconds": 10,
        },
        {
            **common,
            "active_bit_count": 3,
            "unit_status": "ranked",
            "family_status": "FAMILY_COMPLETE",
            "objective_status": "OBJECTIVE_CERTIFIED",
            "read_af_status": "ranked_complete",
            "best_tree_unique": False,
            "best_tree_tie_count": "7",
            "representative_best_morphology": {
                "has_sister": True,
                "has_direct": True,
            },
            "search_nodes": 10,
            "solver_elapsed_microseconds": 20,
        },
        {
            **common,
            "active_bit_count": 2,
            "unit_status": "zero_denominator",
            "family_status": "FAMILY_COMPLETE",
            "objective_status": "OBJECTIVE_CERTIFIED",
            "read_af_status": "zero_denominator",
            "best_tree_unique": None,
            "search_nodes": 3,
            "solver_elapsed_microseconds": 8,
        },
        {
            **common,
            "active_bit_count": 4,
            "unit_status": "family_incomplete",
            "family_status": "ABSTAIN_RESOURCE_LIMIT",
            "objective_status": "ABSTAIN_RESOURCE_LIMIT",
            "read_af_status": "not_evaluable_family_incomplete",
            "best_tree_unique": None,
            "search_nodes": 1000,
            "solver_elapsed_microseconds": 100,
        },
    ]
    topology = sample_root / f"{sample}.topology.jsonl"
    topology.write_text(
        "".join(
            json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n"
            for row in rows
        ),
        encoding="utf-8",
    )
    topology_receipt = sample_root / f"{sample}.topology.receipt.json"
    write_json(
        topology_receipt,
        {
            "schema_name": "intersubmod.exact_ps_cpp_topology_af.receipt",
            "schema_version": "1.0.0",
            "all_pass": True,
            "all_mutation_bearing_families_complete": False,
            "all_units_objective_certified": False,
            "input": identity(mlhp),
            "output": identity(topology),
            "counts": {
                "total_units": 5,
                "malformed_units": 0,
                "mutation_bearing_units": 4,
                "objective_certified_units": 4,
                "family_complete_mutation_bearing_units": 3,
                "ranking_complete_count": 2,
            },
            "status_census": {
                "unit_status": {
                    "no_active_alt": 1,
                    "ranked": 2,
                    "zero_denominator": 1,
                    "family_incomplete": 1,
                },
                "family_status": {
                    "FAMILY_COMPLETE": 4,
                    "ABSTAIN_RESOURCE_LIMIT": 1,
                },
                "objective_status": {
                    "OBJECTIVE_CERTIFIED": 4,
                    "ABSTAIN_RESOURCE_LIMIT": 1,
                },
                "read_af_status": {
                    "not_applicable_no_active_alt": 1,
                    "ranked_complete": 2,
                    "zero_denominator": 1,
                    "not_evaluable_family_incomplete": 1,
                },
            },
            "parameters": {
                "max_family_size": 100000,
                "max_search_nodes": 1000,
            },
            "runtime": {"elapsed_seconds": 0.5},
        },
    )
    pipeline = sample_root / "pipeline_receipt.json"
    write_json(
        pipeline,
        {
            "schema_name": (
                "intersubmod.exact_ps_cpp_topology_cohort.sample_pipeline_receipt"
            ),
            "schema_version": "1.0.0",
            "all_pass": True,
            "technical_all_pass": True,
            "sample": sample,
            "all_mutation_bearing_families_complete": False,
            "all_units_objective_certified": False,
            "counts": {"mlhp_groups": 5, "topology_jsonl_rows": 5},
            "stages": [
                {
                    "stage": "exact_ps_partition_to_mlhp",
                    "wall_seconds": 1.25,
                },
                {"stage": "cpp_topology", "wall_seconds": 0.75},
            ],
            "bound_artifacts": [
                identity(path)
                for path in (mlhp, mlhp_receipt, topology, topology_receipt)
            ],
        },
    )
    write_sidecar(pipeline)
    return sample_root


def build_cohort_receipt(root: Path, sample: str = "SYNTH") -> None:
    pipeline = root / "samples" / sample / "pipeline_receipt.json"
    receipt = root / "cohort_receipt.json"
    write_json(
        receipt,
        {
            "schema_name": "intersubmod.exact_ps_cpp_topology_cohort.cohort_receipt",
            "schema_version": "1.0.0",
            "all_pass": True,
            "technical_all_pass": True,
            "all_mutation_bearing_families_complete": False,
            "all_units_objective_certified": False,
            "scope": {
                "samples": [sample],
                "chromosomes": ["chr1"],
                "canonical_all7": False,
                "autosomes_complete": False,
            },
            "sample_receipts": {
                sample: {
                    "partition": None,
                    "pipeline": identity(pipeline),
                }
            },
        },
    )
    write_sidecar(receipt)


def run_summary(
    root: Path,
    output_root: Path,
    *,
    partial: bool,
) -> subprocess.CompletedProcess[str]:
    command = [
        sys.executable,
        str(SCRIPT),
        "--cohort-root",
        str(root),
        "--output-json",
        str(output_root / "summary.json"),
        "--output-tsv",
        str(output_root / "summary.tsv"),
    ]
    if partial:
        command.extend(["--samples", "SYNTH", "--allow-partial"])
    return subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )


def test_completed_subset_summary_conserves_counts_and_hides_tied_morphology(
    tmp_path: Path,
) -> None:
    root = tmp_path / "cohort"
    build_completed_sample(root)
    completed = run_summary(root, tmp_path / "out", partial=True)
    assert completed.returncode == 0, completed.stderr
    summary = json.loads((tmp_path / "out" / "summary.json").read_text())
    sample = summary["samples"][0]
    assert sample["groups_total"] == 5
    assert sample["mutation_bearing_units"] == 4
    assert sample["mutation_family_complete_units"] == 3
    assert sample["mutation_family_abstain_units"] == 1
    assert sample["ranked_units"] == 2
    assert sample["best_tree_unique_units"] == 1
    assert sample["best_tree_tied_units"] == 1
    assert sample["unique_best_morphology"]["direct_only"] == 1
    assert sample["unique_best_morphology"]["sister_plus_direct"] == 0
    assert sample["morphology_unresolved_tied"] == 1
    assert sample["active_k_distribution"]["4"] == 1
    assert sample["search_nodes"]["max"] == 1000
    assert sample["solver_elapsed_microseconds"]["max"] == 100
    assert summary["scope"]["partial"] is True
    assert summary["validation_evidence_eligible"] is False
    assert all(summary["checks"].values())


def test_completed_cohort_receipt_mode_passes(tmp_path: Path) -> None:
    root = tmp_path / "cohort"
    build_completed_sample(root)
    build_cohort_receipt(root)
    completed = run_summary(root, tmp_path / "out", partial=False)
    assert completed.returncode == 0, completed.stderr
    summary = json.loads((tmp_path / "out" / "summary.json").read_text())
    assert summary["scope"]["partial"] is False
    assert summary["cohort_receipt_identity"]["sha256"] == sha256_path(
        root / "cohort_receipt.json"
    )


def test_tampered_topology_fails_without_outputs(tmp_path: Path) -> None:
    root = tmp_path / "cohort"
    sample_root = build_completed_sample(root)
    topology = sample_root / "SYNTH.topology.jsonl"
    with topology.open("a", encoding="utf-8") as handle:
        handle.write("{}\n")
    completed = run_summary(root, tmp_path / "out", partial=True)
    assert completed.returncode == 2
    assert (
        "SHA256" in completed.stderr
        or "size mismatch" in completed.stderr
        or "schema_name mismatch" in completed.stderr
    )
    assert not (tmp_path / "out" / "summary.json").exists()
    assert not (tmp_path / "out" / "summary.tsv").exists()


def test_unfinished_sample_is_rejected_before_downstream_read(tmp_path: Path) -> None:
    root = tmp_path / "cohort"
    sample_root = root / "samples" / "SYNTH"
    sample_root.mkdir(parents=True)
    # Deliberately malformed downstream artifact: the missing pipeline gate must win.
    (sample_root / "SYNTH.topology.jsonl").write_text("not-json\n", encoding="utf-8")
    completed = run_summary(root, tmp_path / "out", partial=True)
    assert completed.returncode == 2
    assert "pipeline receipt absent" in completed.stderr
    assert "invalid topology" not in completed.stderr
