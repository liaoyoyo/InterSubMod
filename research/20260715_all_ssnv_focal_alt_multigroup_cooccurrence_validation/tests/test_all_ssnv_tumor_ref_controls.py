from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "analyze_all_ssnv_tumor_ref_controls.py"
)
SPEC = importlib.util.spec_from_file_location("analyze_all_ssnv_tumor_ref_controls", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


SITE_FIELDS = [
    "sample",
    "biological_id",
    "truth_label",
    "chrom",
    "pos",
    "ref",
    "alt",
    "region_dir",
    "screen_contract",
    "stable_null_multigroup",
    "component_id",
    "alt_readset_sha256",
    "cluster_sizes",
    "n_alt_raw",
    "n_alt_matrix",
]


def site_row(
    region_dir: Path,
    *,
    sample: str = "S1",
    pos: int = 100,
    stable: bool = True,
) -> dict[str, object]:
    return {
        "sample": sample,
        "biological_id": "BIO-1",
        "truth_label": "TP",
        "chrom": "chr1",
        "pos": pos,
        "ref": "A",
        "alt": "G",
        "region_dir": str(region_dir),
        "screen_contract": MODULE.SCREEN_CONTRACT,
        "stable_null_multigroup": stable,
        "component_id": "component-1",
        "alt_readset_sha256": "readset-1",
        "cluster_sizes": json.dumps({"1": 3, "2": 3}),
        "n_alt_raw": 6,
        "n_alt_matrix": 6,
    }


def assignment_for(row: dict[str, object]) -> dict[str, object]:
    return {
        "schema_name": MODULE.PRIMARY_ASSIGNMENT_SCHEMA,
        "schema_version": "1.0.0",
        "screen_contract": MODULE.SCREEN_CONTRACT,
        "sample": row["sample"],
        "chrom": row["chrom"],
        "pos": row["pos"],
        "region_dir": row["region_dir"],
        "read_ids": [f"alt-{index}" for index in range(6)],
        "labels": ["1", "1", "1", "2", "2", "2"],
        "posthoc": {"ref": row["ref"], "alt": row["alt"]},
    }


def write_gzip_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SITE_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_gzip_jsonl(path: Path, rows: list[dict[str, object]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row) + "\n")


def write_matrix(path: Path, ids: list[str], values: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["read_id", *[f"c{index}" for index in range(values.shape[1])]])
        for read_id, row in zip(ids, values):
            writer.writerow([read_id, *row.tolist()])


def make_region(tmp_path: Path) -> tuple[Path, list[str], set[str], set[str]]:
    region = tmp_path / "region"
    reads_path = region / "reads" / "reads.tsv"
    reads_path.parent.mkdir(parents=True)
    alt_ids = {f"alt-{index}" for index in range(6)}
    ref_ids = {f"ref-{index}" for index in range(6)}
    normal_ids = {"normal-ref"}
    ordered_ids = [
        value
        for index in range(6)
        for value in (f"alt-{index}", f"ref-{index}")
    ] + ["normal-ref"]
    with reads_path.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["read_id", "is_tumor", "alt_support"],
            delimiter="\t",
        )
        writer.writeheader()
        for read_id in ordered_ids:
            writer.writerow(
                {
                    "read_id": read_id,
                    "is_tumor": "0" if read_id in normal_ids else "1",
                    "alt_support": "ALT" if read_id in alt_ids else "REF",
                }
            )
    distance = np.full((len(ordered_ids), len(ordered_ids)), 0.2, dtype=float)
    np.fill_diagonal(distance, 0.0)
    methylation = np.tile(np.asarray([0.1, 0.9, 0.1, 0.9]), (len(ordered_ids), 1))
    write_matrix(region / "distance" / "BERNOULLI" / "matrix.csv", ordered_ids, distance)
    write_matrix(region / "methylation" / "methylation.csv", ordered_ids, methylation)
    return region, ordered_ids, alt_ids, ref_ids


def evaluable_result(read_ids: list[str], labels: list[str], stable: bool = True) -> dict[str, Any]:
    return {
        "status": "evaluable",
        "evaluable": True,
        "not_testable_reason": None,
        "n_raw": len(read_ids),
        "n_matrix": len(read_ids),
        "n_after_peel": len(read_ids),
        "n_peeled": 0,
        "kept_ids": list(read_ids),
        "labels": labels,
        "coarse_ng": 2 if stable else 1,
        "modal_fraction": 1.0,
        "unstable": False,
        "stable_null_multigroup": stable,
        "cluster_sizes": dict(MODULE.Counter(labels)),
    }


def minimal_control_row(**updates: Any) -> dict[str, Any]:
    row: dict[str, Any] = {
        "sample": "S1",
        "biological_id": "BIO-1",
        "truth_label": "TP",
        "chrom": "chr1",
        "pos": 100,
        "ref": "A",
        "alt": "G",
        "component_id": "component-1",
        "alt_readset_sha256": "readset-1",
        "n_tumor_alt": 6,
        "n_tumor_ref": 6,
        "ref_evaluable": True,
        "ref_stable_null_multigroup": False,
        "ref_not_testable_reason": None,
        "joint_evaluable": True,
        "joint_stable_null_multigroup": True,
        "joint_not_testable_reason": None,
        "joint_allele_testable": True,
        "joint_allele_axis_aligned": True,
        "joint_allele_not_testable_reason": None,
    }
    row.update(updates)
    return row


def test_gzip_inputs_load_only_complete_primary_stable_set(tmp_path: Path) -> None:
    stable = site_row(tmp_path / "stable")
    nonstable = site_row(tmp_path / "nonstable", pos=101, stable=False)
    site_path = tmp_path / "all_ssnv_site_results.tsv.gz"
    assignment_path = tmp_path / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
    write_gzip_tsv(site_path, [stable, nonstable])
    write_gzip_jsonl(assignment_path, [assignment_for(stable)])

    fields, tasks = MODULE.load_primary_inputs(site_path, assignment_path)

    assert fields == SITE_FIELDS
    assert len(tasks) == 1
    assert tasks[0]["site"]["pos"] == "100"
    assert tasks[0]["site"]["truth_label"] == "TP"
    assert tasks[0]["primary_assignment_n_core_reads"] == 6
    assert len(tasks[0]["primary_assignment_labels_sha256"]) == 64


def test_site_assignment_key_mismatch_is_hard_failure(tmp_path: Path) -> None:
    stable = site_row(tmp_path / "stable")
    wrong = site_row(tmp_path / "stable", pos=999)
    site_path = tmp_path / "all_ssnv_site_results.tsv.gz"
    assignment_path = tmp_path / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
    write_gzip_tsv(site_path, [stable])
    write_gzip_jsonl(assignment_path, [assignment_for(wrong)])

    with pytest.raises(RuntimeError, match="site/assignment key mismatch"):
        MODULE.load_primary_inputs(site_path, assignment_path)


def test_missing_schema_or_screen_contract_is_hard_failure(tmp_path: Path) -> None:
    stable = site_row(tmp_path / "stable")
    site_path = tmp_path / "all_ssnv_site_results.tsv.gz"
    assignment_path = tmp_path / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"

    assignment = assignment_for(stable)
    assignment.pop("schema_version")
    write_gzip_tsv(site_path, [stable])
    write_gzip_jsonl(assignment_path, [assignment])
    with pytest.raises(RuntimeError, match="Unexpected stable assignment schema"):
        MODULE.load_primary_inputs(site_path, assignment_path)

    assignment = assignment_for(stable)
    stable_without_contract = dict(stable)
    stable_without_contract["screen_contract"] = ""
    write_gzip_tsv(site_path, [stable_without_contract])
    write_gzip_jsonl(assignment_path, [assignment])
    with pytest.raises(RuntimeError, match="screen-contract drift"):
        MODULE.load_primary_inputs(site_path, assignment_path)


def test_duplicate_primary_site_and_assignment_keys_are_rejected(tmp_path: Path) -> None:
    stable = site_row(tmp_path / "stable")
    site_path = tmp_path / "all_ssnv_site_results.tsv.gz"
    assignment_path = tmp_path / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
    write_gzip_tsv(site_path, [stable, stable])
    write_gzip_jsonl(assignment_path, [assignment_for(stable)])
    with pytest.raises(RuntimeError, match="Duplicate primary stable site key"):
        MODULE.load_primary_inputs(site_path, assignment_path)

    write_gzip_tsv(site_path, [stable])
    write_gzip_jsonl(assignment_path, [assignment_for(stable), assignment_for(stable)])
    with pytest.raises(RuntimeError, match="Duplicate stable assignment site key"):
        MODULE.load_primary_inputs(site_path, assignment_path)


def test_tumor_ref_selection_excludes_alt_and_normal_reads() -> None:
    reads = {
        "tumor-ref": {"is_tumor": "1", "alt_support": "REF"},
        "tumor-alt": {"is_tumor": "true", "alt_support": "ALT"},
        "normal-ref": {"is_tumor": "0", "alt_support": "REF"},
        "normal-alt": {"is_tumor": "false", "alt_support": "ALT"},
    }

    assert MODULE.select_tumor_read_ids(reads, "REF") == ["tumor-ref"]
    assert MODULE.select_tumor_read_ids(reads, "ALT") == ["tumor-alt"]


def test_analyze_subset_uses_primary_column_null_contract(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ids = [f"read-{index}" for index in range(6)]
    distance = np.full((6, 6), 0.2, dtype=float)
    np.fill_diagonal(distance, 0.0)
    methylation = np.tile(np.asarray([0.1, 0.9, 0.1]), (6, 1))
    captured: dict[str, Any] = {}

    def fake_analyze_phylo(*args: Any, **kwargs: Any) -> dict[str, Any]:
        captured.update(kwargs)
        return {
            "coarse_labels": ["1"] * 6,
            "coarse_ng": 1,
            "modal_fraction": 1.0,
            "unstable": False,
        }

    monkeypatch.setattr(MODULE.F, "analyze_phylo", fake_analyze_phylo)
    result = MODULE.analyze_subset(ids, ids, distance, ids, methylation, 17)

    assert result["evaluable"] is True
    assert captured == {
        "base_seed": 17,
        "seeds": 10,
        "rnull": MODULE.F.RNULL,
        "null_mode": "column",
        "empirical_alpha": None,
    }


def test_joint_clustering_association_is_posthoc_and_primary_fields_are_unchanged(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    region, ordered_ids, alt_ids, ref_ids = make_region(tmp_path)
    source = site_row(region)
    calls: list[list[str]] = []

    def fake_analyze_subset(read_ids: list[str], *args: Any, **kwargs: Any) -> dict[str, Any]:
        del args, kwargs
        calls.append(list(read_ids))
        if set(read_ids) == ref_ids:
            labels = ["ref-1"] * 3 + ["ref-2"] * 3
        else:
            labels = ["alt-group" if read_id in alt_ids else "ref-group" for read_id in read_ids]
        return evaluable_result(list(read_ids), labels)

    monkeypatch.setattr(MODULE, "analyze_subset", fake_analyze_subset)
    task = {
        "site": {key: str(value) for key, value in source.items()},
        "primary_assignment_n_core_reads": 6,
        "primary_assignment_labels_sha256": "f" * 64,
    }
    result = MODULE.process_site_task(task)

    assert set(calls[0]) == ref_ids
    expected_joint = [read_id for read_id in ordered_ids if read_id in alt_ids | ref_ids]
    assert calls[1] == expected_joint
    assert result["n_tumor_alt"] == 6
    assert result["n_tumor_ref"] == 6
    assert result["ref_stable_null_multigroup"] is True
    assert result["joint_stable_null_multigroup"] is True
    assert result["joint_allele_testable"] is True
    assert result["joint_allele_v"] == pytest.approx(1.0)
    assert result["joint_allele_p_perm"] < 0.05
    assert result["cluster_sizes"] == source["cluster_sizes"]
    assert task["site"]["cluster_sizes"] == source["cluster_sizes"]
    assert result["background_control_interpretation"] == "REF_REPLICATION_WEAKENS_ALT_SPECIFICITY"


def test_joint_association_reports_explicit_not_testable_reason() -> None:
    result = MODULE.compute_joint_allele_association(
        [f"alt-{index}" for index in range(6)],
        ["group-1"] * 3 + ["group-2"] * 3,
        {f"alt-{index}": "ALT" for index in range(6)},
        11,
        permutations=19,
    )

    assert result["testable"] is False
    assert result["v"] is None
    assert result["not_testable_reason"] == (
        "joint_core_lacks_both_focal_alleles_or_multiple_clusters"
    )


def test_bounded_execution_never_submits_all_chunks_at_once() -> None:
    class FakeFuture:
        def __init__(self, executor: Any, value: list[dict[str, Any]]) -> None:
            self.executor = executor
            self.value = value
            self.resolved = False

        def result(self) -> list[dict[str, Any]]:
            if not self.resolved:
                self.executor.active -= 1
                self.resolved = True
            return self.value

    class FakeExecutor:
        def __init__(self) -> None:
            self.active = 0
            self.max_active = 0
            self.submissions = 0

        def submit(self, function: Any, chunk: list[dict[str, Any]]) -> FakeFuture:
            del function
            self.active += 1
            self.submissions += 1
            self.max_active = max(self.max_active, self.active)
            return FakeFuture(self, chunk)

    executor = FakeExecutor()
    chunks = [[{"index": index}] for index in range(7)]
    observed = list(MODULE.bounded_chunk_results(executor, chunks, max_pending=2))

    assert observed == chunks
    assert executor.submissions == 7
    assert executor.max_active == 2


def test_component_and_readset_denominators_are_deduplicated() -> None:
    summary = MODULE.StratumSummary()
    summary.add(minimal_control_row(sample="S1", pos=100))
    summary.add(minimal_control_row(sample="S1-platform2", pos=100))
    summary.add(
        minimal_control_row(
            sample="S2",
            biological_id="BIO-2",
            pos=200,
            component_id="component-2",
            alt_readset_sha256="readset-2",
        )
    )
    payload = summary.payload()

    assert payload["site_weighted"]["n_sites"] == 3
    assert payload["component_deduplicated"]["n_sites"] == 2
    assert payload["readset_deduplicated"]["n_sites"] == 2


def test_synthetic_main_writes_gzip_summary_and_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    region, _, alt_ids, ref_ids = make_region(tmp_path)
    source = site_row(region)
    site_path = tmp_path / "all_ssnv_site_results.tsv.gz"
    assignment_path = tmp_path / "all_ssnv_stable_multigroup_read_assignments.jsonl.gz"
    primary_audit_pre = tmp_path / "primary-artifact-audit-pre.json"
    output_dir = tmp_path / "control-output"
    write_gzip_tsv(site_path, [source])
    write_gzip_jsonl(assignment_path, [assignment_for(source)])
    primary_audit_pre.write_text(
        json.dumps({"schema_name": "synthetic-primary-audit", "pass": True}) + "\n",
        encoding="utf-8",
    )

    def fake_analyze_subset(read_ids: list[str], *args: Any, **kwargs: Any) -> dict[str, Any]:
        del args, kwargs
        if set(read_ids) == ref_ids:
            labels = ["ref-1"] * 3 + ["ref-2"] * 3
        else:
            labels = ["alt-group" if read_id in alt_ids else "ref-group" for read_id in read_ids]
        return evaluable_result(list(read_ids), labels)

    class InlineExecutor:
        def __init__(self, max_workers: int) -> None:
            assert max_workers == 1

        def __enter__(self) -> "InlineExecutor":
            return self

        def __exit__(self, *args: Any) -> None:
            del args

        def submit(self, function: Any, chunk: list[dict[str, Any]]) -> Any:
            future = MODULE.Future()
            try:
                future.set_result(function(chunk))
            except Exception as error:
                future.set_exception(error)
            return future

    monkeypatch.setattr(MODULE, "analyze_subset", fake_analyze_subset)
    monkeypatch.setattr(MODULE, "ProcessPoolExecutor", InlineExecutor)
    monkeypatch.setattr(
        MODULE,
        "parse_args",
        lambda: SimpleNamespace(
            site_results=site_path,
            stable_assignments=assignment_path,
            primary_artifact_audit_pre=primary_audit_pre,
            output_dir=output_dir,
            workers=1,
            chunk_size=1,
            max_pending_chunks=1,
            progress_every=1,
        ),
    )

    MODULE.main()

    table_path = output_dir / "all_ssnv_tumor_ref_control_site_results.tsv.gz"
    summary_path = output_dir / "all_ssnv_tumor_ref_control_summary.json"
    manifest_path = output_dir / "run_manifest.json"
    with gzip.open(table_path, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert len(rows) == 1
    assert rows[0]["truth_label"] == "TP"
    assert rows[0]["ref_evaluable"] == "true"
    assert rows[0]["joint_allele_testable"] == "true"
    assert summary["pooled"]["site_weighted"]["n_sites"] == 1
    assert summary["pooled"]["component_deduplicated"]["n_sites"] == 1
    assert summary["guardrail"]["ref_nonreplication"].endswith("does not prove a subclone.")
    assert manifest["counts"] == {"primary_stable_sites": 1, "processed_sites": 1}
    audit_record = manifest["inputs"]["primary_artifact_audit_pre"]
    assert Path(audit_record["path"]).resolve() == primary_audit_pre.resolve()
    assert audit_record["sha256"] == MODULE.sha256(primary_audit_pre)
    assert manifest["created_at_utc"] == manifest["finished_at_utc"]
    assert datetime.fromisoformat(manifest["started_at_utc"]) <= datetime.fromisoformat(
        manifest["finished_at_utc"]
    )
    assert manifest["pass"] is True


def test_cli_requires_primary_artifact_audit_pre(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(sys, "argv", [str(SCRIPT)])
    with pytest.raises(SystemExit) as error:
        MODULE.parse_args()
    assert error.value.code == 2


def test_existing_output_directory_is_never_overwritten(tmp_path: Path) -> None:
    output = tmp_path / "existing"
    output.mkdir()

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        MODULE.create_output_dir(output)
