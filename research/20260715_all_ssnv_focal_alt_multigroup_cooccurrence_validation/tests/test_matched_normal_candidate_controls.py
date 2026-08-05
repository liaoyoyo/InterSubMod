from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
import pysam
import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str):
    path = TOPIC_ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


RUNNER = load_script("run_matched_normal_candidate_controls")
ANALYZER = load_script("analyze_matched_normal_candidate_controls")


def test_current_release_matched_normal_paths_bind_v7_and_v3() -> None:
    assert RUNNER.CANONICAL_CANDIDATE_TABLE.parent.name == (
        "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_"
        "source_locked_command_parity"
    )
    assert (
        RUNNER.CANONICAL_OUTPUT_ROOT.name
        == "matched_normal_candidate_controls_v3_m2v5_source_authority_v5"
    )
    assert ANALYZER.CANONICAL_PAIRED_OUTPUT_ROOT == RUNNER.CANONICAL_OUTPUT_ROOT
    assert ANALYZER.CANONICAL_OUTPUT_DIR.name == (
        "matched_normal_candidate_control_analysis_v3_m2v5_source_authority_v5"
    )
    expected_prefix = [
        sys.executable,
        "-I",
        "-X",
        f"pycache_prefix={RUNNER.CANONICAL_PYTHON_CACHE_ROOT}",
    ]
    assert RUNNER.canonical_python_prefix() == expected_prefix
    assert ANALYZER.canonical_python_prefix() == expected_prefix
    assert RUNNER.canonical_task_b_command()[:4] == expected_prefix
    assert ANALYZER.canonical_task_b_command()[:4] == expected_prefix


@pytest.fixture(autouse=True)
def test_only_source_authority(monkeypatch: pytest.MonkeyPatch) -> None:
    authority = {"authority_id": "TEST_ONLY_UNSIGNED_AUTHORITY", "pass": True}
    monkeypatch.setattr(
        RUNNER.SOURCE_AUTHORITY,
        "validate_release_source_authority",
        lambda *_args, **_kwargs: authority,
    )


def write_vcf(path: Path, records: list[tuple[str, int, str, str]]) -> None:
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")
    for chrom in sorted({record[0] for record in records}):
        header.contigs.add(chrom)
    with pysam.VariantFile(str(path), "w", header=header) as destination:
        for chrom, pos, ref, alt in records:
            record = destination.new_record(
                contig=chrom,
                start=pos - 1,
                stop=pos - 1 + len(ref),
                alleles=(ref, alt),
            )
            record.filter.add("PASS")
            destination.write(record)


def write_candidate_table(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    if path.name.endswith(".gz"):
        context = gzip.open(path, "wt", encoding="utf-8", newline="")
    else:
        context = path.open("w", encoding="utf-8", newline="")
    with context as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def zero_selection_runner_args(tmp_path: Path) -> tuple[list[str], dict[str, Path]]:
    selection_column = "multi_marker_molecular_haplotype_base_candidate"
    candidate_table = tmp_path / "candidates.tsv"
    write_candidate_table(
        candidate_table,
        ["sample", "chrom", "pos", "ref", "alt", selection_column],
        [
            {
                "sample": "HCC1395",
                "chrom": "chr1",
                "pos": 10,
                "ref": "A",
                "alt": "G",
                selection_column: "false",
            }
        ],
    )
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_name": "intersubmod.all_ssnv_focal_alt_input_manifest",
                "pass": True,
                "workspace": str(tmp_path / "workspace"),
                "samples": [{"sample": sample} for sample in RUNNER.DATASETS],
            }
        ),
        encoding="utf-8",
    )
    normal_audit = tmp_path / "normal_audit.json"
    normal_audit.write_text(
        json.dumps(
            {
                "schema_name": "intersubmod.matched_normal_methyl_tag_audit",
                "pass": True,
                "totals": {"n_samples": 7, "n_normal_control_eligible": 7},
                "samples": [
                    {
                        "sample": sample,
                        "bam": str(tmp_path / f"missing-{sample}.bam"),
                        "index": str(tmp_path / f"missing-{sample}.bam.bai"),
                        "normal_control_eligible": True,
                    }
                    for sample in RUNNER.DATASETS
                ],
            }
        ),
        encoding="utf-8",
    )
    binary = tmp_path / "inter_sub_mod"
    binary.write_bytes(b"synthetic binary\n")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nA\n", encoding="ascii")
    Path(str(reference) + ".fai").write_text("chr1\t1\t6\t1\t2\n", encoding="ascii")
    output_root = tmp_path / "matched-normal-output"
    paths = {
        "candidate_table": candidate_table,
        "all_ssnv_manifest": manifest,
        "normal_audit": normal_audit,
        "binary": binary,
        "reference": reference,
        "runner_script": Path(RUNNER.__file__).resolve(),
        "output_root": output_root,
    }
    args = [
        "runner",
        "--candidate-table",
        str(candidate_table),
        "--selection-column",
        selection_column,
        "--selection-value",
        "true",
        "--manifest",
        str(manifest),
        "--normal-audit",
        str(normal_audit),
        "--binary",
        str(binary),
        "--reference",
        str(reference),
        "--output-root",
        str(output_root),
    ]
    return args, paths


def test_exact_vcf_selection_preserves_only_requested_records(tmp_path: Path) -> None:
    source = tmp_path / "source.vcf"
    write_vcf(
        source,
        [("chr1", 10, "A", "G"), ("chr1", 20, "C", "T"), ("chr1", 30, "G", "A")],
    )
    table = tmp_path / "candidates.tsv.gz"
    write_candidate_table(
        table,
        ["sample", "chrom", "pos", "ref", "alt"],
        [
            {"sample": "S1", "chrom": "chr1", "pos": 10, "ref": "a", "alt": "g"},
            {"sample": "S1", "chrom": "chr1", "pos": 30, "ref": "G", "alt": "A"},
        ],
    )
    candidates = RUNNER.read_candidate_table(table, allowed_samples={"S1"})
    output = tmp_path / "selected.vcf.gz"
    receipt = RUNNER.write_candidate_vcf(source, output, candidates)

    with pysam.VariantFile(str(output)) as selected:
        observed = [
            (record.chrom, int(record.pos), record.ref, record.alts[0]) for record in selected
        ]
    assert observed == [("chr1", 10, "A", "G"), ("chr1", 30, "G", "A")]
    assert receipt["site_count"] == 2
    assert len(receipt["sha256"]) == 64
    assert Path(receipt["index"]).exists()


def test_vcf_selection_rejects_duplicate_missing_and_allele_mismatch(tmp_path: Path) -> None:
    source = tmp_path / "source.vcf"
    write_vcf(source, [("chr1", 10, "A", "G")])
    duplicate = tmp_path / "duplicate.tsv"
    row = {"sample": "S1", "chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"}
    write_candidate_table(
        duplicate,
        ["sample", "chrom", "pos", "ref", "alt"],
        [row, row],
    )
    with pytest.raises(ValueError, match="Duplicate candidate row"):
        RUNNER.read_candidate_table(duplicate, allowed_samples={"S1"})
    with pytest.raises(RuntimeError, match="REF/ALT mismatch"):
        RUNNER.inspect_candidate_matches(
            source, [RUNNER.Candidate("S1", "chr1", 10, "A", "T")]
        )
    with pytest.raises(RuntimeError, match="does not exist"):
        RUNNER.inspect_candidate_matches(
            source, [RUNNER.Candidate("S1", "chr1", 11, "A", "G")]
        )


def test_boolean_selection_supports_text_and_numeric_values(tmp_path: Path) -> None:
    table = tmp_path / "selection.tsv"
    write_candidate_table(
        table,
        ["sample", "chrom", "pos", "ref", "alt", "selected"],
        [
            {"sample": "S1", "chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "selected": "TRUE"},
            {"sample": "S1", "chrom": "chr1", "pos": 20, "ref": "A", "alt": "C", "selected": "false"},
            {"sample": "S1", "chrom": "chr1", "pos": 30, "ref": "C", "alt": "T", "selected": "1"},
            {"sample": "S1", "chrom": "chr1", "pos": 40, "ref": "G", "alt": "A", "selected": "0"},
        ],
    )
    selected_true = RUNNER.read_candidate_table(
        table, "selected", "true", allowed_samples={"S1"}
    )
    selected_false = RUNNER.read_candidate_table(
        table, "selected", "FALSE", allowed_samples={"S1"}
    )
    implicit_true = RUNNER.read_candidate_table(
        table, "selected", allowed_samples={"S1"}
    )
    assert [candidate.pos for candidate in selected_true] == [10, 30]
    assert [candidate.pos for candidate in selected_false] == [20, 40]
    assert implicit_true == selected_true


def test_zero_selection_with_explicit_column_returns_empty(tmp_path: Path) -> None:
    table = tmp_path / "zero-selection.tsv"
    write_candidate_table(
        table,
        ["sample", "chrom", "pos", "ref", "alt", "selected"],
        [
            {
                "sample": "S1",
                "chrom": "chr1",
                "pos": 10,
                "ref": "A",
                "alt": "G",
                "selected": "false",
            }
        ],
    )
    assert RUNNER.read_candidate_table(
        table, "selected", "true", allowed_samples={"S1"}
    ) == []


def test_missing_selection_column_remains_a_hard_failure(tmp_path: Path) -> None:
    table = tmp_path / "missing-selection-column.tsv"
    write_candidate_table(
        table,
        ["sample", "chrom", "pos", "ref", "alt"],
        [{"sample": "S1", "chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"}],
    )
    with pytest.raises(ValueError, match="missing selection column"):
        RUNNER.read_candidate_table(
            table, "selected", "true", allowed_samples={"S1"}
        )


def test_zero_selection_main_writes_only_sha_locked_not_applicable_receipt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    args, paths = zero_selection_runner_args(tmp_path)
    monkeypatch.setattr(sys, "argv", args)

    def reject_subprocess(*_args: object, **_kwargs: object) -> None:
        pytest.fail("Zero-selected-candidate path must not execute C++")

    monkeypatch.setattr(RUNNER.subprocess, "run", reject_subprocess)
    with pytest.raises(SystemExit) as exit_info:
        RUNNER.main()
    assert exit_info.value.code == 0

    output_root = paths["output_root"]
    assert {path.name for path in output_root.iterdir()} == {"not_applicable_receipt.json"}
    assert not (output_root / "run_receipt.json").exists()
    receipt = json.loads((output_root / "not_applicable_receipt.json").read_text())
    assert receipt["schema_name"].startswith(
        "intersubmod.matched_normal_candidate_control"
    )
    assert receipt["schema_version"] == "1.0.0"
    assert receipt["status"] == "NOT_APPLICABLE"
    assert receipt["execution_status"] == "NOT_APPLICABLE"
    assert receipt["reason"] == "ZERO_SELECTED_CANDIDATES"
    assert receipt["n_selected_candidates"] == 0
    assert receipt["selection_column"] == (
        "multi_marker_molecular_haplotype_base_candidate"
    )
    assert receipt["selection_value"] == "true"
    assert receipt["task_type"] == "B_comprehensive_validation"
    assert receipt["pass_semantics"] == (
        "receipt_integrity_only_not_normal_control_execution_or_negative_evidence"
    )
    assert receipt["pass"] is True
    assert receipt["source_authority"]["authority_id"] == "TEST_ONLY_UNSIGNED_AUTHORITY"
    assert receipt["code"] == RUNNER.capture_source_identity()
    assert receipt["source_lock"]["all_sources_read_only_and_unchanged"] is True
    assert receipt["command"][1] == "-I"
    assert oct((output_root / "not_applicable_receipt.json").stat().st_mode & 0o777) == "0o444"
    assert receipt["sample_outputs_created"] is False
    assert receipt["cpp_executed"] is False
    assert receipt["normal_control_executed"] is False
    assert receipt["not_evaluable_is_negative"] is False
    assert "not a negative" in receipt["interpretation"]

    expected_inputs = {
        key: path for key, path in paths.items() if key != "output_root"
    }
    assert set(receipt["inputs"]) == set(expected_inputs)
    for key, path in expected_inputs.items():
        artifact = receipt["inputs"][key]
        assert artifact["path"] == str(path.resolve())
        assert artifact["size_bytes"] == path.stat().st_size
        assert artifact["sha256"] == RUNNER.sha256(path)


def test_zero_selection_main_rejects_existing_output_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    args, paths = zero_selection_runner_args(tmp_path)
    output_root = paths["output_root"]
    output_root.mkdir()
    sentinel = output_root / "existing.txt"
    sentinel.write_text("keep\n", encoding="ascii")
    monkeypatch.setattr(sys, "argv", args)

    with pytest.raises(FileExistsError, match="Refusing to overwrite existing output root"):
        RUNNER.main()
    assert sentinel.read_text(encoding="ascii") == "keep\n"
    assert {path.name for path in output_root.iterdir()} == {"existing.txt"}


def test_normal_audit_gate_requires_all_seven_eligible() -> None:
    rows = [
        {
            "sample": sample,
            "bam": f"/{sample}.bam",
            "index": f"/{sample}.bam.bai",
            "normal_control_eligible": True,
        }
        for sample in RUNNER.DATASETS
    ]
    audit = {
        "schema_name": "intersubmod.matched_normal_methyl_tag_audit",
        "pass": True,
        "totals": {"n_samples": 7, "n_normal_control_eligible": 7},
        "samples": rows,
    }
    assert set(RUNNER.validate_normal_audit(audit, set(RUNNER.DATASETS))) == set(
        RUNNER.DATASETS
    )
    audit["samples"][3]["normal_control_eligible"] = False
    with pytest.raises(RuntimeError, match="ineligible"):
        RUNNER.validate_normal_audit(audit, set(RUNNER.DATASETS))


def test_runner_default_uses_versioned_normal_audit(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", ["runner", "--candidate-table", "candidates.tsv"])
    args = RUNNER.parse_args()
    assert args.normal_audit == RUNNER.CANONICAL_NORMAL_AUDIT


@pytest.mark.parametrize("module", [RUNNER, ANALYZER])
def test_formal_matched_normal_producers_reject_injected_argv(module: object) -> None:
    authority = {
        "authority_id": module.SOURCE_AUTHORITY.AUTHORITY_ID,
        "pass": True,
    }
    with pytest.raises(RuntimeError, match="direct-CLI canonical-process only"):
        module.resolve_release_command(["--injected"], authority)


@pytest.mark.parametrize("module", [RUNNER, ANALYZER])
def test_formal_matched_normal_producers_reject_process_missing_isolation(
    module: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    authority = {
        "authority_id": module.SOURCE_AUTHORITY.AUTHORITY_ID,
        "pass": True,
    }
    expected = module.canonical_task_b_command()
    assert expected[:4] == module.canonical_python_prefix()
    monkeypatch.setattr(sys, "argv", [expected[4], *expected[5:]])
    monkeypatch.setattr(
        module,
        "observed_process_command",
        lambda: [expected[0], expected[4], *expected[5:]],
    )
    with pytest.raises(RuntimeError, match="direct-CLI canonical-process only"):
        module.resolve_release_command(None, authority)


def test_overwrite_guards_reject_existing_roots(tmp_path: Path) -> None:
    existing = tmp_path / "existing"
    existing.mkdir()
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        RUNNER.create_output_root(existing)
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        ANALYZER.create_output_dir(existing)


def test_runner_artifact_counts_are_exact_and_region_errors_fail_closed(tmp_path: Path) -> None:
    output = tmp_path / "output"
    region = output / "sample" / "chr1" / "chr1_10" / "chr1_5_15"
    for relative in (
        "reads/reads.tsv",
        "methylation/methylation.csv",
        "distance/BERNOULLI/matrix.csv",
    ):
        path = region / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("header\n", encoding="ascii")
    (output / "significance_summary.csv").write_text("site\nchr1:10\n", encoding="ascii")
    (output / "significance_statistics.txt").write_text(
        "Total Regions Processed: 1\n", encoding="ascii"
    )
    (output / "region_stratification_status.tsv").write_text(
        "status\treason\nINSUFFICIENT_REGIONS\tBELOW_MIN_REGIONS\n", encoding="ascii"
    )
    log = tmp_path / "execution.log"
    log.write_text("[INFO] complete\n", encoding="ascii")
    candidate = RUNNER.Candidate("S1", "chr1", 10, "A", "G")

    validation = RUNNER.validate_output(output, [candidate], [log])
    assert validation["pass"] is True
    assert validation["exact_artifact_counts"] is True
    assert validation["region_stratification_status_accepted"] is True
    assert (
        validation["reads_files"],
        validation["methylation_files"],
        validation["bernoulli_matrix_files"],
    ) == (1, 1, 1)
    assert len(validation["artifact_set_sha256"]) == 64

    log.write_text("[ERROR] Region 1 (chr1:10) failed: synthetic\n", encoding="ascii")
    failed = RUNNER.validate_output(output, [candidate], [log])
    assert failed["pass"] is False
    assert failed["region_error_count"] == 1

    log.write_text("[INFO] complete\n", encoding="ascii")
    (output / "region_stratification_status.tsv").write_text(
        "status\treason\nFAILED\tSYNTHETIC\n", encoding="ascii"
    )
    failed_status = RUNNER.validate_output(output, [candidate], [log])
    assert failed_status["pass"] is False
    assert failed_status["region_stratification_status"] == "FAILED"
    assert failed_status["region_stratification_status_accepted"] is False


def make_read(
    read_id: str,
    is_tumor: bool,
    allele: str,
    start: int,
    hp: str = "poison-normal-hp",
) -> dict[str, str]:
    return {
        "read_id": read_id,
        "read_name": f"read-{read_id}",
        "chr": "chr1",
        "start": str(start),
        "end": str(start + 100),
        "mapq": "60",
        "hp": hp,
        "alt_support": allele,
        "is_tumor": "1" if is_tumor else "0",
        "strand": "+" if start % 2 == 0 else "-",
    }


def primary_core(row: dict[str, str], label: str) -> dict[str, object]:
    return {
        "read_name": row["read_name"],
        "chrom": row["chr"],
        "start": int(row["start"]),
        "end": int(row["end"]),
        "mapq": int(row["mapq"]),
        "strand": row["strand"],
        "label": label,
    }


def fake_phylo(distance: np.ndarray, methylation: np.ndarray, **_: object) -> dict[str, object]:
    midpoint = distance.shape[0] // 2
    labels = ["1-1"] * midpoint + ["1-2"] * (distance.shape[0] - midpoint)
    return {
        "coarse_ng": 2,
        "modal_fraction": 1.0,
        "unstable": False,
        "coarse_labels": labels,
    }


def fake_subset_factory(calls: list[list[str]]):
    def fake_subset(
        ids: list[str],
        _distance_ids: list[str],
        _distance: np.ndarray,
        _methylation_ids: list[str],
        _methylation: np.ndarray,
        _seed: int,
    ) -> dict[str, object]:
        calls.append(list(ids))
        evaluable = len(ids) >= 2 * ANALYZER.F.MIN_SIZE
        return {
            "status": "evaluable" if evaluable else "insufficient_reads",
            "evaluable": evaluable,
            "n_raw": len(ids),
            "n_matrix": len(ids) if evaluable else 0,
            "n_after_peel": len(ids) if evaluable else 0,
            "stable_multigroup": bool(evaluable and ids and ids[0].startswith("r")),
            "group_sizes": {"1": len(ids)} if evaluable else {},
        }

    return fake_subset


def test_normal_and_tumor_ref_subsets_are_separate_and_normal_hp_is_unused(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reads: dict[str, dict[str, str]] = {}
    normal_alleles = ["REF"] * 10 + ["ALT", "UNKNOWN"]
    for index, allele in enumerate(normal_alleles):
        reads[f"n{index}"] = make_read(f"n{index}", False, allele, 100 + index)
    for index in range(6):
        reads[f"r{index}"] = make_read(f"r{index}", True, "REF", 200 + index)
    for index in range(3):
        reads[f"a{index}"] = make_read(f"a{index}", True, "ALT", 300 + index)
    identifiers = list(reads)
    distance = np.zeros((len(identifiers), len(identifiers)), dtype=float)
    methylation = np.tile(np.asarray([0.05, 0.10, 0.90, 0.95]), (len(identifiers), 1))
    assignment = {
        "core_reads": [
            primary_core(reads["a0"], "group-1"),
            primary_core(reads["a1"], "group-1"),
            primary_core(reads["a2"], "group-2"),
        ]
    }
    candidate = {"sample": "S1", "chrom": "chr1", "pos": 123, "ref": "A", "alt": "G"}
    monkeypatch.setattr(ANALYZER.F, "analyze_phylo", fake_phylo)

    first = ANALYZER.analyze_site_from_data(
        candidate, reads, identifiers, distance, identifiers, methylation, assignment
    )
    changed_hp = deepcopy(reads)
    for read_id in [value for value in changed_hp if value.startswith("n")]:
        changed_hp[read_id]["hp"] = f"different-{read_id}"
    second = ANALYZER.analyze_site_from_data(
        candidate, changed_hp, identifiers, distance, identifiers, methylation, assignment
    )

    assert first == second
    assert first["normal_n_raw"] == 10
    assert first["normal_ref_methyl_n_raw"] == 10
    assert first["tumor_ref_n_raw"] == 6
    assert first["normal_stable_multigroup"] is True
    assert first["normal_ref_methyl_stable_multigroup"] is True
    assert first["tumor_ref_stable_multigroup"] is True
    assert first["normal_genetic_alt_support_count"] == 1
    assert first["normal_focal_allele_counts"] == {"ALT": 1, "REF": 10, "UNKNOWN": 1}
    assert first["tumor_focal_allele_counts"] == {"ALT": 3, "REF": 6, "UNKNOWN": 0}
    assert first["normal_called_reads"] == 11
    assert first["normal_focal_callability_gate"] is True
    assert first["normal_control_evaluable"] is True
    assert first["normal_genetic_alt_absence_gate"] is False
    assert first["normal_ref_methyl_nonreplication_gate"] is False
    assert first["normal_methyl_background_population"] == (
        "is_tumor=0_and_focal_call=REF_only"
    )
    assert first["normal_hp_used"] is False
    assert first["normal_hp_policy"] == "PROHIBITED_NOT_USED"
    assert first["primary_group_assignment_coverage"] == 1.0


def test_all_unknown_normal_reads_never_pass_callability_or_evaluability(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reads = {
        f"n{index}": make_read(f"n{index}", False, "UNKNOWN", 100 + index)
        for index in range(12)
    }
    reads["a0"] = make_read("a0", True, "ALT", 300)
    assignment = {"core_reads": [primary_core(reads["a0"], "group-1")]}
    identifiers = list(reads)
    distance = np.zeros((len(identifiers), len(identifiers)), dtype=float)
    methylation = np.zeros((len(identifiers), 4), dtype=float)
    calls: list[list[str]] = []
    monkeypatch.setattr(ANALYZER, "analyze_subset", fake_subset_factory(calls))

    row = ANALYZER.analyze_site_from_data(
        {"sample": "S1", "chrom": "chr1", "pos": 123, "ref": "A", "alt": "G"},
        reads,
        identifiers,
        distance,
        identifiers,
        methylation,
        assignment,
    )
    assert calls[0] == []
    assert row["normal_called_reads"] == 0
    assert row["normal_ref_reads"] == 0
    assert row["normal_alt_reads"] == 0
    assert row["normal_unknown_reads"] == 12
    assert row["normal_focal_callability_gate"] is False
    assert row["normal_control_evaluable"] is False
    assert row["normal_ref_methyl_nonreplication_gate"] is None
    assert row["normal_control_status"] == "NOT_EVALUABLE"
    assert row["normal_control_not_evaluable_reason"] == "NORMAL_CALLED_READS_LT_10"


def test_callability_does_not_convert_insufficient_ref_methyl_depth_into_negative(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reads: dict[str, dict[str, str]] = {}
    for index in range(5):
        reads[f"r{index}"] = make_read(f"r{index}", False, "REF", 100 + index)
        reads[f"n{index}"] = make_read(f"n{index}", False, "ALT", 120 + index)
    reads["a0"] = make_read("a0", True, "ALT", 300)
    assignment = {"core_reads": [primary_core(reads["a0"], "group-1")]}
    identifiers = list(reads)
    distance = np.zeros((len(identifiers), len(identifiers)), dtype=float)
    methylation = np.zeros((len(identifiers), 4), dtype=float)
    monkeypatch.setattr(ANALYZER, "analyze_subset", fake_subset_factory([]))

    row = ANALYZER.analyze_site_from_data(
        {"sample": "S1", "chrom": "chr1", "pos": 123, "ref": "A", "alt": "G"},
        reads,
        identifiers,
        distance,
        identifiers,
        methylation,
        assignment,
    )
    assert row["normal_called_reads"] == 10
    assert row["normal_ref_reads"] == 5
    assert row["normal_focal_callability_gate"] is True
    assert row["normal_ref_methyl_evaluable"] is False
    assert row["normal_control_evaluable"] is False
    assert row["normal_ref_methyl_nonreplication_gate"] is None
    assert row["normal_control_not_evaluable_reason"] == (
        "NORMAL_REF_METHYL_INSUFFICIENT_READS"
    )


def test_summary_keeps_not_evaluable_separate_from_negative(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reads = {
        f"n{index}": make_read(f"n{index}", False, "UNKNOWN", 100 + index)
        for index in range(12)
    }
    reads["a0"] = make_read("a0", True, "ALT", 300)
    identifiers = list(reads)
    assignment = {"core_reads": [primary_core(reads["a0"], "group-1")]}
    monkeypatch.setattr(ANALYZER, "analyze_subset", fake_subset_factory([]))
    row = ANALYZER.analyze_site_from_data(
        {"sample": "S1", "chrom": "chr1", "pos": 123, "ref": "A", "alt": "G"},
        reads,
        identifiers,
        np.zeros((len(identifiers), len(identifiers)), dtype=float),
        identifiers,
        np.zeros((len(identifiers), 4), dtype=float),
        assignment,
    )

    summary = ANALYZER.summarize([row])
    assert summary["n_normal_control_evaluable"] == 0
    assert summary["n_normal_control_not_evaluable"] == 1
    assert summary["n_normal_ref_methyl_nonreplicating"] == 0


def test_normal_negative_guardrail_explicitly_preserves_asm_and_cell_state() -> None:
    assert "allele-specific methylation (ASM)" in ANALYZER.GUARDRAIL
    assert "cell-state-dependent methylation" in ANALYZER.GUARDRAIL
    assert "does not prove" in ANALYZER.GUARDRAIL


def test_primary_identity_collision_and_missing_are_hard_failures() -> None:
    first = make_read("a", True, "ALT", 100)
    duplicate = make_read("b", True, "ALT", 100)
    duplicate["read_name"] = first["read_name"]
    duplicate["strand"] = first["strand"]
    core = [primary_core(first, "group-1")]
    with pytest.raises(RuntimeError, match="identity collision"):
        ANALYZER.join_primary_group_assignments(core, {"a": first, "b": duplicate})
    with pytest.raises(RuntimeError, match="Missing paired tumor identities"):
        ANALYZER.join_primary_group_assignments(core, {})


def test_primary_identity_distinguishes_mapq() -> None:
    primary = make_read("a", True, "ALT", 100)
    paired = make_read("b", True, "ALT", 100)
    paired["read_name"] = primary["read_name"]
    paired["strand"] = primary["strand"]
    paired["mapq"] = "59"
    with pytest.raises(RuntimeError, match="Missing paired tumor identities"):
        ANALYZER.join_primary_group_assignments(
            [primary_core(primary, "group-1")], {"b": paired}
        )


def test_primary_assignment_parser_reads_gzip_and_rejects_ref_alt_drift(tmp_path: Path) -> None:
    path = tmp_path / "assignments.jsonl.gz"
    assignment = {
        "sample": "S1",
        "chrom": "chr1",
        "pos": 10,
        "posthoc": {"ref": "A", "alt": "G"},
        "core_reads": [{"read_name": "r", "chrom": "chr1", "start": 1, "end": 2, "strand": "+", "label": "1"}],
    }
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps(assignment) + "\n")
    expected = {("S1", "chr1", 10, "A", "G")}
    assert set(ANALYZER.load_primary_assignments(path, expected)) == expected
    with pytest.raises(RuntimeError, match="REF/ALT mismatch"):
        ANALYZER.load_primary_assignments(path, {("S1", "chr1", 10, "A", "T")})


def test_paired_analyzer_recomputes_runner_artifact_set_digest(tmp_path: Path) -> None:
    paired_root = tmp_path / "paired"
    output_dir = paired_root / "S1" / "output"
    region = output_dir / "chr1_100" / "region_1"
    artifact_paths = [
        region / "reads" / "reads.tsv",
        region / "distance" / "BERNOULLI" / "matrix.csv",
        region / "methylation" / "methylation.csv",
        output_dir / "significance_summary.csv",
        output_dir / "significance_statistics.txt",
        output_dir / "region_stratification_status.tsv",
    ]
    for path in artifact_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"{path.name}\n", encoding="utf-8")
    digest = ANALYZER.file_set_sha256(artifact_paths, output_dir)
    receipt = {
        "schema_name": "intersubmod.matched_normal_candidate_control_run",
        "schema_version": "1.0.0",
        "receipts": [
            {
                "sample": "S1",
                "pass": True,
                "outputs": {"output_dir": str(output_dir.resolve())},
                "candidate_sites": [
                    {
                        "sample": "S1",
                        "chrom": "chr1",
                        "pos": 100,
                        "ref": "A",
                        "alt": "C",
                    }
                ],
                "validation": {
                    "exact_artifact_counts": True,
                    "artifact_set_sha256": digest,
                },
            }
        ],
        "pass": True,
    }
    paired_root.mkdir(parents=True, exist_ok=True)
    (paired_root / "run_receipt.json").write_text(
        json.dumps(receipt), encoding="utf-8"
    )

    tasks, _, _, validations = ANALYZER.load_paired_tasks(paired_root)
    assert len(tasks) == 1
    assert validations["S1"]["artifact_set_sha256"] == digest
    assert validations["S1"]["artifacts_verified"] == 6

    artifact_paths[0].write_text("tampered\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="artifact_set_sha256 mismatch"):
        ANALYZER.load_paired_tasks(paired_root)
