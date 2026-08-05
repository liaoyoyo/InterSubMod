from pathlib import Path
import csv
import gzip
import sys

import pytest


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from summarize_hcc1395_full import (  # noqa: E402
    ContractError,
    _verify_extract_child,
    block_evidence_gate,
    file_identity,
    parse_gnu_time_v,
    validate_outer_time_contract,
    validate_partition_internal_timings,
    validate_partition_output,
    verify_sha256_sidecar,
    write_summary_outputs,
)


CHR22_PROBE = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260718_k_gt8_read_supported_segmentation/probes/HCC1395_chr22/"
    "partition_v2"
)
CHR22_EXTRACTION_PROBE = CHR22_PROBE.parent / "extraction_v2"
CHR6_FULL_STAGE = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260718_k_gt8_read_supported_segmentation/full/HCC1395_chr1_22_v1/"
    "chromosomes/chr6"
)


@pytest.mark.parametrize(
    ("patterns", "weight", "active_sites", "expected"),
    [
        (1, 1, 2, "TREE_READY_LOCAL"),
        (0, 0, 2, "ABSTAIN_ZERO_LOCAL_EVIDENCE"),
        (1, 1, 1, "ABSTAIN_ZERO_LOCAL_EVIDENCE"),
        (0, 0, 0, "ABSTAIN_ZERO_LOCAL_EVIDENCE"),
    ],
)
def test_block_evidence_gate_is_exact_three_condition_and(
    patterns, weight, active_sites, expected
):
    assert block_evidence_gate(patterns, weight, active_sites) == expected


@pytest.mark.parametrize(
    "timings",
    [
        {},
        {"load_and_aggregate_primary_patterns": 1.0},
        {
            "load_and_aggregate_primary_patterns": -1.0,
            "ordered_hypergraph_dp": 0.1,
        },
        {
            "load_and_aggregate_primary_patterns": 1.0,
            "ordered_hypergraph_dp": float("nan"),
        },
    ],
)
def test_partition_internal_timing_contract_fails_closed(timings):
    with pytest.raises(ContractError, match="timing"):
        validate_partition_internal_timings(
            {"timings_seconds": timings},
            "chrFixture",
        )


def test_chr22_authenticated_partition_probe_dry_run():
    if not CHR22_PROBE.is_dir():
        pytest.skip(f"authenticated chr22 probe is unavailable: {CHR22_PROBE}")

    metrics, component_rows, evidence = validate_partition_output(
        CHR22_PROBE,
        "chr22",
        {
            "ssnv_sites": 543,
            "k_gt8_components": 6,
            "k_gt8_sites": 98,
        },
    )

    assert len(component_rows) == 6
    assert metrics["components"] == 6
    assert metrics["sites"] == 98
    assert metrics["old_selected_sites"] == 48
    assert metrics["old_excluded_sites"] == 50
    assert metrics["new_blocks"] == 16
    assert metrics["new_retained_sites"] == 98
    assert metrics["primary_active_sites"] == 97
    assert metrics["exact_patterns"] == 371
    assert metrics["raw_total_molecule_weight"] == 1151
    assert metrics["new_retained_molecule_weight"] == 1037
    assert metrics["new_lost_molecule_weight"] == 114
    assert metrics["old_densest8_retained_molecule_weight"] == 787
    assert metrics["unavoidable_patterns"] == 87
    assert metrics["unavoidable_molecule_weight"] == 93
    assert metrics["unavoidable_n_fixed_ra_gt8_patterns"] == 64
    assert metrics["unavoidable_n_fixed_ra_gt8_molecule_weight"] == 68
    assert metrics["unavoidable_n_fixed_ra_lte8_span_gt8_patterns"] == 23
    assert (
        metrics["unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight"]
        == 25
    )
    assert (
        metrics["unavoidable_n_fixed_ra_gt8_patterns"]
        + metrics["unavoidable_n_fixed_ra_lte8_span_gt8_patterns"]
        == metrics["unavoidable_patterns"]
    )
    assert metrics["partition_pattern_load_aggregate_seconds"] == pytest.approx(
        1.2350103419739753
    )
    assert metrics["partition_ordered_hypergraph_dp_seconds"] == pytest.approx(
        0.1026568750385195
    )
    assert metrics["weight_stable_components"] == 4
    assert metrics["zero_evidence_blocks"] == 4
    assert metrics["zero_evidence_block_sites"] == 25
    assert metrics["tree_ready_blocks"] == 12
    assert metrics["tree_ready_block_sites"] == 73
    assert metrics["abstain_blocks"] == 4
    assert metrics["abstain_block_sites"] == 25
    assert metrics["tree_ready_blocks"] + metrics["abstain_blocks"] == 16
    assert metrics["tree_ready_block_sites"] + metrics["abstain_block_sites"] == 98
    assert metrics["block_span_bp_distribution"]["median"] == "51937.50"
    assert (
        metrics["component_retention_delta_distribution"]["median"]
        == "0.1946428571428571428571428572"
    )
    assert metrics["block_k_distribution"] == {
        "1": 1,
        "2": 2,
        "4": 2,
        "6": 1,
        "7": 1,
        "8": 9,
    }
    assert evidence["partition_receipt"]["sha256"]
    assert set(evidence["outputs"]) == {
        "legacy_components",
        "blocks",
        "site_membership",
        "cut_constraints",
    }
    gates = [row["evidence_gate"] for row in evidence["_block_all_rows"]]
    assert gates.count("TREE_READY_LOCAL") == 12
    assert gates.count("ABSTAIN_ZERO_LOCAL_EVIDENCE") == 4


def test_chr22_authenticated_extraction_probe_diagnostics():
    receipt = CHR22_EXTRACTION_PROBE / "receipt.json"
    if not receipt.is_file():
        pytest.skip(
            f"authenticated chr22 extraction probe is unavailable: {receipt}"
        )

    verified = _verify_extract_child(
        CHR22_EXTRACTION_PROBE,
        chrom="chr22",
        inventory={"ssnv_sites": 543},
        stage_receipt={"child_receipt": file_identity(receipt)},
    )

    assert verified["diagnostics"] == {
        "raw_overlapping_alignments": 69262,
        "eligible_alignments_pre_identity_collapse": 63424,
        "duplicate_alignment_identities_collapsed": 9660,
        "canonical_unique_molecules": 53764,
        "molecule_sparse_rows": 53764,
        "sparse_site_calls": 70163,
        "known_ps_hp12_molecule_rows": 44531,
        "known_ps_hp1_molecule_rows": 25304,
        "known_ps_hp2_molecule_rows": 19227,
        "extraction_failed_checks": 0,
        "duplicate_identity_conflicts": 0,
    }


def test_chr6_full_stage_abstains_on_3207_zero_evidence_blocks():
    partition = CHR6_FULL_STAGE / "partition"
    extraction = CHR6_FULL_STAGE / "extract"
    if not (partition / "receipt.json").is_file():
        pytest.skip(f"authenticated chr6 full stage is unavailable: {partition}")

    metrics, _, evidence = validate_partition_output(
        partition,
        "chr6",
        {
            "ssnv_sites": 27099,
            "k_gt8_components": 83,
            "k_gt8_sites": 25657,
        },
        expected_input_paths={
            "site_catalog": extraction / "HCC1395.chr6.site_catalog.tsv.gz",
            "molecule_calls": (
                extraction / "HCC1395.chr6.molecule_sparse_calls.tsv.gz"
            ),
        },
    )

    assert metrics["new_blocks"] == 3252
    assert metrics["tree_ready_blocks"] == 45
    assert metrics["tree_ready_block_sites"] == 310
    assert metrics["abstain_blocks"] == 3207
    assert metrics["abstain_block_sites"] == 25347
    assert metrics["zero_evidence_blocks"] == 3207
    assert metrics["unavoidable_n_fixed_ra_gt8_patterns"] == 827
    assert metrics["unavoidable_n_fixed_ra_gt8_molecule_weight"] == 845
    assert metrics["unavoidable_n_fixed_ra_lte8_span_gt8_patterns"] == 150
    assert (
        metrics["unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight"]
        == 154
    )
    assert metrics["tree_ready_blocks"] + metrics["abstain_blocks"] == 3252
    assert (
        metrics["tree_ready_block_sites"] + metrics["abstain_block_sites"]
        == 25657
    )
    assert (
        sum(
            row["evidence_gate"] == "ABSTAIN_ZERO_LOCAL_EVIDENCE"
            for row in evidence["_block_all_rows"]
        )
        == 3207
    )


def test_chr22_probe_fails_closed_on_inventory_drift():
    if not CHR22_PROBE.is_dir():
        pytest.skip(f"authenticated chr22 probe is unavailable: {CHR22_PROBE}")

    with pytest.raises(ContractError, match="canonical target mismatch"):
        validate_partition_output(
            CHR22_PROBE,
            "chr22",
            {
                "ssnv_sites": 543,
                "k_gt8_components": 7,
                "k_gt8_sites": 98,
            },
        )


def test_probe_summary_outputs_are_sha_bound_and_never_overwritten(tmp_path):
    if not CHR22_PROBE.is_dir():
        pytest.skip(f"authenticated chr22 probe is unavailable: {CHR22_PROBE}")
    metrics, component_rows, evidence = validate_partition_output(
        CHR22_PROBE,
        "chr22",
        {
            "ssnv_sites": 543,
            "k_gt8_components": 6,
            "k_gt8_sites": 98,
        },
    )
    metrics.pop("_block_span_bp_values")
    metrics.pop("_component_retention_delta_values")
    metrics["partition_stage_status"] = "PROBE"
    payload = {
        "per_chromosome": [metrics],
        "totals": {**metrics, "chrom": "ALL"},
    }
    output_dir = tmp_path / "summary"

    block_rows = evidence["_block_all_rows"]
    outputs = write_summary_outputs(
        output_dir, payload, component_rows, block_rows
    )

    assert set(outputs) == {
        "summary.json",
        "summary.tsv",
        "component_all.tsv.gz",
        "block_all.tsv.gz",
    }
    for name in outputs:
        assert verify_sha256_sidecar(output_dir / name) == outputs[name]["sha256"]
    with (output_dir / "summary.tsv").open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        written_summary = list(csv.DictReader(handle, delimiter="\t"))
    assert written_summary[0]["tree_ready_blocks"] == "12"
    assert written_summary[0]["tree_ready_block_sites"] == "73"
    assert written_summary[0]["abstain_blocks"] == "4"
    assert written_summary[0]["abstain_block_sites"] == "25"
    with gzip.open(
        output_dir / "block_all.tsv.gz",
        "rt",
        encoding="utf-8",
        newline="",
    ) as handle:
        written_blocks = list(csv.DictReader(handle, delimiter="\t"))
    assert len(written_blocks) == 16
    assert sum(row["evidence_gate"] == "TREE_READY_LOCAL" for row in written_blocks) == 12
    assert (
        sum(
            row["evidence_gate"] == "ABSTAIN_ZERO_LOCAL_EVIDENCE"
            for row in written_blocks
        )
        == 4
    )
    with pytest.raises(ContractError, match="refusing to overwrite"):
        write_summary_outputs(output_dir, payload, component_rows, block_rows)


def test_parse_gnu_time_v_supports_hour_elapsed_format(tmp_path):
    resource = tmp_path / "time_v.txt"
    resource.write_text(
        "\n".join(
            [
                'Command being timed: "example"',
                "User time (seconds): 10.25",
                "System time (seconds): 1.75",
                "Percent of CPU this job got: 50%",
                "Elapsed (wall clock) time (h:mm:ss or m:ss): 1:02:03.50",
                "Maximum resident set size (kbytes): 4096",
                "File system inputs: 8",
                "File system outputs: 16",
                "Exit status: 0",
                "",
            ]
        ),
        encoding="utf-8",
    )

    parsed = parse_gnu_time_v(resource)

    assert parsed["elapsed_seconds"] == 3723.5
    assert parsed["max_rss_kb"] == 4096
    assert parsed["exit_status"] == 0
    assert parsed["command_text"] == "example"
    assert parsed["command_argv"] == ["example"]


def test_parse_gnu_time_v_requires_command_identity(tmp_path):
    resource = tmp_path / "time_v_without_command.txt"
    resource.write_text(
        "\n".join(
            [
                "User time (seconds): 1",
                "System time (seconds): 0",
                "Percent of CPU this job got: 50%",
                "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.00",
                "Maximum resident set size (kbytes): 4096",
                "File system inputs: 0",
                "File system outputs: 0",
                "Exit status: 0",
            ]
        ),
        encoding="utf-8",
    )

    with pytest.raises(ContractError, match="Command being timed"):
        parse_gnu_time_v(resource)


def _outer_contract_fixture(tmp_path):
    runner = tmp_path / "run_hcc1395_full_segmentation.py"
    manifest = tmp_path / "canonical_manifest.json"
    output_root = tmp_path / "HCC1395_chr1_22_v1"
    argv = [
        "/usr/bin/nice",
        "-n",
        "10",
        "/usr/bin/python3",
        str(runner),
        "--manifest",
        str(manifest),
        "--output-root",
        str(output_root),
    ]
    return runner, manifest, output_root, argv


def test_outer_time_contract_binds_fresh_command_and_reports_overhead(tmp_path):
    runner, manifest, output_root, argv = _outer_contract_fixture(tmp_path)

    result = validate_outer_time_contract(
        {"command_argv": argv, "elapsed_seconds": 120.0},
        full_root=output_root,
        runner_path=runner,
        manifest_path=manifest,
        sequential_stage_wall_seconds=100.0,
    )

    assert result["command_binding_verified"] is True
    assert result["fresh_non_resume_command_verified"] is True
    assert result["sequential_stage_wall_seconds"] == 100.0
    assert result["outer_minus_sequential_stage_wall_seconds"] == 20.0
    assert result["runner_overhead_seconds"] == 20.0


@pytest.mark.parametrize(
    ("mutation", "error"),
    [
        ("runner", "contracted runner"),
        ("manifest", "--manifest does not match"),
        ("output_root", "--output-root does not match"),
        ("resume", "without --resume"),
        ("duplicate_manifest", "exactly one --manifest"),
    ],
)
def test_outer_time_contract_fails_closed_on_command_drift(
    tmp_path, mutation, error
):
    runner, manifest, output_root, argv = _outer_contract_fixture(tmp_path)
    if mutation == "runner":
        argv[4] = str(tmp_path / "wrong_runner.py")
    elif mutation == "manifest":
        argv[6] = str(tmp_path / "wrong_manifest.json")
    elif mutation == "output_root":
        argv[8] = str(tmp_path / "wrong_output")
    elif mutation == "resume":
        argv.append("--resume")
    elif mutation == "duplicate_manifest":
        argv.extend(["--manifest", str(manifest)])

    with pytest.raises(ContractError, match=error):
        validate_outer_time_contract(
            {"command_argv": argv, "elapsed_seconds": 120.0},
            full_root=output_root,
            runner_path=runner,
            manifest_path=manifest,
            sequential_stage_wall_seconds=100.0,
        )


def test_outer_time_contract_fails_when_elapsed_is_below_stage_sum(tmp_path):
    runner, manifest, output_root, argv = _outer_contract_fixture(tmp_path)

    with pytest.raises(ContractError, match="shorter than"):
        validate_outer_time_contract(
            {"command_argv": argv, "elapsed_seconds": 99.70},
            full_root=output_root,
            runner_path=runner,
            manifest_path=manifest,
            sequential_stage_wall_seconds=100.0,
        )


def test_outer_time_contract_allows_small_gnu_display_rounding_gap(tmp_path):
    runner, manifest, output_root, argv = _outer_contract_fixture(tmp_path)

    result = validate_outer_time_contract(
        {"command_argv": argv, "elapsed_seconds": 99.80},
        full_root=output_root,
        runner_path=runner,
        manifest_path=manifest,
        sequential_stage_wall_seconds=100.0,
    )

    assert result["outer_minus_sequential_stage_wall_seconds"] == pytest.approx(-0.2)
    assert result["runner_overhead_seconds"] == 0.0
