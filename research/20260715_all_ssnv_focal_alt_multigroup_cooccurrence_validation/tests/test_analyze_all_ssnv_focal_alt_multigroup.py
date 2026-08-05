from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pysam
import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "analyze_all_ssnv_focal_alt_multigroup.py"
)
SPEC = importlib.util.spec_from_file_location("analyze_all_ssnv_focal_alt_multigroup", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_gzip_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def make_sidecar(tmp_path: Path, rows: list[str]) -> Path:
    plain = tmp_path / "tags.tsv"
    plain.write_text(
        "#CHROM\tSTART0\tEND0\tQNAME\tFLAG\tMAPQ\tCIGAR_B2\tHP\tPS\n"
        + "\n".join(rows)
        + "\n",
        encoding="ascii",
    )
    compressed = tmp_path / "tags.tsv.gz"
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(
        str(compressed),
        seq_col=0,
        start_col=1,
        end_col=2,
        zerobased=True,
        force=True,
    )
    return compressed


def make_reads(tmp_path: Path, hp: str = "2") -> Path:
    path = tmp_path / "reads.tsv"
    path.write_text(
        "read_id\tread_name\tchr\tstart\tend\tmapq\thp\talt_support\tis_tumor\tstrand\n"
        f"0\tread-a\tchr1\t90\t150\t60\t{hp}\tALT\t1\t+\n",
        encoding="ascii",
    )
    return path


def test_truth_and_ledger_join_are_posthoc(tmp_path: Path) -> None:
    truth_path = tmp_path / "truth.tsv.gz"
    write_gzip_tsv(
        truth_path,
        ["sample", "chrom", "pos", "ref", "alt", "truth_label"],
        [
            {
                "sample": "S1",
                "chrom": "chr1",
                "pos": 100,
                "ref": "A",
                "alt": "G",
                "truth_label": "TP",
            }
        ],
    )
    ledger_path = tmp_path / "ledger.tsv.gz"
    ledger_fields = [
        "sample",
        "chrom",
        "pos",
        "ref",
        "alt",
        "longphase_recalibrated_filter",
        "caller_filter",
        "longphase_filter_transition",
        "ssnv_branch",
        "selected_for_read_census",
        "component_id",
        "component_size",
        "selected_group_id",
        "selected_group_size",
        "pooled_ref_reads",
        "pooled_alt_reads",
        "family_coverage_json",
    ]
    write_gzip_tsv(
        ledger_path,
        ledger_fields,
        [
            {
                "sample": "S1",
                "chrom": "chr1",
                "pos": 100,
                "ref": "A",
                "alt": "G",
                "longphase_recalibrated_filter": "PASS",
                "caller_filter": "LowQual",
                "longphase_filter_transition": "LowQual->PASS",
                "ssnv_branch": "retained",
                "selected_for_read_census": "true",
                "component_id": "component-7",
                "component_size": 9,
                "selected_group_id": "group-2",
                "selected_group_size": 4,
                "pooled_ref_reads": 12,
                "pooled_alt_reads": 8,
                "family_coverage_json": json.dumps({"1-1": [5, 3], "2-1": [4, 0]}),
            }
        ],
    )

    truth = MODULE.load_truth_labels(truth_path, "S1")
    ledger = MODULE.load_ledger_context(ledger_path, "S1")
    posthoc = MODULE.join_posthoc_metadata(
        {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G", "filter": "PASS"},
        truth,
        ledger,
        "BIO-1",
    )

    assert posthoc["truth_label"] == "TP"
    assert posthoc["biological_id"] == "BIO-1"
    assert posthoc["ssnv_branch"] == "retained"
    assert posthoc["component_id"] == "component-7"
    assert posthoc["selected_group_id"] == "group-2"
    assert json.loads(posthoc["layered_alt_families"]) == ["1-1"]


def test_latest_tags_replace_reads_tsv_hp_and_add_ps(tmp_path: Path) -> None:
    sidecar = make_sidecar(tmp_path, ["chr1\t90\t150\tread-a\t0\t60\tcigar-a\t1-2\t101"])
    try:
        reads, stats = MODULE.load_latest_joined_reads(
            make_reads(tmp_path, hp="2"),
            str(sidecar),
            str(sidecar) + ".tbi",
            "chr1",
            100,
        )
    finally:
        MODULE.close_sidecar_cache()

    assert reads["0"]["source_hp"] == "2"
    assert reads["0"]["hp"] == "1-2"
    assert reads["0"]["ps"] == 101
    assert reads["0"]["latest_tag_full_identity_count"] == 1
    assert stats == {
        "latest_tag_rows_fetched": 1,
        "latest_tag_rows_eligible": 1,
        "latest_tag_reads_joined": 1,
        "latest_tag_ps_present": 1,
        "latest_tag_projection_multimatch_reads": 0,
        "source_hp_replaced_reads": 1,
    }


def test_latest_tag_projection_requires_one_full_alignment_identity(tmp_path: Path) -> None:
    sidecar = make_sidecar(
        tmp_path,
        [
            "chr1\t90\t150\tread-a\t0\t60\tcigar-a\t1-2\t101",
            "chr1\t90\t150\tread-a\t0\t60\tcigar-b\t1-2\t101",
        ],
    )
    try:
        with pytest.raises(RuntimeError, match="not a unique full alignment identity"):
            MODULE.load_latest_joined_reads(
                make_reads(tmp_path, hp="2"),
                str(sidecar),
                str(sidecar) + ".tbi",
                "chr1",
                100,
            )
    finally:
        MODULE.close_sidecar_cache()


def test_screen_summary_terminal_latest_hp_ps_join_audit() -> None:
    summary = MODULE.ScreenSummary()

    def row(*, joined: int, ps: int, replaced: int) -> dict[str, object]:
        return {
            "analysis_status": "evaluable",
            "evidence_tier": "E1_EVALUABLE_OR_SILHOUETTE_ONLY",
            "ssnv_branch": "retained",
            "stable_null_multigroup": False,
            "residual_unexplained_multigroup": False,
            "phase_anchored_robust_epigenetic_candidate": False,
            "strict_confirm_candidate": False,
            "strict_confirm_status": "NOT_RUN",
            "latest_tag_join_status": "PASS",
            "latest_tag_rows_fetched": joined + 2,
            "latest_tag_rows_eligible": joined + 1,
            "latest_tag_reads_joined": joined,
            "latest_tag_ps_present": ps,
            "latest_tag_projection_multimatch_reads": 0,
            "source_hp_replaced_reads": replaced,
            "n_reads_total": joined,
        }

    summary.add(row(joined=3, ps=2, replaced=1))
    summary.add(row(joined=0, ps=0, replaced=0))
    audit = summary.latest_tag_payload()

    assert audit == {
        "authoritative_tag_source": "same_run_LongPhase_S_external_HP_PS_sidecar",
        "embedded_reads_tsv_hp_used_for_analysis": False,
        "join_occurs_before_focal_ALT_selection": True,
        "counting_unit": "site_read_observation_not_globally_unique_read",
        "n_sites": 2,
        "join_status_counts": {"PASS": 2},
        "all_sites_pass": True,
        "n_reads_tsv_site_rows": 3,
        "n_exact_hp_ps_site_read_joins": 3,
        "every_reads_tsv_row_joined": True,
        "n_ps_present_site_read_joins": 2,
        "n_source_hp_replaced_site_read_joins": 1,
        "n_sidecar_rows_fetched_site_observations": 7,
        "n_sidecar_rows_eligible_site_observations": 5,
        "n_projection_multimatch_site_reads": 0,
        "all_projection_identities_unique": True,
        "n_sites_with_zero_reads_tsv_rows": 1,
        "pass": True,
    }

    bad = row(joined=2, ps=0, replaced=0)
    bad["n_reads_total"] = 3
    with pytest.raises(RuntimeError, match="consume every reads.tsv row"):
        MODULE.ScreenSummary().add(bad)


def test_focal_key_contract_accepts_longphase_pass_ledger_superset() -> None:
    focal_a = ("chr1", 100, "A", "C")
    focal_b = ("chr2", 200, "G", "T")
    ledger_extra = ("chrX", 300, "C", "A")
    result = MODULE.validate_focal_key_contract(
        "S1",
        2,
        {focal_a, focal_b},
        {focal_a: "TP", focal_b: "FP"},
        {focal_a: {}, focal_b: {}, ledger_extra: {}},
        {("chr1", 100), ("chr2", 200)},
        {("chr1", 100): {}, ("chr2", 200): {}},
    )

    assert result["longphase_pass_ledger_keys"] == 3
    assert result["longphase_pass_ledger_extra_context_keys"] == 1
    assert result["all_focal_keys_in_longphase_pass_ledger"] is True
    assert result["pass"] is True


def test_focal_key_contract_rejects_missing_longphase_pass_key() -> None:
    focal = ("chr1", 100, "A", "C")
    with pytest.raises(RuntimeError, match="missing from the LongPhase-S PASS ledger"):
        MODULE.validate_focal_key_contract(
            "S1",
            1,
            {focal},
            {focal: "TP"},
            {},
            {("chr1", 100)},
            {("chr1", 100): {}},
        )


def test_clustering_receives_no_truth_or_ledger_fields(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: list[dict[str, object]] = []

    def fake_screen(screen: dict[str, object]):
        captured.append(dict(screen))
        return {
            "sample": screen["sample"],
            "chrom": screen["chrom"],
            "pos": screen["pos"],
            "analysis_status": "insufficient_alt_reads",
            "screen_token": "identical",
        }, None

    monkeypatch.setattr(MODULE, "screen_site", fake_screen)
    screen = {"sample": "S1", "chrom": "chr1", "pos": 100, "reads_path": "/synthetic"}
    common = {
        "biological_id": "BIO-1",
        "chrom": "chr1",
        "pos": 100,
        "ref": "A",
        "alt": "G",
        "ssnv_branch": "retained",
        "component_id": "component-7",
        "selected_group_id": "group-2",
    }
    first, _ = MODULE.process_site_task(
        {"screen": dict(screen), "posthoc": {**common, "truth_label": "TP"}}
    )
    second, _ = MODULE.process_site_task(
        {"screen": dict(screen), "posthoc": {**common, "truth_label": "FP"}}
    )

    assert captured == [screen, screen]
    assert not MODULE.CLUSTERING_FORBIDDEN_FIELDS.intersection(captured[0])
    assert first["screen_token"] == second["screen_token"] == "identical"
    assert (first["truth_label"], second["truth_label"]) == ("TP", "FP")


def test_clustering_hard_rejects_label_leakage() -> None:
    with pytest.raises(RuntimeError, match="forbidden in clustering input"):
        MODULE.reject_clustering_label_leakage(
            {"sample": "S1", "chrom": "chr1", "pos": 100, "truth_label": "TP"}
        )


@pytest.mark.parametrize(
    ("coarse_ng", "unstable", "ari_min", "expected"),
    [
        (2, False, 0.80, True),
        (3, False, 1.00, True),
        (1, False, 1.00, False),
        (2, True, 1.00, False),
        (2, False, 0.799999, False),
        (2, False, None, False),
        (2, False, float("nan"), False),
    ],
)
def test_m1_stable_gate_and_legacy_strict_alias_stay_synchronized(
    coarse_ng: int,
    unstable: bool,
    ari_min: float | None,
    expected: bool,
) -> None:
    output: dict[str, object] = {}
    observed = MODULE.apply_m1_screen_flags(
        output,
        {
            "coarse_ng": coarse_ng,
            "unstable": unstable,
            "modal_assignment_ari_min": ari_min,
        },
    )

    assert observed is expected
    assert output["stable_null_multigroup"] is expected
    assert output["strict_confirm_candidate"] is expected
    assert output["strict_confirm_candidate_is_formal_r1_claim"] is False


def test_stable_detail_has_identity_latest_tags_labels_and_artifact_hashes(tmp_path: Path) -> None:
    region_dir = tmp_path / "chr1" / "100"
    artifact_contents = {
        region_dir / "reads" / "reads.tsv": b"reads\n",
        region_dir / "distance" / "BERNOULLI" / "matrix.csv": b"distance\n",
        region_dir / "methylation" / "methylation.csv": b"methylation\n",
    }
    for path, content in artifact_contents.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    reads = {
        "0": {
            "read_name": "read-a",
            "chr": "chr1",
            "start": "90",
            "end": "150",
            "mapq": "60",
            "strand": "+",
            "hp": "1-2",
            "ps": 101,
            "latest_tag_full_identity_count": 1,
        }
    }
    detail = MODULE.stable_cluster_detail(
        {"sample": "S1", "chrom": "chr1", "pos": 100, "region_dir": str(region_dir)},
        reads,
        ["0"],
        ["group_1"],
        ["0"],
        ["group_1"],
        {},
        {
            "coarse_ng": 2,
            "unstable": False,
            "modal_assignment_ari_min": 0.8,
            "coarse_split_trace": [],
            "fine_split_trace": [],
        },
    )

    assert detail["strict_confirm_status"] == "NOT_RUN"
    assert detail["strict_confirm_candidate"] is True
    assert detail["strict_confirm_candidate_is_formal_r1_claim"] is False
    assert detail["strict_methyl_partition_robustness_status"] == "NOT_EVALUATED_AT_M1_SCREEN"
    assert detail["artifact_identity_contract"] == "sha256_size_path_v1"
    assert set(detail["primary_artifacts"]) == {"reads", "distance_matrix", "methylation_matrix"}
    for name, record in detail["primary_artifacts"].items():
        path = Path(record["path"])
        assert path.is_file(), name
        assert record["size_bytes"] == path.stat().st_size
        assert record["sha256"] == MODULE.sha256(path)
    assert detail["core_reads"] == [
        {
            "read_id": "0",
            "read_name": "read-a",
            "chrom": "chr1",
            "start": 90,
            "end": 150,
            "mapq": 60,
            "strand": "+",
            "latest_hp": "1-2",
            "latest_ps": 101,
            "latest_tag_full_identity_count": 1,
            "label": "group_1",
        }
    ]


def test_stable_detail_rejects_assignment_below_ari_gate(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError, match="does not pass the M1 gate"):
        MODULE.stable_cluster_detail(
            {"sample": "S1", "chrom": "chr1", "pos": 100, "region_dir": str(tmp_path)},
            {},
            [],
            [],
            [],
            [],
            {},
            {
                "coarse_ng": 2,
                "unstable": False,
                "modal_assignment_ari_min": 0.79,
            },
        )


def test_completed_receipt_allows_summary_omission_but_requires_site_artifacts(tmp_path: Path) -> None:
    sample_root = tmp_path / "S1"
    sample_root.mkdir()
    sidecar = tmp_path / "latest.tags.tsv.gz"
    sidecar.write_bytes(b"sidecar")
    entry = {
        "sample": "S1",
        "counts": {"all_ssnv": 1000},
        "all_ssnv_vcf": {"sha256": "frozen-vcf"},
        "latest_read_tag_sidecar": {"path": str(sidecar)},
    }
    receipt = {
        "schema_name": "intersubmod.all_ssnv_site_run",
        "sample": "S1",
        "pass": True,
        "exit_code": 0,
        "input_vcf_sha256": "frozen-vcf",
        "output_dir": str(sample_root),
        "latest_read_tag_sidecar": {"path": str(sidecar)},
        "validation": {
            "expected_vcf_sites": 1000,
            "summary_rows": 881,
            "summary_rows_not_emitted": 119,
            "reads_files": 1000,
            "bernoulli_matrix_files": 1000,
            "methylation_files": 1000,
            "pass": True,
        },
    }
    receipt_path = sample_root / "run_receipt.json"
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")

    validated = MODULE.validate_completed_sample_run(sample_root, entry)
    assert validated["validation"]["summary_rows_not_emitted"] == 119
    receipt["validation"]["reads_files"] = 999
    receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
    with pytest.raises(RuntimeError, match="reads_files"):
        MODULE.validate_completed_sample_run(sample_root, entry)


def test_chunking_submits_groups_not_one_future_per_site() -> None:
    tasks = [{"index": index} for index in range(10)]
    chunks = list(MODULE.chunked(tasks, 3))
    assert [len(chunk) for chunk in chunks] == [3, 3, 3, 1]
    assert len(chunks) < len(tasks)


def test_seed_parallel_phylo_is_exactly_equal_to_serial() -> None:
    methylation = np.asarray(
        [
            [0.05, 0.10, 0.15, 0.10, 0.05, 0.10],
            [0.10, 0.05, 0.10, 0.15, 0.10, 0.05],
            [0.15, 0.10, 0.05, 0.10, 0.15, 0.10],
            [0.08, 0.12, 0.10, 0.08, 0.12, 0.10],
            [0.90, 0.85, 0.95, 0.90, 0.85, 0.95],
            [0.85, 0.95, 0.90, 0.85, 0.95, 0.90],
            [0.95, 0.90, 0.85, 0.95, 0.90, 0.85],
            [0.88, 0.92, 0.90, 0.88, 0.92, 0.90],
        ],
        dtype=float,
    )
    distance = MODULE.F.bernoulli_distance(methylation)
    serial = MODULE.F.analyze_phylo(distance, methylation, seeds=3, rnull=3)
    parallel = MODULE.analyze_phylo_parallel(
        distance,
        methylation,
        workers=4,
        seeds=3,
        rnull=3,
    )

    assert parallel == serial


def test_phylo_parallel_configuration_is_fail_closed() -> None:
    with pytest.raises(ValueError, match="at least the clustering minimum"):
        MODULE.configure_phylo_parallel(2, 5)
    with pytest.raises(ValueError, match="must be positive"):
        MODULE.configure_phylo_parallel(0, 0)
    MODULE.configure_phylo_parallel(1, 0)


def test_existing_output_directory_is_never_overwritten(tmp_path: Path) -> None:
    output = tmp_path / "existing"
    output.mkdir()
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        MODULE.create_output_dir(output)
