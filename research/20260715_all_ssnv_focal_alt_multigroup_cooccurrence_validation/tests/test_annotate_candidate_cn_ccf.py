from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "annotate_candidate_cn_ccf.py"
)
SPEC = importlib.util.spec_from_file_location("annotate_candidate_cn_ccf", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


GOLDEN = MODULE.VariantKey("chr1", 47163, "T", "C")
H1437_NEAR = MODULE.VariantKey("chr1", 700000, "A", "G")
HCC_PAIR = MODULE.VariantKey("chr1", 877772, "G", "C")
H2009_SITE = MODULE.VariantKey("chr1", 100, "A", "C")
HCC1954_SITE = MODULE.VariantKey("chr1", 200, "C", "T")


def test_current_release_cn_ccf_paths_bind_v7_and_v3() -> None:
    assert MODULE.CANONICAL_INPUT.parent.name == (
        "methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_"
        "source_locked_command_parity"
    )
    assert (
        MODULE.CANONICAL_OUTPUT_DIR.name
        == "candidate_cn_ccf_annotations_v3_m2v5_source_authority_v5"
    )


@dataclass(frozen=True)
class FixtureAuthority:
    paths: MODULE.AuthorityPaths
    config: Path
    h1437_cn: Path


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        handle_context = gzip.open(path, "wt", encoding="utf-8", newline="")
    else:
        handle_context = path.open("w", encoding="utf-8", newline="")
    with handle_context as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def file_receipt(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": MODULE.sha256_file(path),
    }


def cn_row(
    start: int,
    end: int,
    segment_id: str,
    total: float,
    minor: float,
) -> dict[str, object]:
    return {
        "chromosome": "chr1",
        "start": start,
        "end": end,
        "segment_id": segment_id,
        "copyNumber": total,
        "minorAlleleCopyNumber": minor,
    }


CN_FIELDS = [
    "chromosome",
    "start",
    "end",
    "segment_id",
    "copyNumber",
    "minorAlleleCopyNumber",
]


def segment_for_variant(sample: str, variant: MODULE.VariantKey) -> MODULE.Segment:
    if sample in ("HCC1395", "HCC1395_DORADO"):
        return MODULE.Segment("chr1", 520001, 930000, "chr1_seg3", 2.2157, 0.1876, 2, 2, 0, True)
    if sample == "H1437" and variant == GOLDEN:
        return MODULE.Segment("chr1", 10001, 660000, "chr1_seg1", 3.4338, 0.793, 3, 2, 1, False)
    if sample == "H1437":
        return MODULE.Segment("chr1", 660001, 1000000, "chr1_seg2", 2.1, 0.9, 2, 1, 1, True)
    if sample == "H2009":
        return MODULE.Segment("chr1", 1, 1000000, "chr1_seg1", 2.0, 1.0, 2, 1, 1, True)
    return MODULE.Segment("chr1", 1, 1000000, "chr1_seg1", 3.0, 1.0, 3, 2, 1, True)


def variants_for_role(role: MODULE.BundleRole) -> list[MODULE.VariantKey]:
    if "hcc1395" in role.bundle_id:
        return [HCC_PAIR]
    sample = role.samples[0]
    if sample == "H1437":
        return [H1437_NEAR] if role.near_integer_only else [GOLDEN, H1437_NEAR]
    if sample == "H2009":
        return [H2009_SITE]
    if sample == "HCC1954":
        return [HCC1954_SITE]
    raise AssertionError(role)


def result_values(sample: str, variant: MODULE.VariantKey) -> tuple[str, str, str, str]:
    if sample == "HCC1395":
        return "2", "1.0000", "0.0000", "0.7499"
    if sample == "HCC1395_DORADO":
        return "0", "1.0000", "0.0000", "0.8599"
    if sample == "H1437" and variant == GOLDEN:
        return "0", "1.0000", "0.0000", "0.9659"
    return "1", "0.7000", "0.0500", "0.9000"


def make_bundle(
    root: Path,
    role: MODULE.BundleRole,
) -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
    input_root = root / "data" / "pyclone_inputs"
    run_root = root / "runs" / "pyclone_vi" / role.bundle_id
    input_path = input_root / f"{role.bundle_id}.pyclone_input.tsv.gz"
    metadata_path = input_root / f"{role.bundle_id}.site_metadata.tsv.gz"
    qa_path = input_root / f"{role.bundle_id}.qa.json"
    results_path = run_root / "results.tsv.gz"
    status_path = run_root / "status.json"
    variants = variants_for_role(role)

    input_rows: list[dict[str, object]] = []
    metadata_rows: list[dict[str, object]] = []
    result_rows: list[dict[str, object]] = []
    for variant in variants:
        segment = segment_for_variant(role.samples[0], variant)
        metadata_rows.append(
            {
                "mutation_id": variant.mutation_id,
                "chrom": variant.chrom,
                "pos": variant.pos,
                "ref": variant.ref,
                "alt": variant.alt,
                "segment_id": segment.segment_id,
                "segment_start": segment.start,
                "segment_end": segment.end,
                "total_cn_raw": segment.total_raw,
                "minor_cn_raw": segment.minor_raw,
                "total_cn_discrete": segment.total_discrete,
                "major_cn": segment.major_discrete,
                "minor_cn": segment.minor_discrete,
                "near_integer": int(segment.near_integer),
            }
        )
        for sample in role.samples:
            input_rows.append(
                {
                    "mutation_id": variant.mutation_id,
                    "sample_id": sample,
                    "ref_counts": 10,
                    "alt_counts": 10,
                    "normal_cn": 2,
                    "major_cn": segment.major_discrete,
                    "minor_cn": segment.minor_discrete,
                    "tumour_content": 0.95,
                }
            )
            cluster, prevalence, std, probability = result_values(sample, variant)
            result_rows.append(
                {
                    "mutation_id": variant.mutation_id,
                    "sample_id": sample,
                    "cluster_id": cluster,
                    "cellular_prevalence": prevalence,
                    "cellular_prevalence_std": std,
                    "cluster_assignment_prob": probability,
                }
            )

    write_tsv(
        input_path,
        [
            "mutation_id",
            "sample_id",
            "ref_counts",
            "alt_counts",
            "normal_cn",
            "major_cn",
            "minor_cn",
            "tumour_content",
        ],
        input_rows,
    )
    write_tsv(
        metadata_path,
        [
            "mutation_id",
            "chrom",
            "pos",
            "ref",
            "alt",
            "segment_id",
            "segment_start",
            "segment_end",
            "total_cn_raw",
            "minor_cn_raw",
            "total_cn_discrete",
            "major_cn",
            "minor_cn",
            "near_integer",
        ],
        metadata_rows,
    )
    write_tsv(
        results_path,
        [
            "mutation_id",
            "sample_id",
            "cluster_id",
            "cellular_prevalence",
            "cellular_prevalence_std",
            "cluster_assignment_prob",
        ],
        result_rows,
    )
    counters = {
        "universe_exact_keys": len(variants),
        "excluded_missing_counts_any_sample": 0,
        "excluded_no_measured_cn_segment": 0,
        "excluded_not_near_integer_cn": 0,
        "included_mutations": len(variants),
        "included_rows": len(result_rows),
        "count_conservation_fail": 0,
        "duplicate_mutation_sample_rows": 0,
    }
    artifacts = {
        "pyclone_input": str(input_path.resolve()),
        "pyclone_input_sha256": MODULE.sha256_file(input_path),
        "site_metadata": str(metadata_path.resolve()),
        "site_metadata_sha256": MODULE.sha256_file(metadata_path),
    }
    qa = {
        "schema_name": "intersubmod.pyclone_input_qa",
        "schema_version": "1.0.0",
        "bundle_id": role.bundle_id,
        "status": "PASS",
        "failures": [],
        "samples": list(role.samples),
        "near_integer_only": role.near_integer_only,
        "counters": counters,
        "artifacts": artifacts,
    }
    qa_path.parent.mkdir(parents=True, exist_ok=True)
    qa_path.write_text(json.dumps(qa, sort_keys=True), encoding="utf-8")
    provenance_artifacts = dict(artifacts, qa_json=str(qa_path.resolve()))
    bundle = {
        **qa,
        "artifacts": provenance_artifacts,
    }
    status = {
        "schema_name": "intersubmod.pyclone_vi_fit_receipt",
        "schema_version": "1.0.0",
        "bundle_id": role.bundle_id,
        "status": "PASS",
        "fit_exit_code": 0,
        "write_results_exit_code": 0,
        "input": {
            "path": str(input_path.resolve()),
            "sha256": MODULE.sha256_file(input_path),
            "qa_path": str(qa_path.resolve()),
            "qa_sha256": MODULE.sha256_file(qa_path),
        },
        "outputs": {
            "results": str(results_path.resolve()),
            "results_sha256": MODULE.sha256_file(results_path),
        },
        "result_shape": {
            "rows": len(result_rows),
            "mutations": len(variants),
            "samples": list(role.samples),
        },
    }
    status_path.parent.mkdir(parents=True, exist_ok=True)
    status_path.write_text(json.dumps(status, sort_keys=True), encoding="utf-8")
    summary_receipt = {
        "results_path": str(results_path.resolve()),
        "results_sha256": MODULE.sha256_file(results_path),
        "status_path": str(status_path.resolve()),
    }
    fit_summary = {
        "bundle_id": role.bundle_id,
        "rows": len(result_rows),
        "mutations": len(variants),
        "samples": list(role.samples),
    }
    return bundle, summary_receipt, fit_summary


def make_authority(tmp_path: Path) -> FixtureAuthority:
    root = tmp_path / "authority"
    config_path = root / "config" / "pyclone_validation_config.json"
    provenance_path = root / "data" / "pyclone_inputs" / "provenance.json"
    summary_path = root / "runs" / "pyclone_vi" / "analysis" / "analysis_summary.json"
    manifest_path = root / "manifest.json"
    hcc_cn = root / "cn" / "HCC1395.tsv"
    h1437_cn = root / "cn" / "H1437.tsv"
    h2009_cn = root / "cn" / "H2009.tsv"
    hcc1954_cn = root / "cn" / "HCC1954.tsv"
    write_tsv(hcc_cn, CN_FIELDS, [cn_row(520001, 930000, "chr1_seg3", 2.2157, 0.1876)])
    write_tsv(
        h1437_cn,
        CN_FIELDS,
        [
            cn_row(10001, 660000, "chr1_seg1", 3.4338, 0.793),
            cn_row(660001, 1000000, "chr1_seg2", 2.1, 0.9),
        ],
    )
    write_tsv(h2009_cn, CN_FIELDS, [cn_row(1, 1000000, "chr1_seg1", 2.0, 1.0)])
    write_tsv(hcc1954_cn, CN_FIELDS, [cn_row(1, 1000000, "chr1_seg1", 3.0, 1.0)])
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text("{}\n", encoding="ascii")
    config = {
        "schema_name": "intersubmod.pyclone_vi_validation_config",
        "schema_version": "1.0.0",
        "analysis_scope": "GRCh38 autosomes chr1-chr22, biallelic SNVs, exact joins",
        "source_manifest": str(manifest_path.resolve()),
        "hcc1395_pair": {
            "samples": ["HCC1395", "HCC1395_DORADO"],
            "copy_number_tsv": str(hcc_cn.resolve()),
            "purity": 0.96,
            "dorado_copy_number_policy": "shared-CN sensitivity",
        },
        "individual_samples": {
            "H1437": {"copy_number_tsv": str(h1437_cn.resolve()), "purity": 0.95},
            "H2009": {"copy_number_tsv": str(h2009_cn.resolve()), "purity": 0.95},
            "HCC1954": {"copy_number_tsv": str(hcc1954_cn.resolve()), "purity": 0.66},
        },
        "blocked_samples": {
            "COLO829": "BLOCKED_CN_MISFIT",
            "HCC1937": "BLOCKED_CN_MISFIT",
        },
    }
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(config, sort_keys=True), encoding="utf-8")
    roles = MODULE.expected_bundle_roles(config)
    bundles: list[dict[str, object]] = []
    summary_receipts: dict[str, object] = {}
    fit_summaries: dict[str, object] = {}
    for bundle_id in sorted(roles):
        bundle, receipt, fit_summary = make_bundle(root, roles[bundle_id])
        bundles.append(bundle)
        summary_receipts[bundle_id] = receipt
        fit_summaries[bundle_id] = fit_summary

    sources = {
        "config": file_receipt(config_path),
        "source_manifest": file_receipt(manifest_path),
        str(hcc_cn.resolve()): file_receipt(hcc_cn),
        str(h1437_cn.resolve()): file_receipt(h1437_cn),
        str(h2009_cn.resolve()): file_receipt(h2009_cn),
        str(hcc1954_cn.resolve()): file_receipt(hcc1954_cn),
    }
    provenance = {
        "schema_name": "intersubmod.pyclone_input_build_provenance",
        "schema_version": "1.0.0",
        "config": config,
        "source_receipts": sources,
        "bundles": bundles,
        "blocked_samples": config["blocked_samples"],
    }
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.write_text(json.dumps(provenance, sort_keys=True), encoding="utf-8")
    summary = {
        "schema_name": "intersubmod.pyclone_vi_conditional_clustering_analysis",
        "schema_version": "1.0.0",
        "all_ready": True,
        "expected_bundle_count": len(roles),
        "pass_bundle_count": len(roles),
        "pending_or_failed_bundles": [],
        "source_receipts": summary_receipts,
        "fit_summaries": fit_summaries,
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, sort_keys=True), encoding="utf-8")
    return FixtureAuthority(
        MODULE.AuthorityPaths(config_path, provenance_path, summary_path),
        config_path,
        h1437_cn,
    )


def write_candidates(path: Path, rows: list[dict[str, object]]) -> None:
    fields = ["sample", "chrom", "pos", "ref", "alt"]
    if any("selected" in row for row in rows):
        fields.append("selected")
    write_tsv(path, fields, rows)


def candidate_row(sample: str, variant: MODULE.VariantKey, **extra: object) -> dict[str, object]:
    return {
        "sample": sample,
        "chrom": variant.chrom,
        "pos": variant.pos,
        "ref": variant.ref,
        "alt": variant.alt,
        **extra,
    }


def read_annotations(output_dir: Path) -> list[dict[str, str]]:
    with gzip.open(output_dir / "candidate_cn_ccf_annotations.tsv.gz", "rt", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_exact_allele_join_rejects_same_position_other_alt(tmp_path: Path) -> None:
    authority = make_authority(tmp_path)
    candidates = tmp_path / "candidates.tsv.gz"
    mismatch = MODULE.VariantKey("chr1", GOLDEN.pos, GOLDEN.ref, "A")
    write_candidates(candidates, [candidate_row("H1437", mismatch)])
    output_dir = tmp_path / "output"

    MODULE.annotate_candidates(candidates, output_dir, authority_paths=authority.paths)
    row = read_annotations(output_dir)[0]
    assert row["cn_status"] == "AVAILABLE_EXACT_SEGMENT"
    assert row["pyclone_status"] == "NOT_IN_FIT_UNIVERSE"
    assert row["pyclone_fit_local_cluster_id"] == ""


def test_segment_index_uses_inclusive_boundaries_gap_and_rejects_overlap(tmp_path: Path) -> None:
    cn_path = tmp_path / "segments.tsv"
    write_tsv(
        cn_path,
        CN_FIELDS,
        [cn_row(10, 20, "left", 2.0, 1.0), cn_row(30, 40, "right", 3.0, 1.0)],
    )
    index = MODULE.SegmentIndex(cn_path)
    assert index.find(MODULE.VariantKey("chr1", 10, "A", "C")).segment_id == "left"
    assert index.find(MODULE.VariantKey("chr1", 20, "A", "C")).segment_id == "left"
    assert index.find(MODULE.VariantKey("chr1", 21, "A", "C")) is None
    assert index.find(MODULE.VariantKey("chr1", 30, "A", "C")).segment_id == "right"
    assert index.find(MODULE.VariantKey("chr1", 40, "A", "C")).segment_id == "right"

    overlap_path = tmp_path / "overlap.tsv"
    write_tsv(
        overlap_path,
        CN_FIELDS,
        [cn_row(10, 20, "first", 2.0, 1.0), cn_row(20, 30, "second", 2.0, 1.0)],
    )
    with pytest.raises(MODULE.ContractError, match="Overlapping CN segments"):
        MODULE.SegmentIndex(overlap_path)


def test_h1437_real_golden_exact_metadata_and_fit_result() -> None:
    root = Path(__file__).resolve().parents[2] / "20260712_vaf_ccf_subclone_crosssoftware_validation"
    cn_path = Path(
        "/big7_disk/liaoyoyo2001/cnv_sv_work/H1437/savana_wgs/cna_normalhet/"
        "H1437_segmented_absolute_copy_number.tsv"
    )
    metadata_path = root / "data/pyclone_inputs/H1437_individual_main.site_metadata.tsv.gz"
    results_path = root / "runs/pyclone_vi/H1437_individual_main/results.tsv.gz"
    qa_path = root / "data/pyclone_inputs/H1437_individual_main.qa.json"
    if not all(path.is_file() for path in (cn_path, metadata_path, results_path, qa_path)):
        pytest.skip("Canonical 20260712 H1437 authority is unavailable")
    segment = MODULE.SegmentIndex(cn_path).find(GOLDEN)
    assert segment is not None
    assert (segment.segment_id, segment.total_raw, segment.minor_raw) == (
        "chr1_seg1",
        pytest.approx(3.4338),
        pytest.approx(0.793),
    )
    qa = json.loads(qa_path.read_text())
    role = MODULE.BundleRole("H1437_individual_main", ("H1437",), False, "INDIVIDUAL_MAIN_PRIMARY")
    key = MODULE.CandidateKey("H1437", GOLDEN)
    hits = MODULE.scan_bundle(
        role,
        metadata_path,
        results_path,
        frozenset((key,)),
        qa["counters"]["included_mutations"],
        qa["counters"]["included_rows"],
    )
    assert hits[key].metadata["segment_id"] == "chr1_seg1"
    assert hits[key].result["cellular_prevalence"] == "1.0000"
    assert hits[key].result["cluster_assignment_prob"] == "0.9659"


def test_hcc_pair_shared_cn_and_fit_local_clusters_are_not_cross_fit_ids(tmp_path: Path) -> None:
    authority = make_authority(tmp_path)
    candidates = tmp_path / "candidates.tsv"
    write_candidates(
        candidates,
        [candidate_row("HCC1395", HCC_PAIR), candidate_row("HCC1395_DORADO", HCC_PAIR)],
    )
    output_dir = tmp_path / "output"
    MODULE.annotate_candidates(candidates, output_dir, authority_paths=authority.paths)
    rows = {row["sample"]: row for row in read_annotations(output_dir)}

    assert rows["HCC1395"]["cn_segment_id"] == rows["HCC1395_DORADO"]["cn_segment_id"] == "chr1_seg3"
    assert rows["HCC1395"]["cn_source_sha256"] == rows["HCC1395_DORADO"]["cn_source_sha256"]
    assert rows["HCC1395"]["cn_status"] == "AVAILABLE_EXACT_SEGMENT"
    assert rows["HCC1395_DORADO"]["cn_status"] == "SHARED_CN_SENSITIVITY"
    assert rows["HCC1395_DORADO"]["independent_cn"] == "false"
    assert rows["HCC1395"]["pyclone_primary_bundle_id"] != rows["HCC1395_DORADO"]["pyclone_primary_bundle_id"]
    assert rows["HCC1395"]["pyclone_fit_local_cluster_id"] == "2"
    assert rows["HCC1395_DORADO"]["pyclone_fit_local_cluster_id"] == "0"
    assert rows["HCC1395"]["pyclone_primary_model_role"] == "SEPARATE_MAIN_PRIMARY"


def test_blocked_samples_never_receive_cn2_or_pyclone_values(tmp_path: Path) -> None:
    authority = make_authority(tmp_path)
    candidates = tmp_path / "blocked.tsv"
    unknown = MODULE.VariantKey("chr1", 500, "A", "T")
    write_candidates(
        candidates,
        [candidate_row("COLO829", unknown), candidate_row("HCC1937", unknown)],
    )
    output_dir = tmp_path / "output"
    MODULE.annotate_candidates(candidates, output_dir, authority_paths=authority.paths)
    rows = read_annotations(output_dir)
    assert len(rows) == 2
    for row in rows:
        assert row["cn_status"] == "BLOCKED_CN_MISFIT"
        assert row["pyclone_status"] == "BLOCKED_CN_MISFIT"
        assert row["savana_total_cn_discrete"] == ""
        assert row["purity_model_value"] == ""
        assert row["pyclone_fit_local_cluster_id"] == ""


def test_not_in_fit_universe_preserves_one_output_row(tmp_path: Path) -> None:
    authority = make_authority(tmp_path)
    candidates = tmp_path / "candidate.tsv"
    outside_fit = MODULE.VariantKey("chr1", 300, "G", "T")
    write_candidates(candidates, [candidate_row("H2009", outside_fit)])
    output_dir = tmp_path / "output"
    receipt = MODULE.annotate_candidates(candidates, output_dir, authority_paths=authority.paths)
    rows = read_annotations(output_dir)
    assert len(rows) == receipt["conservation"]["rows_in"] == receipt["conservation"]["rows_out"] == 1
    assert receipt["pass"] is True
    assert receipt["status"] == "PASS"
    assert receipt["execution_status"] == "EXECUTION_PASS"
    assert receipt["reason"] is None
    assert receipt["n_selected_candidates"] == 1
    assert receipt["scientific_interpretation"]["c1_formed"] is False
    assert rows[0]["cn_status"] == "AVAILABLE_EXACT_SEGMENT"
    assert rows[0]["pyclone_status"] == "NOT_IN_FIT_UNIVERSE"
    assert rows[0]["claim_ceiling"] == MODULE.CLAIM_CEILING


def test_source_hash_drift_fails_before_formal_output(tmp_path: Path) -> None:
    authority = make_authority(tmp_path)
    with authority.h1437_cn.open("a", encoding="utf-8") as handle:
        handle.write("\n")
    candidates = tmp_path / "candidate.tsv"
    write_candidates(candidates, [candidate_row("H1437", GOLDEN)])
    output_dir = tmp_path / "must_not_exist"
    with pytest.raises(MODULE.ContractError, match="SHA drift"):
        MODULE.annotate_candidates(candidates, output_dir, authority_paths=authority.paths)
    assert not output_dir.exists()


def test_fit_status_and_provenance_count_drift_fail_closed(tmp_path: Path) -> None:
    status_authority = make_authority(tmp_path / "status_case")
    summary = json.loads(status_authority.paths.analysis_summary.read_text())
    status_path = Path(summary["source_receipts"]["H1437_individual_main"]["status_path"])
    status = json.loads(status_path.read_text())
    status["status"] = "FAIL"
    status_path.write_text(json.dumps(status, sort_keys=True), encoding="utf-8")
    candidates = tmp_path / "status.tsv"
    write_candidates(candidates, [candidate_row("H1437", GOLDEN)])
    with pytest.raises(MODULE.ContractError, match="not terminal PASS"):
        MODULE.annotate_candidates(
            candidates,
            tmp_path / "status_output",
            authority_paths=status_authority.paths,
        )
    assert not (tmp_path / "status_output").exists()

    count_authority = make_authority(tmp_path / "count_case")
    provenance = json.loads(count_authority.paths.input_provenance.read_text())
    target = next(
        bundle
        for bundle in provenance["bundles"]
        if bundle["bundle_id"] == "H1437_individual_main"
    )
    target["counters"]["included_mutations"] += 1
    count_authority.paths.input_provenance.write_text(
        json.dumps(provenance, sort_keys=True), encoding="utf-8"
    )
    with pytest.raises(MODULE.ContractError, match="QA counter drift"):
        MODULE.annotate_candidates(
            candidates,
            tmp_path / "count_output",
            authority_paths=count_authority.paths,
        )
    assert not (tmp_path / "count_output").exists()


def test_result_key_must_first_match_exact_site_metadata(tmp_path: Path) -> None:
    metadata_path = tmp_path / "metadata.tsv.gz"
    results_path = tmp_path / "results.tsv.gz"
    segment = segment_for_variant("H1437", GOLDEN)
    write_tsv(
        metadata_path,
        [
            "mutation_id",
            "chrom",
            "pos",
            "ref",
            "alt",
            "segment_id",
            "segment_start",
            "segment_end",
            "total_cn_raw",
            "minor_cn_raw",
            "total_cn_discrete",
            "major_cn",
            "minor_cn",
            "near_integer",
        ],
        [
            {
                "mutation_id": GOLDEN.mutation_id,
                "chrom": GOLDEN.chrom,
                "pos": GOLDEN.pos,
                "ref": GOLDEN.ref,
                "alt": GOLDEN.alt,
                "segment_id": segment.segment_id,
                "segment_start": segment.start,
                "segment_end": segment.end,
                "total_cn_raw": segment.total_raw,
                "minor_cn_raw": segment.minor_raw,
                "total_cn_discrete": segment.total_discrete,
                "major_cn": segment.major_discrete,
                "minor_cn": segment.minor_discrete,
                "near_integer": 0,
            }
        ],
    )
    mismatched = MODULE.VariantKey(GOLDEN.chrom, GOLDEN.pos, GOLDEN.ref, "A")
    write_tsv(
        results_path,
        [
            "mutation_id",
            "sample_id",
            "cluster_id",
            "cellular_prevalence",
            "cellular_prevalence_std",
            "cluster_assignment_prob",
        ],
        [
            {
                "mutation_id": mismatched.mutation_id,
                "sample_id": "H1437",
                "cluster_id": 0,
                "cellular_prevalence": 1.0,
                "cellular_prevalence_std": 0.0,
                "cluster_assignment_prob": 0.9,
            }
        ],
    )
    role = MODULE.BundleRole(
        "H1437_individual_main",
        ("H1437",),
        False,
        "INDIVIDUAL_MAIN_PRIMARY",
    )
    with pytest.raises(MODULE.ContractError, match="lacks exact metadata allele"):
        MODULE.scan_bundle(
            role,
            metadata_path,
            results_path,
            frozenset(),
            expected_mutations=1,
            expected_rows=1,
        )


def test_output_schema_has_no_prohibited_claim_fields() -> None:
    lowered = [column.lower() for column in MODULE.OUTPUT_COLUMNS]
    forbidden = (
        "true_cn",
        "validated_ccf",
        "global_clone",
        "subclone_confirmed",
        "orthogonal",
        "independent_confirmation",
    )
    for token in forbidden:
        assert all(token not in column for column in lowered)
    assert "pyclone_fit_local_cluster_id" in MODULE.OUTPUT_COLUMNS
    assert MODULE.CLAIM_CEILING == "conditional_annotation_only_not_orthogonal_confirmation"


def test_zero_selected_candidates_writes_fixed_schema_and_passing_receipt(tmp_path: Path) -> None:
    authority = make_authority(tmp_path)
    candidates = tmp_path / "candidate.tsv.gz"
    write_candidates(candidates, [candidate_row("H1437", GOLDEN, selected="false")])
    output_dir = tmp_path / "output"
    receipt = MODULE.annotate_candidates(
        candidates,
        output_dir,
        selection_column="selected",
        selection_value="true",
        authority_paths=authority.paths,
    )
    assert read_annotations(output_dir) == []
    assert receipt["pass"] is True
    assert receipt["status"] == "NOT_APPLICABLE"
    assert receipt["execution_status"] == "NOT_APPLICABLE"
    assert receipt["reason"] == "ZERO_SELECTED_CANDIDATES"
    assert receipt["n_selected_candidates"] == 0
    assert receipt["input"]["rows_read_total"] == 1
    assert receipt["conservation"] == {
        "rows_in": 0,
        "rows_out": 0,
        "rows_in_equals_rows_out": True,
    }
    on_disk = json.loads((output_dir / "receipt.json").read_text())
    with gzip.open(output_dir / "candidate_cn_ccf_annotations.tsv.gz", "rt", encoding="utf-8") as handle:
        assert len(handle.readlines()) == 1
    assert on_disk["output"]["columns"] == list(MODULE.OUTPUT_COLUMNS)
    assert on_disk["output"]["sha256"] == MODULE.sha256_file(
        output_dir / "candidate_cn_ccf_annotations.tsv.gz"
    )
    assert on_disk["authority"]["all_source_hashes"]
    assert on_disk["scientific_interpretation"] == {
        "is_negative_result": False,
        "c1_formed": False,
        "statement": (
            "ZERO_SELECTED_CANDIDATES is not a biological negative; "
            "C1 was not evaluated or formed."
        ),
    }
    assert on_disk["pass_semantics"] == MODULE.PASS_SEMANTICS
    assert on_disk["source_authority"] == {
        "authority_id": "TEST_ONLY_UNSIGNED_AUTHORITY",
        "pass": True,
    }
    assert on_disk["code"] == MODULE.capture_source_identity()
    assert on_disk["source_lock"]["all_sources_read_only_and_unchanged"] is True
    assert oct((output_dir / "receipt.json").stat().st_mode & 0o777) == "0o444"
    assert oct(
        (output_dir / "candidate_cn_ccf_annotations.tsv.gz").stat().st_mode & 0o777
    ) == "0o444"


def test_formal_cn_producer_rejects_injected_argv() -> None:
    authority = {"authority_id": MODULE.SOURCE_AUTHORITY.AUTHORITY_ID, "pass": True}
    with pytest.raises(MODULE.ContractError, match="direct-CLI canonical-process only"):
        MODULE.resolve_release_command(["--injected"], authority)


def test_formal_cn_producer_rejects_process_missing_isolation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    authority = {"authority_id": MODULE.SOURCE_AUTHORITY.AUTHORITY_ID, "pass": True}
    expected = MODULE.canonical_task_b_command()
    assert expected[:4] == MODULE.canonical_python_prefix()
    monkeypatch.setattr(sys, "argv", [expected[4], *expected[5:]])
    monkeypatch.setattr(
        MODULE,
        "observed_process_command",
        lambda: [expected[0], expected[4], *expected[5:]],
    )
    with pytest.raises(MODULE.ContractError, match="direct-CLI canonical-process only"):
        MODULE.resolve_release_command(None, authority)


def test_duplicate_candidate_key_and_existing_output_fail_closed(tmp_path: Path) -> None:
    authority = make_authority(tmp_path)
    candidates = tmp_path / "duplicate.tsv"
    row = candidate_row("H1437", GOLDEN)
    write_candidates(candidates, [row, row])
    with pytest.raises(MODULE.ContractError, match="Duplicate candidate key"):
        MODULE.annotate_candidates(candidates, tmp_path / "output", authority_paths=authority.paths)

    unique = tmp_path / "unique.tsv"
    write_candidates(unique, [row])
    existing = tmp_path / "existing"
    existing.mkdir()
    with pytest.raises(MODULE.ContractError, match="must not exist"):
        MODULE.annotate_candidates(unique, existing, authority_paths=authority.paths)
