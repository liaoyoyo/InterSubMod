from __future__ import annotations

import copy
import importlib.util
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = (
    REPO_ROOT
    / "research/20260723_production_exact_ps_strict_read_linkage/scripts"
)
HCC_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260723_production_exact_ps_strict_read_linkage/hcc1395_strict_regions_v2"
)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


DATA_BUILDER = load_module(
    "all7_region_data_builder",
    SCRIPT_ROOT / "build_all7_region_report_data.py",
)
ARTIFACT_BUILDER = load_module(
    "all7_region_artifact_builder",
    SCRIPT_ROOT / "build_all7_region_artifact.py",
)


def write_data_package(tmp_path: Path, data: dict) -> Path:
    data_dir = tmp_path / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    audit_dir = tmp_path / "topology_status_audit"
    audit_dir.mkdir(parents=True, exist_ok=True)
    audit_path = audit_dir / "all7_topology_completion_audit.json"
    DATA_BUILDER.write_json(
        audit_path,
        {
            "schema_name": "intersubmod.all7_topology_completion_audit",
            "schema_version": "1.0.0",
            "all_pass": True,
        },
    )
    audit_receipt_path = audit_dir / "receipt.json"
    DATA_BUILDER.write_json(
        audit_receipt_path,
        {
            "schema_name": "intersubmod.all7_topology_completion_audit.receipt",
            "schema_version": "1.0.0",
            "all_pass": True,
        },
    )
    audit_sidecar_path = audit_dir / "receipt.json.sha256"
    audit_sidecar_path.write_text(
        f"{ARTIFACT_BUILDER.sha256_path(audit_receipt_path)}  receipt.json\n",
        encoding="ascii",
    )
    data["supporting_inputs"] = {
        "topology_completion_audit": ARTIFACT_BUILDER.file_identity(audit_path),
        "topology_completion_audit_receipt": ARTIFACT_BUILDER.file_identity(
            audit_receipt_path
        ),
        "topology_completion_audit_receipt_sidecar": (
            ARTIFACT_BUILDER.file_identity(audit_sidecar_path)
        ),
    }
    data_path = data_dir / "all7_report_data.json"
    DATA_BUILDER.write_json(data_path, data)
    receipt = {
        "schema_name": "intersubmod.strict_region_all7_report_data_receipt",
        "schema_version": DATA_BUILDER.SCHEMA_VERSION,
        "all_pass": True,
        "checks": data["checks"],
        "outputs": {data_path.name: ARTIFACT_BUILDER.file_identity(data_path)},
    }
    receipt_path = data_dir / "receipt.json"
    DATA_BUILDER.write_json(receipt_path, receipt)
    (data_dir / "receipt.json.sha256").write_text(
        f"{ARTIFACT_BUILDER.sha256_path(receipt_path)}  receipt.json\n",
        encoding="ascii",
    )
    return data_path


def make_extraction_receipt_fixture(
    tmp_path: Path,
    *,
    dataset: str = "HCC1395",
    chrom: str = "chr1",
    write_receipt: bool = True,
    checks: dict | None = None,
) -> tuple[dict, Path]:
    extraction_dir = tmp_path / "extraction"
    extraction_dir.mkdir(parents=True)
    molecule_path = extraction_dir / f"{dataset}.{chrom}.molecule_sparse_calls.tsv.gz"
    catalog_path = extraction_dir / f"{dataset}.{chrom}.site_catalog.tsv.gz"
    extractor_path = tmp_path / "extractor.py"
    manifest_path = tmp_path / "manifest.json"
    molecule_path.write_bytes(b"molecules")
    catalog_path.write_bytes(b"catalog")
    extractor_path.write_text("# fixture\n", encoding="utf-8")
    manifest_path.write_text("{}\n", encoding="utf-8")
    molecule_identity = DATA_BUILDER.identity(molecule_path)
    catalog_identity = DATA_BUILDER.identity(catalog_path)
    strict_inputs = {
        "molecule_calls": molecule_identity,
        "site_catalog": catalog_identity,
    }
    receipt_path = extraction_dir / "receipt.json"
    if write_receipt:
        receipt = {
            "schema_name": "intersubmod.lossless_read_linkage_chromosome_receipt",
            "schema_version": "1.3.0",
            "all_pass": True,
            "scope": {"dataset": dataset, "chrom": chrom},
            "checks": (
                checks
                if checks is not None
                else {
                    name: True
                    for name in DATA_BUILDER.REQUIRED_EXTRACTION_RECEIPT_CHECKS
                }
            ),
            "provenance": {
                "extractor": {
                    "path": str(extractor_path.resolve()),
                    "sha256": DATA_BUILDER.sha256_path(extractor_path),
                },
                "manifest": {
                    "path": str(manifest_path.resolve()),
                    "sha256": DATA_BUILDER.sha256_path(manifest_path),
                },
            },
            "outputs": {
                molecule_path.name: molecule_identity,
                catalog_path.name: catalog_identity,
            },
        }
        DATA_BUILDER.write_json(receipt_path, receipt)
        (extraction_dir / "receipt.json.sha256").write_text(
            f"{DATA_BUILDER.sha256_path(receipt_path)}  receipt.json\n",
            encoding="ascii",
        )
    return strict_inputs, receipt_path


@pytest.mark.skipif(not (HCC_ROOT / "summary/summary.json").is_file(), reason="HCC strict SoT absent")
def test_hcc_distance_metrics_are_independently_recomputed():
    row = DATA_BUILDER.load_dataset("HCC1395", HCC_ROOT)["dataset_row"]

    assert row["candidate_loci_S"] == 79_687
    assert row["tree_eligible_W"] == 11_462
    assert row["W_span_gt_50kb"] == 1_064
    assert row["W_adjacent_gap_gt_50kb"] == 4
    assert row["W_span_gt_50kb_without_adjacent_gap_gt_50kb"] == 1_060
    assert row["direct_edges_gt_50kb"] == 47
    assert row["direct_edges_gt_50kb_support_3"] == 35
    assert row["direct_edges_gt_50kb_support_4"] == 12
    assert row["direct_edges_gt_50kb_support_ge_5"] == 0
    assert row["direct_edges_gt_50kb_RR_only"] == 24
    assert row["direct_edges_gt_50kb_alt_informative"] == 23
    assert row["W_with_direct_edge_gt_50kb"] == 22
    assert row["W_partition_changed_if_50kb_cut"] == 4
    assert row["W_fully_lost_if_50kb_cut"] == 0
    assert row["W_reduced_to_one_W_if_50kb_cut"] == 4
    assert row["W_split_to_multiple_W_if_50kb_cut"] == 0
    assert row["counterfactual_W_after_50kb_cut"] == 11_462
    assert row["counterfactual_W_delta_after_50kb_cut"] == 0
    assert row["W_memberships_lost_to_singletons_if_50kb_cut"] == 4


@pytest.mark.skipif(not (HCC_ROOT / "summary/summary.json").is_file(), reason="HCC strict SoT absent")
def test_artifact_contract_has_svg_charts_tables_and_full_154_row_scope(tmp_path: Path):
    hcc = DATA_BUILDER.load_dataset("HCC1395", HCC_ROOT)
    datasets = list(DATA_BUILDER.CANONICAL_DATASETS)

    dataset_rows = []
    chromosome_rows = []
    k_band_rows = []
    exact_k_rows = []
    span_band_rows = []
    threshold_rows = []
    for dataset in datasets:
        dataset_row = copy.deepcopy(hcc["dataset_row"])
        dataset_row["dataset"] = dataset
        dataset_rows.append(dataset_row)
        for source_row in hcc["chromosome_rows"]:
            row = copy.deepcopy(source_row)
            row["dataset"] = dataset
            chromosome_rows.append(row)
        for source_row in hcc["k_band_rows"]:
            row = copy.deepcopy(source_row)
            row["dataset"] = dataset
            k_band_rows.append(row)
        for source_row in hcc["k_rows"]:
            row = copy.deepcopy(source_row)
            row["dataset"] = dataset
            exact_k_rows.append(row)
        for source_row in hcc["span_band_rows"]:
            row = copy.deepcopy(source_row)
            row["dataset"] = dataset
            span_band_rows.append(row)
        for source_row in hcc["threshold_rows"]:
            row = copy.deepcopy(source_row)
            row["dataset"] = dataset
            threshold_rows.append(row)

    data = {
        "schema_name": DATA_BUILDER.SCHEMA_NAME,
        "schema_version": DATA_BUILDER.SCHEMA_VERSION,
        "created_at_utc": "2026-07-23T00:00:00Z",
        "all_pass": True,
        "scope": {
            "task_type": "B_comprehensive_validation",
            "datasets": datasets,
            "chromosomes": list(DATA_BUILDER.AUTOSOMES),
            "primary_minimum_distinct_molecules": 3,
        },
        "datasets": dataset_rows,
        "chromosomes": chromosome_rows,
        "k_bands": k_band_rows,
        "exact_k": exact_k_rows,
        "span_bands": span_band_rows,
        "threshold_sensitivity": threshold_rows,
        "topology_completion": [
            {
                "dataset": dataset,
                "strict_linkage_status": "COMPLETE_22_OF_22",
                "strict_W": next(
                    row["tree_eligible_W"]
                    for row in dataset_rows
                    if row["dataset"] == dataset
                ),
                "legacy_candidate_tree_status": "LEGACY_REFERENCE_ONLY",
                "exploratory_exact_ps_topology_status": (
                    "HCC1395_TECHNICAL_PILOT_NOT_ELIGIBLE"
                    if dataset == "HCC1395"
                    else "NOT_AVAILABLE"
                ),
                "strict_directed_topology_status": (
                    "NOT_RUN_LATEST_V4_PARTITION_ONLY"
                    if dataset == "HCC1395"
                    else "NOT_RUN"
                ),
                "strict_topology_reason": "fixture",
                "vaf_or_read_likelihood_ranking_status": "NOT_PRODUCTION_VALIDATED",
                "cellular_clone_count_status": "NOT_VALIDATED",
                "exact_cellular_parent_child_status": "NOT_VALIDATED",
                "cross_hp_or_technical_fusion_status": "NOT_VALIDATED",
            }
            for dataset in datasets
        ],
        "topology_completion_summary": {
            "strict_linkage_complete_datasets": 7,
            "strict_directed_topology_production_complete_datasets": 0,
            "cellular_clone_count_validated_datasets": 0,
            "exact_cellular_parent_child_validated_datasets": 0,
            "fused_tree_validated_datasets": 0,
        },
        "checks": {
            "fixture_all_pass": True,
            "all_154_extraction_receipt_identities_verified": True,
        },
        "inputs": {
            dataset: {
                "summary": {"path": f"/fixture/{dataset}/summary.json"},
                "chromosome_receipts": {
                    chrom: {"path": f"/fixture/{dataset}/{chrom}/receipt.json"}
                    for chrom in DATA_BUILDER.AUTOSOMES
                },
                "extraction_receipts": {
                    chrom: {
                        "path": f"/fixture/{dataset}/{chrom}/extraction/receipt.json",
                        "size_bytes": 1,
                        "sha256": "1" * 64,
                    }
                    for chrom in DATA_BUILDER.AUTOSOMES
                },
            }
            for dataset in datasets
        },
    }
    output_path = tmp_path / "artifact.json"
    data_path = write_data_package(tmp_path, data)
    artifact = ARTIFACT_BUILDER.build_artifact(
        data,
        data_path,
        output_path,
    )
    DATA_BUILDER.write_json(output_path, artifact)

    manifest = artifact["manifest"]
    snapshot = artifact["snapshot"]
    source_ids = {item["id"] for item in manifest["sources"]}
    chart_ids = {item["id"] for item in manifest["charts"]}
    chart_block_ids = {
        item["chartId"] for item in manifest["blocks"] if item["type"] == "chart"
    }
    assert manifest["title"] == ARTIFACT_BUILDER.TITLE
    assert manifest["blocks"][0]["body"].startswith(f"# {manifest['title']}\n")
    assert chart_ids == chart_block_ids
    assert len(chart_ids) == 11
    assert "<svg" in next(
        item["body"] for item in manifest["blocks"] if item["id"] == "method_flow"
    )
    for item in manifest["charts"]:
        assert item["sourceId"] in source_ids
        assert item["sourceId"].startswith("source_")
        assert item["encodings"]["x"]["field"]
        assert item["encodings"]["y"]["field"]
    for item in manifest["tables"]:
        declared = {column["field"] for column in item["columns"]}
        assert item["sourceId"] in source_ids
        assert item["sourceId"].startswith("source_")
        assert item["defaultSort"]["field"] in declared
        assert item["defaultSort"]["direction"] in {"asc", "desc"}
    assert snapshot["status"] == "ready"
    assert len(snapshot["datasets"]["dataset_summary"]) == 7
    assert len(snapshot["datasets"]["chromosome_detail"]) == 154
    assert len(snapshot["datasets"]["topology_completion"]) == 7
    assert {
        row["strict_directed_topology_status"]
        for row in snapshot["datasets"]["topology_completion"]
    } == {"NOT_RUN", "NOT_RUN_LATEST_V4_PARTITION_ONLY"}
    assert all(
        row["cellular_clone_count_status"] == "NOT_VALIDATED"
        and row["exact_cellular_parent_child_status"] == "NOT_VALIDATED"
        and row["cross_hp_or_technical_fusion_status"] == "NOT_VALIDATED"
        for row in snapshot["datasets"]["topology_completion"]
    )
    assert "topology_completion_table" in {
        item["id"] for item in manifest["tables"]
    }
    summary_text = next(
        item["body"]
        for item in manifest["blocks"]
        if item["id"] == "summary_text"
    )
    assert "strict directed topology 為 0/7" in summary_text
    assert "clone 數" in summary_text
    assert len(snapshot["datasets"]["exact_k"]) == 7 * len(hcc["k_rows"])
    for dataset in datasets:
        molecule_shares = [
            row["fraction_of_canonical_rows"]
            for row in snapshot["datasets"]["molecule_row_outcome"]
            if row["dataset"] == dataset
        ]
        component_shares = [
            row["fraction_of_components"]
            for row in snapshot["datasets"]["component_outcome"]
            if row["dataset"] == dataset
        ]
        k_shares = [
            row["fraction_of_W"]
            for row in snapshot["datasets"]["k_bands"]
            if row["dataset"] == dataset
        ]
        span_shares = [
            row["fraction_of_W"]
            for row in snapshot["datasets"]["span_bands"]
            if row["dataset"] == dataset
        ]
        assert sum(molecule_shares) == pytest.approx(1.0)
        assert sum(component_shares) == pytest.approx(1.0)
        assert sum(k_shares) == pytest.approx(1.0)
        assert sum(span_shares) == pytest.approx(1.0)
        assert [
            row["k_band"]
            for row in snapshot["datasets"]["k_bands"]
            if row["dataset"] == dataset
        ] == ["2", "3", "4", "5", "6–8", "9–12", ">12"]
        assert [
            row["span_band"]
            for row in snapshot["datasets"]["span_bands"]
            if row["dataset"] == dataset
        ] == ["≤10 kb", "10–25 kb", "25–50 kb", "50–100 kb", ">100 kb"]
    grouped_composition_charts = {
        "molecule_row_outcome_chart",
        "component_outcome_chart",
        "k_distribution_chart",
        "span_distribution_chart",
        "long_edge_support_chart",
        "long_edge_state_chart",
    }
    assert {
        item["id"]
        for item in manifest["charts"]
        if item["options"].get("grouping") == "grouped"
    }.issuperset(grouped_composition_charts)
    assert all(
        isinstance(rows, list) and len(rows) <= 2_000
        for rows in snapshot["datasets"].values()
    )
    assert (tmp_path / "data/report.sqlite").is_file()
    assert (tmp_path / "queries/dataset_summary.sql").read_text(
        encoding="utf-8"
    ).startswith("SELECT ")
    dataset_source = next(
        source
        for source in manifest["sources"]
        if source["id"] == "source_dataset_summary"
    )
    assert "data/report.sqlite" in dataset_source["query"]["description"]
    assert {
        row["validation_basis"]
        for row in threshold_rows
        if row["minimum_distinct_molecules"] == 3
    } == {"independent_row_recomputation_and_receipt_reconciliation"}
    assert {
        row["validation_basis"]
        for row in threshold_rows
        if row["minimum_distinct_molecules"] in {1, 2, 5}
    } == {"chromosome_receipt_reaggregation_and_summary_reconciliation"}

    published_root = tmp_path / "published"
    published_data_path = write_data_package(published_root, data)
    published_artifact_path = published_root / "artifact.json"
    assert ARTIFACT_BUILDER.main(
        [
            "--data",
            str(published_data_path),
            "--output",
            str(published_artifact_path),
        ]
    ) == 0
    assert published_artifact_path.is_file()
    artifact_receipt_path = published_root / "artifact_receipt.json"
    artifact_sidecar_path = published_root / "artifact_receipt.json.sha256"
    assert artifact_receipt_path.is_file()
    assert artifact_sidecar_path.read_text(encoding="ascii").strip() == (
        f"{ARTIFACT_BUILDER.sha256_path(artifact_receipt_path)}  "
        "artifact_receipt.json"
    )
    artifact_receipt = ARTIFACT_BUILDER.read_json(artifact_receipt_path)
    assert artifact_receipt["all_pass"] is True
    assert all(artifact_receipt["checks"].values())
    assert artifact_receipt["outputs"]["artifact.json"] == (
        ARTIFACT_BUILDER.file_identity(published_artifact_path)
    )
    assert (published_root / "data/report.sqlite").is_file()
    assert (published_root / "queries/dataset_summary.sql").is_file()
    assert not list(published_root.glob(".artifact.json.pending-*"))

    empty_checks = copy.deepcopy(data)
    empty_checks["checks"] = {}
    with pytest.raises(ValueError, match="missing or failed check"):
        ARTIFACT_BUILDER.validate_report_data_package(empty_checks, data_path)

    missing_extraction_grid = copy.deepcopy(data)
    del missing_extraction_grid["inputs"]["HCC1395"]["extraction_receipts"]["chr1"]
    with pytest.raises(ValueError, match="input provenance is incomplete"):
        ARTIFACT_BUILDER.validate_report_data_package(
            missing_extraction_grid,
            data_path,
        )

    wrong_threshold_basis = copy.deepcopy(data)
    wrong_threshold_basis["threshold_sensitivity"][0]["validation_basis"] = "placeholder"
    with pytest.raises(ValueError, match="validation_basis"):
        ARTIFACT_BUILDER.validate_report_data_package(
            wrong_threshold_basis,
            data_path,
        )

    overstated_topology = copy.deepcopy(data)
    overstated_topology["topology_completion"][0][
        "strict_directed_topology_status"
    ] = "COMPLETE_PRODUCTION"
    with pytest.raises(ValueError, match="topology completion boundary"):
        ARTIFACT_BUILDER.validate_report_data_package(
            overstated_topology,
            data_path,
        )

    data_path.write_text(
        data_path.read_text(encoding="utf-8") + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="identity mismatch"):
        ARTIFACT_BUILDER.validate_report_data_package(data, data_path)


def test_primary_container_validator_rejects_invalid_basis_and_missing_ps(tmp_path: Path):
    valid = {
        "dataset": "HCC1395",
        "chrom": "chr1",
        "linkage_basis": "PS_HP1",
        "phase_set": "123",
    }
    assert DATA_BUILDER.require_exact_primary_container(
        dataset="HCC1395",
        chrom="chr1",
        row=valid,
        path=tmp_path / "row.tsv.gz",
    ) == ("PS_HP1", "123")
    invalid_basis = {**valid, "linkage_basis": "HP3"}
    with pytest.raises(ValueError, match="linkage_basis"):
        DATA_BUILDER.require_exact_primary_container(
            dataset="HCC1395",
            chrom="chr1",
            row=invalid_basis,
            path=tmp_path / "row.tsv.gz",
        )
    missing_ps = {**valid, "phase_set": "."}
    with pytest.raises(ValueError, match="phase_set"):
        DATA_BUILDER.require_exact_primary_container(
            dataset="HCC1395",
            chrom="chr1",
            row=missing_ps,
            path=tmp_path / "row.tsv.gz",
        )


def test_extraction_receipt_provenance_accepts_complete_fixture(tmp_path: Path):
    strict_inputs, receipt_path = make_extraction_receipt_fixture(tmp_path)
    assert DATA_BUILDER.validate_extraction_receipt(
        dataset="HCC1395",
        chrom="chr1",
        strict_inputs=strict_inputs,
    ) == DATA_BUILDER.identity(receipt_path)


def test_extraction_receipt_provenance_rejects_missing_receipt(tmp_path: Path):
    strict_inputs, _ = make_extraction_receipt_fixture(
        tmp_path,
        write_receipt=False,
    )
    with pytest.raises(ValueError, match="receipt or checksum sidecar is missing"):
        DATA_BUILDER.validate_extraction_receipt(
            dataset="HCC1395",
            chrom="chr1",
            strict_inputs=strict_inputs,
        )


def test_extraction_receipt_provenance_rejects_tampered_receipt(tmp_path: Path):
    strict_inputs, receipt_path = make_extraction_receipt_fixture(tmp_path)
    receipt_path.write_text(
        receipt_path.read_text(encoding="utf-8") + " ",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="checksum mismatch"):
        DATA_BUILDER.validate_extraction_receipt(
            dataset="HCC1395",
            chrom="chr1",
            strict_inputs=strict_inputs,
        )


def test_extraction_receipt_provenance_rejects_placeholder_checks(tmp_path: Path):
    strict_inputs, _ = make_extraction_receipt_fixture(
        tmp_path,
        checks={"placeholder": True},
    )
    with pytest.raises(ValueError, match="canonical extraction checks"):
        DATA_BUILDER.validate_extraction_receipt(
            dataset="HCC1395",
            chrom="chr1",
            strict_inputs=strict_inputs,
        )


def test_extraction_receipt_provenance_rejects_strict_input_identity_drift(
    tmp_path: Path,
):
    strict_inputs, _ = make_extraction_receipt_fixture(tmp_path)
    strict_inputs["molecule_calls"] = {
        **strict_inputs["molecule_calls"],
        "sha256": "0" * 64,
    }
    with pytest.raises(ValueError, match="differs from strict input"):
        DATA_BUILDER.validate_extraction_receipt(
            dataset="HCC1395",
            chrom="chr1",
            strict_inputs=strict_inputs,
        )


@pytest.mark.skipif(not (HCC_ROOT / "summary/summary.json").is_file(), reason="HCC strict SoT absent")
def test_dataset_loader_rejects_summary_to_chromosome_receipt_identity_mismatch(
    monkeypatch: pytest.MonkeyPatch,
):
    original_identity = DATA_BUILDER.identity

    def tampered_identity(path: Path):
        observed = original_identity(path)
        if path.name == "receipt.json" and path.parent.name == "chr1":
            observed = {**observed, "sha256": "0" * 64}
        return observed

    monkeypatch.setattr(DATA_BUILDER, "identity", tampered_identity)
    with pytest.raises(ValueError, match="input identity differs"):
        DATA_BUILDER.load_dataset("HCC1395", HCC_ROOT)


@pytest.mark.skipif(not (HCC_ROOT / "summary/summary.json").is_file(), reason="HCC strict SoT absent")
@pytest.mark.parametrize("target", ["summary", "receipt"])
def test_dataset_loader_rejects_placeholder_check_maps(
    monkeypatch: pytest.MonkeyPatch,
    target: str,
):
    original_read_json = DATA_BUILDER.read_json

    def tampered_read_json(path: Path):
        observed = original_read_json(path)
        is_summary = path.name == "summary.json"
        is_first_receipt = path.name == "receipt.json" and path.parent.name == "chr1"
        if (target == "summary" and is_summary) or (
            target == "receipt" and is_first_receipt
        ):
            observed = {**observed, "checks": {"placeholder": True}}
        return observed

    monkeypatch.setattr(DATA_BUILDER, "read_json", tampered_read_json)
    expected = "canonical summary checks" if target == "summary" else "canonical chromosome checks"
    with pytest.raises(ValueError, match=expected):
        DATA_BUILDER.load_dataset("HCC1395", HCC_ROOT)
