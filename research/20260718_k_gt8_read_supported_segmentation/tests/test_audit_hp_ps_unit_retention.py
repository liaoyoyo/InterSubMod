import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys

import pytest


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from audit_hp_ps_unit_retention import (  # noqa: E402
    AuditError,
    aggregate_component_pairs,
    aggregate_units,
    load_constraints,
    load_membership,
    quantile_type7,
    run,
    summarize_rows,
    verify_partition_source,
)


CHR22_PROBE = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260718_k_gt8_read_supported_segmentation/probes/HCC1395_chr22/"
    "partition_v2"
)
CHR6_FULL_PARTITION = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260718_k_gt8_read_supported_segmentation/full/HCC1395_chr1_22_v1/"
    "chromosomes/chr6/partition"
)


def _sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _identity(path):
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": _sha(path),
    }


def _write_gzip_tsv(path, fields, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _constraint(
    component,
    constraint_id,
    hp,
    phase_set,
    positions,
    weight,
    disposition="retained",
):
    return {
        "dataset": "HCC1395",
        "chrom": "chr22",
        "legacy_component_id": component,
        "constraint_id": constraint_id,
        "unit_id": f"HP{hp}|PS{phase_set}",
        "hp_family": hp,
        "phase_set": phase_set,
        "positions": ",".join(map(str, positions)),
        "call_codes": "R" * len(positions),
        "n_fixed_ra": len(positions),
        "span_sites": positions[-1] - positions[0] + 1,
        "molecule_weight": weight,
        "disposition": disposition,
        "crossed_cut_count": 0 if disposition == "retained" else 1,
        "retained_block_index": 1 if disposition == "retained" else "",
    }


def _fixture_partition(tmp_path):
    partition = tmp_path / "partition"
    partition.mkdir()
    component1 = "HCC1395:chr22:legacy_gap_50000:1:100-108"
    component2 = "HCC1395:chr22:legacy_gap_50000:2:200-208"
    membership_fields = (
        "dataset",
        "chrom",
        "legacy_component_id",
        "site_index",
        "component_local_index",
        "pos1",
        "ref",
        "alt",
        "block_id",
        "primary_linkage_observed",
        "old_densest8_disposition",
    )
    active = {100, 101, 102, 103, 104, 105, 106, 108}
    membership_rows = []
    site_index = 0
    for component, start in ((component1, 100), (component2, 200)):
        for local_index in range(9):
            pos = start + local_index
            membership_rows.append(
                {
                    "dataset": "HCC1395",
                    "chrom": "chr22",
                    "legacy_component_id": component,
                    "site_index": site_index,
                    "component_local_index": local_index,
                    "pos1": pos,
                    "ref": "A",
                    "alt": "G",
                    "block_id": f"{component}:B1",
                    "primary_linkage_observed": str(
                        component == component1 and pos in active
                    ).lower(),
                    "old_densest8_disposition": (
                        "selected" if local_index < 8 else "cap_excluded"
                    ),
                }
            )
            site_index += 1
    membership_path = partition / "HCC1395.chr22.site_membership.tsv.gz"
    _write_gzip_tsv(membership_path, membership_fields, membership_rows)

    constraints = [
        _constraint(component1, "P01", "1", "1000", (100, 101), 10),
        _constraint(component1, "P02", "1", "1000", (101, 102), 5),
        _constraint(
            component1, "P03", "1", "1000", (102, 103), 3, "cut"
        ),
        _constraint(
            component1,
            "P04",
            "1",
            "1000",
            (100, 108),
            1,
            "unavoidable_span_gt_max_block_size",
        ),
        _constraint(component1, "P05", "1", "1000", (103, 104), 1),
        _constraint(component1, "P06", "2", "1000", (100, 101), 5),
        _constraint(component1, "P07", "2", "1000", (101, 102), 5),
        _constraint(component1, "P08", "2", "1000", (102, 103), 5),
        _constraint(component1, "P09", "2", "1000", (103, 104), 5),
        _constraint(component1, "P10", "2", "1000", (104, 105), 5),
        _constraint(component1, "P11", "1", "2000", (105, 106), 1),
    ]
    constraint_fields = tuple(constraints[0])
    constraint_path = partition / "HCC1395.chr22.cut_constraints.tsv.gz"
    _write_gzip_tsv(constraint_path, constraint_fields, constraints)

    legacy_fields = (
        "dataset",
        "chrom",
        "legacy_component_id",
        "pre_cap_k",
        "primary_active_site_count",
        "exact_pattern_count",
        "raw_total_molecule_weight",
        "raw_retained_molecule_weight",
        "raw_lost_molecule_weight",
        "retained_exact_pattern_count",
        "lost_exact_pattern_count",
        "unavoidable_pattern_count",
    )
    legacy_path = partition / "HCC1395.chr22.legacy_components.tsv.gz"
    _write_gzip_tsv(
        legacy_path,
        legacy_fields,
        [
            {
                "dataset": "HCC1395",
                "chrom": "chr22",
                "legacy_component_id": component1,
                "pre_cap_k": 9,
                "primary_active_site_count": 8,
                "exact_pattern_count": 11,
                "raw_total_molecule_weight": 46,
                "raw_retained_molecule_weight": 42,
                "raw_lost_molecule_weight": 4,
                "retained_exact_pattern_count": 9,
                "lost_exact_pattern_count": 2,
                "unavoidable_pattern_count": 1,
            },
            {
                "dataset": "HCC1395",
                "chrom": "chr22",
                "legacy_component_id": component2,
                "pre_cap_k": 9,
                "primary_active_site_count": 0,
                "exact_pattern_count": 0,
                "raw_total_molecule_weight": 0,
                "raw_retained_molecule_weight": 0,
                "raw_lost_molecule_weight": 0,
                "retained_exact_pattern_count": 0,
                "lost_exact_pattern_count": 0,
                "unavoidable_pattern_count": 0,
            },
        ],
    )
    blocks_path = partition / "HCC1395.chr22.blocks.tsv.gz"
    _write_gzip_tsv(
        blocks_path,
        ("dataset", "chrom", "legacy_component_id", "k"),
        [
            {
                "dataset": "HCC1395",
                "chrom": "chr22",
                "legacy_component_id": component1,
                "k": 8,
            },
            {
                "dataset": "HCC1395",
                "chrom": "chr22",
                "legacy_component_id": component2,
                "k": 8,
            },
        ],
    )
    site_catalog = tmp_path / "site_catalog.tsv.gz"
    molecule_calls = tmp_path / "molecule_calls.tsv.gz"
    _write_gzip_tsv(site_catalog, ("x",), [{"x": 1}])
    _write_gzip_tsv(molecule_calls, ("x",), [{"x": 1}])
    receipt = {
        "schema_name": "intersubmod.k_gt8_read_supported_segmentation",
        "schema_version": "0.1.0",
        "all_pass": True,
        "scope": {
            "dataset": "HCC1395",
            "chrom": "chr22",
            "site_catalog_sites": 18,
        },
        "parameters": {
            "primary_hp_families": ["1", "2"],
            "require_exact_known_phase_set": True,
            "max_block_size": 8,
        },
        "counts": {
            "target_components": 2,
            "target_sites": 18,
            "exact_patterns": 11,
            "raw_total_molecule_weight": 46,
            "raw_retained_molecule_weight": 42,
            "raw_lost_molecule_weight": 4,
            "unavoidable_patterns": 1,
            "primary_active_sites_component_sum": 8,
        },
        "checks": {
            "constraint_molecule_mass_conserved": True,
            "constraint_rows_equal_exact_patterns": True,
            "hp_ps_columns_nonempty_and_isolated": True,
            "target_component_count_matches_expected": True,
            "target_site_count_matches_expected": True,
            "target_sites_assigned_once": True,
        },
        "inputs": {
            "site_catalog": _identity(site_catalog),
            "molecule_calls": _identity(molecule_calls),
        },
        "outputs": {
            "legacy_components": _identity(legacy_path),
            "blocks": _identity(blocks_path),
            "site_membership": _identity(membership_path),
            "cut_constraints": _identity(constraint_path),
        },
    }
    receipt_path = partition / "receipt.json"
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (partition / "receipt.json.sha256").write_text(
        f"{_sha(receipt_path)}  receipt.json\n", encoding="utf-8"
    )
    return partition, component1, component2


def _load_fixture(tmp_path):
    partition, component1, component2 = _fixture_partition(tmp_path)
    source = verify_partition_source(partition, {"mode": "probe"})
    components = load_membership(source)
    constraints, active = load_constraints(source, components)
    units = aggregate_units(components, constraints)
    return source, components, constraints, active, units, component1, component2


def test_unit_aggregation_is_hp_ps_isolated_and_mass_conserving(tmp_path):
    _, components, constraints, _, units, component1, _ = _load_fixture(
        tmp_path
    )
    assert len(constraints) == 11
    assert len(units) == 3
    by_key = {(row["hp_family"], row["phase_set"]): row for row in units}
    hp1 = by_key[("1", "1000")]
    assert hp1["legacy_component_id"] == component1
    assert hp1["component_k"] == 9
    assert hp1["total_pattern_rows"] == 5
    assert hp1["retained_pattern_rows"] == 3
    assert hp1["cut_lost_pattern_rows"] == 1
    assert hp1["unavoidable_pattern_rows"] == 1
    assert hp1["total_molecule_component_incidence_weight"] == 20
    assert hp1["retained_molecule_component_incidence_weight"] == 16
    assert hp1["retention_ratio"] == "0.800000000000"
    assert hp1["support_stratum"] == "20-49"
    assert hp1["eligible_headline"] == "true"
    assert components[component1].k == 9


def test_no_constraint_component_is_not_fabricated_as_zero_or_one(tmp_path):
    _, components, _, _, units, _, component2 = _load_fixture(tmp_path)
    pairs = aggregate_component_pairs(units)
    summary = summarize_rows(
        units,
        pairs,
        source_components={"chr22": set(components)},
        skipped_chromosomes=[],
    )
    assert all(row["legacy_component_id"] != component2 for row in units)
    assert summary["counts"]["components_in_partition_scope"] == 2
    assert summary["counts"]["components_with_observed_constraint_units"] == 1
    assert summary["counts"]["components_without_observed_constraint_units"] == 1
    assert summary["component_hp_unit_coverage"][
        "components_hp1_and_hp2"
    ] == 1
    assert summary["component_hp_unit_coverage"][
        "components_without_observed_unit"
    ] == 1
    assert summary["component_hp_unit_coverage"]["by_chromosome"]["chr22"][
        "components_hp1_and_hp2"
    ] == 1
    assert summary["by_chromosome"]["chr22"][
        "observed_constraint_units"
    ] == 3
    assert summary["by_hp_family"]["HP1"]["observed_constraint_units"] == 2
    assert summary["by_hp_family"]["HP2"]["observed_constraint_units"] == 1
    assert summary["by_support_stratum"]["20-49"][
        "observed_constraint_units"
    ] == 2
    assert summary["scope_contract"]["scope_ceiling"] == (
        "observed_constraint_units_only"
    )


def test_component_pair_sums_phase_sets_within_hp_before_delta(tmp_path):
    _, _, _, _, units, component1, _ = _load_fixture(tmp_path)
    pairs = aggregate_component_pairs(units)
    assert len(pairs) == 1
    pair = pairs[0]
    assert pair["legacy_component_id"] == component1
    assert pair["hp1_phase_set_unit_count"] == 2
    assert pair["hp1_total_pattern_rows"] == 6
    assert pair["hp1_total_molecule_component_incidence_weight"] == 21
    assert pair["hp1_retained_molecule_component_incidence_weight"] == 17
    assert pair["hp2_total_molecule_component_incidence_weight"] == 25
    assert pair["hp1_minus_hp2_retention_delta"] == "-0.190476190476"
    assert pair["both_hp_headline_eligible"] == "true"


def test_unknown_disposition_fails_closed(tmp_path):
    partition, _, _ = _fixture_partition(tmp_path)
    source = verify_partition_source(partition, {"mode": "probe"})
    components = load_membership(source)
    constraints_path = partition / "HCC1395.chr22.cut_constraints.tsv.gz"
    with gzip.open(constraints_path, "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows[0]["disposition"] = "mystery"
    _write_gzip_tsv(
        tmp_path / "bad.tsv.gz",
        tuple(rows[0]),
        rows,
    )
    bad_source = source.__class__(
        dataset=source.dataset,
        chrom=source.chrom,
        partition_dir=source.partition_dir,
        receipt=source.receipt,
        receipt_identity=source.receipt_identity,
        legacy_components_path=source.legacy_components_path,
        membership_path=source.membership_path,
        constraints_path=tmp_path / "bad.tsv.gz",
        source_context=source.source_context,
    )
    with pytest.raises(AuditError, match="unknown constraint disposition"):
        load_constraints(bad_source, components)


def test_declared_output_hash_drift_fails_before_parsing(tmp_path):
    partition, _, _ = _fixture_partition(tmp_path)
    membership = partition / "HCC1395.chr22.site_membership.tsv.gz"
    with membership.open("ab") as handle:
        handle.write(b"drift")
    with pytest.raises(AuditError, match="identity drift"):
        verify_partition_source(partition, {"mode": "probe"})


def test_type7_quantile_is_exact_linear_interpolation():
    from decimal import Decimal

    values = [Decimal("0"), Decimal("0.5"), Decimal("1")]
    assert quantile_type7(values, Decimal("0.25")) == Decimal("0.25")
    assert quantile_type7(values, Decimal("0.75")) == Decimal("0.75")


def test_end_to_end_probe_writes_receipt_and_refuses_existing_output(tmp_path):
    partition, _, _ = _fixture_partition(tmp_path)
    args = type(
        "Args",
        (),
        {
            "mode": "probe",
            "source_root": None,
            "partition_dir": [partition],
            "output_dir": tmp_path / "audit",
        },
    )()
    receipt = run(args)
    assert receipt["all_pass"] is True
    assert receipt["summary"]["counts"]["observed_constraint_units"] == 3
    assert (args.output_dir / "HCC1395.hp_ps_observed_units.tsv.gz").is_file()
    assert (args.output_dir / "HCC1395.hp1_hp2_paired_components.tsv.gz").is_file()
    assert (args.output_dir / "summary.json").is_file()
    assert (args.output_dir / "summary.tsv").is_file()
    assert (args.output_dir / "receipt.json.sha256").is_file()
    with pytest.raises(AuditError, match="must not exist"):
        run(args)


@pytest.mark.parametrize("partition_dir", [CHR22_PROBE, CHR6_FULL_PARTITION])
def test_authenticated_real_partition_source_passes_independent_contract(
    partition_dir,
):
    if not partition_dir.is_dir():
        pytest.skip(f"authenticated source unavailable: {partition_dir}")
    source = verify_partition_source(partition_dir, {"mode": "probe"})
    components = load_membership(source)
    constraints, active = load_constraints(source, components)
    units = aggregate_units(components, constraints)
    assert constraints
    assert units
    assert sum(len(values) for values in active.values()) == int(
        source.receipt["counts"]["primary_active_sites_component_sum"]
    )
    assert all(row["ratio_status"] == "OBSERVED_CONSTRAINT_DENOMINATOR" for row in units)
