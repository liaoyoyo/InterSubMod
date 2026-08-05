from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import pytest


MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "build_strict_ps_hp_regions.py"
SPEC = importlib.util.spec_from_file_location("test_build_strict_ps_hp_regions", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
BUILDER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = BUILDER
SPEC.loader.exec_module(BUILDER)

PRE_HOTFIX_BUILDER_SHA256 = (
    "7260a7631b30cbb5e4878159583b8a5b27a153de07e8d001303417dd2f29aedd"
)
PRE_HOTFIX_BUILDER_SNAPSHOT = (
    MODULE_PATH.parents[1]
    / "research/20260723_production_exact_ps_strict_read_linkage/source_snapshots"
    / f"{PRE_HOTFIX_BUILDER_SHA256}_build_strict_ps_hp_regions.py"
)


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_pre_hotfix_builder_snapshot_identity_is_preserved() -> None:
    assert hashlib.sha256(PRE_HOTFIX_BUILDER_SNAPSHOT.read_bytes()).hexdigest() == (
        PRE_HOTFIX_BUILDER_SHA256
    )


SITE_FIELDS = ("site_index", "chrom", "pos1", "ref", "alt")
MOLECULE_FIELDS = (
    "dataset",
    "chrom",
    "molecule_id",
    "hp_family",
    "phase_set",
    "site_indices",
    "positions1",
    "call_codes",
)


def site_rows() -> list[dict[str, object]]:
    return [
        {"site_index": index, "chrom": "chr1", "pos1": pos1, "ref": "C", "alt": "T"}
        for index, pos1 in enumerate((100, 200, 300, 100_500))
    ]


def molecule_row(
    molecule_id: str,
    phase_set: str,
    indices: str,
    positions: str,
    codes: str,
    *,
    hp: str = "1",
) -> dict[str, object]:
    return {
        "dataset": "SAMPLE",
        "chrom": "chr1",
        "molecule_id": molecule_id,
        "hp_family": hp,
        "phase_set": phase_set,
        "site_indices": indices,
        "positions1": positions,
        "call_codes": codes,
    }


def synthetic_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    rows.extend(
        molecule_row(f"ps10-ab-{index}", "10", "0,1", "100,200", "AR")
        for index in range(3)
    )
    rows.extend(
        molecule_row(f"ps10-bc-{index}", "10", "1,2", "200,300", "AR")
        for index in range(3)
    )
    # PS20: A-C has direct support while B is observed on a different molecule.
    # B must remain a singleton; X cannot create A-B or B-C edges.
    rows.extend(
        molecule_row(f"ps20-ac-{index}", "20", "0,1,2", "100,200,300", "AXA")
        for index in range(3)
    )
    rows.append(molecule_row("ps20-b", "20", "1", "200", "R"))
    # Distance is not a boundary when three molecules directly observe endpoints.
    rows.extend(
        molecule_row(f"ps30-far-{index}", "30", "0,3", "100,100500", "RA")
        for index in range(3)
    )
    for index, missing_token in enumerate(("", ".", "NA", "N/A", "NaN", "None", " NULL ")):
        rows.append(
            molecule_row(
                f"missing-ps-{index}",
                missing_token,
                "0,1",
                "100,200",
                "AA",
            )
        )
    rows.append(molecule_row("hp3", "10", "0,1", "100,200", "AA", hp="3"))
    return rows


def test_builds_exact_ps_components_from_direct_endpoint_reads(tmp_path: Path) -> None:
    sites = tmp_path / "sites.tsv.gz"
    molecules = tmp_path / "molecules.tsv.gz"
    output = tmp_path / "output"
    write_tsv(sites, SITE_FIELDS, site_rows())
    write_tsv(molecules, MOLECULE_FIELDS, synthetic_rows())

    receipt = BUILDER.build_regions(
        dataset="SAMPLE",
        chrom="chr1",
        site_catalog=sites,
        molecule_calls=molecules,
        output_dir=output,
        thresholds=(1, 2, 3, 5),
        primary_threshold=3,
    )

    assert receipt["all_pass"] is True
    assert receipt["counts"]["exact_ps_hp_containers"] == 3
    assert receipt["counts"]["excluded_missing_ps_rows"] == 7
    assert receipt["counts"]["excluded_nonprimary_hp_rows"] == 1
    summary = receipt["summary_by_threshold"]["3"]
    assert summary["regions"] == 4
    assert summary["k1_regions"] == 1
    assert summary["k_gt1_regions"] == 3

    components = read_tsv(output / "SAMPLE.chr1.components.tsv.gz")
    primary = [row for row in components if row["threshold"] == "3"]
    by_ps = {}
    for row in primary:
        by_ps.setdefault(row["phase_set"], []).append(row)
    assert sorted(int(row["k"]) for row in by_ps["10"]) == [3]
    assert sorted(int(row["k"]) for row in by_ps["20"]) == [1, 2]
    assert sorted(int(row["k"]) for row in by_ps["30"]) == [2]
    assert by_ps["30"][0]["contains_gap_gt_50kb"] == "true"
    assert all(row["linkage_rule"] == "strict_fixed_ra_endpoint_pair" for row in primary)
    singleton = next(row for row in by_ps["20"] if row["k"] == "1")
    assert singleton["inference_role"] == "ABSTAIN_SINGLETON_UNLINKED"
    assert singleton["tree_eligible"] == "false"
    assert singleton["solver_route"] == "EXCLUDE_SINGLETON_NO_READ_LINKAGE"
    linked = [row for row in primary if int(row["k"]) > 1]
    assert all(row["inference_role"] == "PRIMARY_PS_AWARE" for row in linked)
    assert all(row["tree_eligible"] == "true" for row in linked)

    edges = read_tsv(output / "SAMPLE.chr1.endpoint_edges.tsv.gz")
    ps20_pairs = {
        (row["site_i_index"], row["site_j_index"])
        for row in edges
        if row["phase_set"] == "20" and row["passes_primary_threshold"] == "true"
    }
    assert ps20_pairs == {("0", "2")}
    assert all(
        int(row["support_total"])
        == sum(int(row[key]) for key in ("support_RR", "support_RA", "support_AR", "support_AA"))
        for row in edges
    )

    persisted = json.loads((output / "receipt.json").read_text(encoding="utf-8"))
    assert persisted["checks"]["distance_not_used_for_connectivity"] is True
    assert (output / "receipt.json.sha256").is_file()


def test_duplicate_molecule_id_fails_closed_even_across_containers(tmp_path: Path) -> None:
    sites = tmp_path / "sites.tsv.gz"
    molecules = tmp_path / "molecules.tsv.gz"
    write_tsv(sites, SITE_FIELDS, site_rows())
    rows = [
        molecule_row("same", "10", "0,1", "100,200", "AA"),
        molecule_row("same", "20", "0,1", "100,200", "AA"),
    ]
    write_tsv(molecules, MOLECULE_FIELDS, rows)

    with pytest.raises(ValueError, match="duplicate molecule_id"):
        BUILDER.build_regions(
            dataset="SAMPLE",
            chrom="chr1",
            site_catalog=sites,
            molecule_calls=molecules,
            output_dir=tmp_path / "output",
            thresholds=(3,),
            primary_threshold=3,
        )
