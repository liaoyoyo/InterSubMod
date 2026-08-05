from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys


REPO = Path(__file__).resolve().parents[3]
SCRIPT_DIR = (
    REPO
    / "docs"
    / "methodology"
    / "_assets"
    / "20260627_subclone_4axis_teaching"
    / "scripts"
)
ADAPTER = SCRIPT_DIR / "exact_ps_partition_to_mlhp.py"
LAYERED = SCRIPT_DIR / "layered_tree_reconstruction.py"
REGION_VIEW = SCRIPT_DIR / "build_region_view.py"

ADAPTER_SPEC = importlib.util.spec_from_file_location(
    "exact_ps_partition_to_mlhp", ADAPTER
)
assert ADAPTER_SPEC is not None and ADAPTER_SPEC.loader is not None
ADAPTER_MODULE = importlib.util.module_from_spec(ADAPTER_SPEC)
ADAPTER_SPEC.loader.exec_module(ADAPTER_MODULE)


BLOCK_FIELDS = (
    "dataset", "chrom", "unit_id", "component_id", "linkage_basis",
    "hp_family", "phase_set", "threshold", "block_id", "block_index",
    "start1", "end1", "k", "positions1", "retained_molecule_weight",
    "retained_pattern_count",
)
DISPOSITION_FIELDS = (
    "dataset", "chrom", "unit_id", "constraint_id", "hp_family",
    "phase_set", "positions1", "call_codes", "molecule_weight",
    "pattern_count", "disposition", "span_sites", "crossed_cut_count",
    "retained_block_index",
)
MOLECULE_FIELDS = (
    "dataset", "chrom", "molecule_id", "hp_family", "phase_set",
    "positions1", "call_codes",
)


def test_receipt_size_bytes_accepts_both_schema_spellings() -> None:
    assert ADAPTER_MODULE.receipt_size_bytes({"size_bytes": 17}) == 17
    assert ADAPTER_MODULE.receipt_size_bytes({"bytes": 17}) == 17
    assert ADAPTER_MODULE.receipt_size_bytes({"size_bytes": True}) is None
    assert ADAPTER_MODULE.receipt_size_bytes({}) is None


def write_tsv_gz(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fixture(root: Path, *, wrong_ps: bool = False) -> Path:
    chrom_dir = root / "chromosomes" / "chr22"
    partition = chrom_dir / "python_partition"
    partition.mkdir(parents=True)
    (chrom_dir / "comparison.json").write_text(
        json.dumps({"all_pass": True, "mismatch_count": 0}), encoding="utf-8"
    )
    (partition / "receipt.json").write_text(
        json.dumps(
            {
                "all_pass": True,
                "scope": {"dataset": "HCC1395", "chrom": "chr22"},
                "checks": {
                    "cross_ps_zero": True,
                    "cross_hp_zero": True,
                    "constraint_mass_conserved": True,
                },
            }
        ),
        encoding="utf-8",
    )
    blocks = [
        {
            "dataset": "HCC1395", "chrom": "chr22", "unit_id": "U1",
            "component_id": "CP1", "linkage_basis": "PS_HP1", "hp_family": "1",
            "phase_set": "100", "threshold": 3, "block_id": "U1:B0001",
            "block_index": 1, "start1": 10, "end1": 30, "k": 3,
            "positions1": "10,20,30", "retained_molecule_weight": 9,
            "retained_pattern_count": 3,
        },
        {
            "dataset": "HCC1395", "chrom": "chr22", "unit_id": "U2",
            "component_id": "CP2", "linkage_basis": "PS_HP1", "hp_family": "1",
            "phase_set": "200", "threshold": 3, "block_id": "U2:B0001",
            "block_index": 1, "start1": 10, "end1": 30, "k": 3,
            "positions1": "10,20,30", "retained_molecule_weight": 3,
            "retained_pattern_count": 1,
        },
    ]
    dispositions = [
        {
            "dataset": "HCC1395", "chrom": "chr22", "unit_id": "U1",
            "constraint_id": "C1", "hp_family": "1",
            "phase_set": "999" if wrong_ps else "100", "positions1": "10,20",
            "call_codes": "AA", "molecule_weight": 3, "pattern_count": 1,
            "disposition": "retained", "span_sites": 2, "crossed_cut_count": 0,
            "retained_block_index": 1,
        },
        {
            "dataset": "HCC1395", "chrom": "chr22", "unit_id": "U1",
            "constraint_id": "C2", "hp_family": "1", "phase_set": "100",
            "positions1": "10,20,30", "call_codes": "RAA", "molecule_weight": 4,
            "pattern_count": 1, "disposition": "retained", "span_sites": 3,
            "crossed_cut_count": 0, "retained_block_index": 1,
        },
        {
            "dataset": "HCC1395", "chrom": "chr22", "unit_id": "U1",
            "constraint_id": "C3", "hp_family": "1", "phase_set": "100",
            "positions1": "10,30", "call_codes": "RR", "molecule_weight": 2,
            "pattern_count": 1, "disposition": "retained", "span_sites": 3,
            "crossed_cut_count": 0, "retained_block_index": 1,
        },
        {
            "dataset": "HCC1395", "chrom": "chr22", "unit_id": "U2",
            "constraint_id": "C4", "hp_family": "1", "phase_set": "200",
            "positions1": "10,20,30", "call_codes": "ARR", "molecule_weight": 3,
            "pattern_count": 1, "disposition": "retained", "span_sites": 3,
            "crossed_cut_count": 0, "retained_block_index": 1,
        },
    ]
    write_tsv_gz(partition / "blocks.tsv.gz", BLOCK_FIELDS, blocks)
    write_tsv_gz(
        partition / "dispositions.tsv.gz", DISPOSITION_FIELDS, dispositions
    )
    molecule_rows = []
    specifications = (
        ("m_aax", "100", "10,20", "AA", 3),
        ("m_raa", "100", "10,20,30", "RAA", 4),
        ("m_rr", "100", "10,30", "RR", 2),
        ("m_arr", "200", "10,20,30", "ARR", 3),
    )
    for prefix, phase_set, positions, codes, count in specifications:
        for index in range(count):
            molecule_rows.append(
                {
                    "dataset": "HCC1395",
                    "chrom": "chr22",
                    "molecule_id": f"{prefix}_{index}",
                    "hp_family": "1",
                    "phase_set": phase_set,
                    "positions1": positions,
                    "call_codes": codes,
                }
            )
    write_tsv_gz(
        chrom_dir / "extraction" / "HCC1395.chr22.molecule_sparse_calls.tsv.gz",
        MOLECULE_FIELDS,
        molecule_rows,
    )
    (root / "run_receipt.json").write_text(
        json.dumps(
            {
                "all_pass": True,
                "sample": "HCC1395",
                "metrics": {
                    "S": 3,
                    "unique_sites": 3,
                    "units": 2,
                    "unit_memberships": 6,
                    "blocks": 2,
                    "python_cpp_mismatch_count": 0,
                },
            }
        ),
        encoding="utf-8",
    )
    return root


def k13_with_singleton_child_fixture(root: Path) -> Path:
    """Model one k=13 source component partitioned into k=12 and k=1 children."""
    chrom_dir = root / "chromosomes" / "chr22"
    partition = chrom_dir / "python_partition"
    partition.mkdir(parents=True)
    (chrom_dir / "comparison.json").write_text(
        json.dumps({"all_pass": True, "mismatch_count": 0}), encoding="utf-8"
    )
    (partition / "receipt.json").write_text(
        json.dumps(
            {
                "all_pass": True,
                "scope": {"dataset": "HCC1395", "chrom": "chr22"},
                "checks": {
                    "cross_ps_zero": True,
                    "cross_hp_zero": True,
                    "constraint_mass_conserved": True,
                },
            }
        ),
        encoding="utf-8",
    )

    source_positions = tuple(range(100, 230, 10))
    retained_positions = source_positions[:12]
    singleton_position = source_positions[12]
    blocks = [
        {
            "dataset": "HCC1395", "chrom": "chr22", "unit_id": "U_BIG",
            "component_id": "CP_BIG", "linkage_basis": "PS_HP1", "hp_family": "1",
            "phase_set": "300", "threshold": 3, "block_id": "U_BIG:B0001",
            "block_index": 1, "start1": retained_positions[0],
            "end1": retained_positions[-1], "k": 12,
            "positions1": ",".join(map(str, retained_positions)),
            "retained_molecule_weight": 3, "retained_pattern_count": 1,
        },
        {
            "dataset": "HCC1395", "chrom": "chr22", "unit_id": "U_BIG",
            "component_id": "CP_BIG", "linkage_basis": "PS_HP1", "hp_family": "1",
            "phase_set": "300", "threshold": 3, "block_id": "U_BIG:B0002",
            "block_index": 2, "start1": singleton_position, "end1": singleton_position,
            "k": 1, "positions1": str(singleton_position),
            "retained_molecule_weight": 0, "retained_pattern_count": 0,
        },
    ]
    write_tsv_gz(partition / "blocks.tsv.gz", BLOCK_FIELDS, blocks)
    write_tsv_gz(partition / "dispositions.tsv.gz", DISPOSITION_FIELDS, [])

    molecule_rows = [
        {
            "dataset": "HCC1395",
            "chrom": "chr22",
            "molecule_id": f"m_big_{index}",
            "hp_family": "1",
            "phase_set": "300",
            "positions1": ",".join(map(str, source_positions)),
            "call_codes": "A" * len(source_positions),
        }
        for index in range(3)
    ]
    write_tsv_gz(
        chrom_dir / "extraction" / "HCC1395.chr22.molecule_sparse_calls.tsv.gz",
        MOLECULE_FIELDS,
        molecule_rows,
    )
    (root / "run_receipt.json").write_text(
        json.dumps(
            {
                "all_pass": True,
                "sample": "HCC1395",
                "metrics": {
                    "S": 13,
                    "unique_sites": 13,
                    "units": 1,
                    "unit_memberships": 13,
                    "blocks": 2,
                    "python_cpp_mismatch_count": 0,
                },
            }
        ),
        encoding="utf-8",
    )
    return root


def run_adapter(root: Path, output: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(ADAPTER),
            "--partition-root",
            str(root),
            "--output",
            str(output),
            "--sample",
            "HCC1395",
            "--chroms",
            "chr22",
            "--min-read",
            "3",
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def test_adapter_keeps_same_hp_different_ps_separate(tmp_path: Path) -> None:
    root = fixture(tmp_path / "input")
    output = tmp_path / "exact_mlhp.json"
    completed = run_adapter(root, output)
    assert completed.returncode == 0, completed.stderr
    document = json.loads(output.read_text(encoding="utf-8"))
    assert document["schema_name"] == "intersubmod.exact_ps_layered_tree_input"
    assert len(document["groups"]) == 2
    assert {group["phase_set"] for group in document["groups"]} == {"100", "200"}
    assert len({group["region_id"] for group in document["groups"]}) == 2
    ps100 = next(group for group in document["groups"] if group["phase_set"] == "100")
    assert ps100["populations_by_hp"] == {"1": {"RAA": 4}}
    assert ps100["subread_groups_by_hp"] == {"1": {"AAX": 3}}
    assert ps100["tree_supported_molecule_block_incidences"] == 7
    assert ps100["projected_molecule_block_incidences"] == 9
    assert ps100["retained_segmentation_constraint_weight"] == 9
    assert ps100["vaf_eligible"] is False
    assert document["input_funnel"]["check_constraint_weight_conservation"] is True
    assert document["input_funnel"]["constraint_weight_total"] == 12
    assert document["input_funnel"]["projected_molecule_block_incidences"] == 12


def test_adapter_rejects_constraint_with_different_ps(tmp_path: Path) -> None:
    root = fixture(tmp_path / "input", wrong_ps=True)
    completed = run_adapter(root, tmp_path / "bad.json")
    assert completed.returncode != 0
    assert "crosses HP/PS" in completed.stderr


def test_adapter_excludes_k1_child_of_multisite_source_component(
    tmp_path: Path,
) -> None:
    root = k13_with_singleton_child_fixture(tmp_path / "input")
    output = tmp_path / "exact_mlhp.json"
    completed = run_adapter(root, output)
    assert completed.returncode == 0, completed.stderr

    document = json.loads(output.read_text(encoding="utf-8"))
    assert document["input_funnel"]["bounded_blocks"] == 2
    assert document["input_funnel"]["k1_blocks_not_tree_eligible"] == 1
    assert document["input_funnel"]["tree_input_groups"] == 1
    assert document["input_funnel"]["tree_input_memberships"] == 12
    assert document["input_funnel"]["projected_molecule_block_incidences"] == 6
    assert document["input_funnel"]["tree_supported_molecule_block_incidences"] == 3
    assert [group["block_id"] for group in document["groups"]] == ["U_BIG:B0001"]
    assert document["groups"][0]["n_sSNV"] == 12
    assert all(group["n_sSNV"] >= 2 for group in document["groups"])

    receipt = json.loads(
        output.with_suffix(".json.receipt.json").read_text(encoding="utf-8")
    )
    assert receipt["all_pass"] is True
    assert receipt["counts"]["blocks_total"] == 2
    assert receipt["counts"]["blocks_k_1"] == 1
    assert receipt["counts"]["blocks_k1_not_tree_eligible"] == 1
    assert receipt["counts"]["groups_tree_input"] == 1


def test_exact_adapter_runs_layered_and_region_view_without_identity_collision(
    tmp_path: Path,
) -> None:
    root = fixture(tmp_path / "input")
    mlhp = tmp_path / "exact_mlhp.json"
    assert run_adapter(root, mlhp).returncode == 0
    layered = tmp_path / "layered.json"
    env = {
        **dict(os.environ),
        "PYTHONDONTWRITEBYTECODE": "1",
        "SM_ML": str(mlhp),
        "SM_OUT": str(layered),
        "SM_VERIFY_EVERY": "1",
        "SM_ANALYSIS_TREE_CAP": "0",
        "SM_DISPLAY_TREE_CAP": "0",
    }
    completed = subprocess.run(
        [sys.executable, str(LAYERED)], env=env, text=True, capture_output=True, check=False
    )
    assert completed.returncode == 0, completed.stderr
    layered_doc = json.loads(layered.read_text(encoding="utf-8"))
    assert layered_doc["schema_version"] == "2.1"
    assert layered_doc["sample"] == "HCC1395"
    assert (
        layered_doc["analysis_contract"]["PS"]
        == "exact non-missing primary evidence boundary; cross-PS pooling forbidden"
    )
    assert "no cross-PS join" in layered_doc["analysis_contract"]["region_rule"]
    assert layered_doc["read_tag_census"]["check_partition_receipts_all_pass"] is True
    assert layered_doc["read_tag_census"]["check_exact_sidecar_join"] is None
    assert len(layered_doc["detail"]) == 2
    assert {unit["phase_set"] for unit in layered_doc["detail"]} == {"100", "200"}

    region_view = tmp_path / "region_view.json"
    env.update(
        {
            "SM_LAYERED": str(layered),
            "SM_OUT": str(region_view),
            "SM_SAMPLE": "HCC1395",
            "SM_ML_GLOB": str(mlhp),
        }
    )
    completed = subprocess.run(
        [sys.executable, str(REGION_VIEW)],
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    region_doc = json.loads(region_view.read_text(encoding="utf-8"))
    assert len(region_doc["regions"]) == 2
    assert len({region["region"] for region in region_doc["regions"]}) == 2
    assert {region["phase_set"] for region in region_doc["regions"]} == {"100", "200"}
