import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    REPO_ROOT
    / "research/20260718_k_gt8_read_supported_segmentation"
    / "scripts/build_k_gt8_partitions_span_cap.py"
)
K_ONLY_SCRIPT = (
    REPO_ROOT
    / "research/20260718_k_gt8_read_supported_segmentation"
    / "scripts/build_k_gt8_partitions.py"
)
DATASET = "SYNTHETIC_SPAN"
CHROM = "chr1"
SITE_FIELDS = ("site_index", "chrom", "pos1", "ref", "alt")
CALL_FIELDS = (
    "dataset",
    "chrom",
    "molecule_id",
    "hp_family",
    "phase_set",
    "site_indices",
    "call_codes",
)


def _write_tsv_gz(path, fieldnames, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv_gz(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_cli_applies_physical_span_cap_and_audits_outputs(tmp_path):
    positions = (100, 200, 300, 400, 10_000, 10_100, 10_200, 10_300, 10_400)
    site_catalog = tmp_path / "sites.tsv.gz"
    molecule_calls = tmp_path / "calls.tsv.gz"
    output_dir = tmp_path / "output"
    _write_tsv_gz(
        site_catalog,
        SITE_FIELDS,
        (
            {
                "site_index": index,
                "chrom": CHROM,
                "pos1": position,
                "ref": "C",
                "alt": "T",
            }
            for index, position in enumerate(positions)
        ),
    )
    _write_tsv_gz(
        molecule_calls,
        CALL_FIELDS,
        (
            {
                "dataset": DATASET,
                "chrom": CHROM,
                "molecule_id": molecule_id,
                "hp_family": "1",
                "phase_set": "100",
                "site_indices": ",".join(map(str, indices)),
                "call_codes": "A" * len(indices),
            }
            for molecule_id, indices in (
                ("left", (0, 1, 2, 3)),
                ("bridge", (3, 4)),
                ("right", (4, 5, 6, 7, 8)),
            )
        ),
    )

    command = [
        sys.executable,
        str(SCRIPT),
        "--dataset",
        DATASET,
        "--chrom",
        CHROM,
        "--site-catalog",
        str(site_catalog),
        "--molecule-calls",
        str(molecule_calls),
        "--output-dir",
        str(output_dir),
        "--max-block-span-bp",
        "500",
        "--expected-target-components",
        "1",
        "--expected-target-sites",
        "9",
    ]
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 0, (
        f"command failed: {' '.join(command)}\n"
        f"stdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}"
    )

    receipt = json.loads((output_dir / "receipt.json").read_text())
    assert receipt["schema_name"].endswith("_span_cap")
    assert receipt["parameters"]["max_block_size"] == 8
    assert receipt["parameters"]["max_block_span_bp"] == 500
    assert receipt["all_pass"] is True
    assert all(receipt["checks"].values())
    assert receipt["checks"]["every_block_span_lte_max"] is True
    assert receipt["counts"]["target_components"] == 1
    assert receipt["counts"]["target_sites"] == 9
    assert receipt["counts"]["new_blocks"] == 2
    assert receipt["counts"]["raw_total_molecule_weight"] == 3
    assert receipt["counts"]["raw_retained_molecule_weight"] == 2
    assert receipt["counts"]["raw_lost_molecule_weight"] == 1
    assert receipt["counts"]["unavoidable_patterns"] == 1
    assert receipt["counts"]["unavoidable_size_patterns"] == 0
    assert receipt["counts"]["unavoidable_span_cap_patterns"] == 1
    assert receipt["counts"]["unavoidable_both_limits_patterns"] == 0

    blocks = _read_tsv_gz(
        output_dir / f"{DATASET}.{CHROM}.blocks.tsv.gz"
    )
    assert [int(row["k"]) for row in blocks] == [4, 5]
    assert [int(row["span_bp"]) for row in blocks] == [300, 400]
    assert all(int(row["span_bp"]) <= 500 for row in blocks)

    constraints = _read_tsv_gz(
        output_dir / f"{DATASET}.{CHROM}.cut_constraints.tsv.gz"
    )
    bridge = next(row for row in constraints if row["positions"] == "400,10000")
    assert bridge["span_sites"] == "2"
    assert bridge["span_bp"] == "9600"
    assert (
        bridge["disposition"]
        == "unavoidable_bp_span_gt_max_block_span_bp"
    )
    assert bridge["crossed_cut_count"] == "1"
    assert bridge["retained_block_index"] == ""


def test_cli_without_span_cap_matches_current_k_only_semantics(tmp_path):
    positions = tuple(range(100, 1_000, 100))
    site_catalog = tmp_path / "parity.sites.tsv.gz"
    molecule_calls = tmp_path / "parity.calls.tsv.gz"
    _write_tsv_gz(
        site_catalog,
        SITE_FIELDS,
        (
            {
                "site_index": index,
                "chrom": CHROM,
                "pos1": position,
                "ref": "C",
                "alt": "T",
            }
            for index, position in enumerate(positions)
        ),
    )
    _write_tsv_gz(
        molecule_calls,
        CALL_FIELDS,
        (
            {
                "dataset": DATASET,
                "chrom": CHROM,
                "molecule_id": molecule_id,
                "hp_family": "1",
                "phase_set": "100",
                "site_indices": ",".join(map(str, indices)),
                "call_codes": "A" * len(indices),
            }
            for molecule_id, indices in (
                ("left", (0, 1, 2, 3, 4)),
                ("right", (4, 5, 6, 7, 8)),
            )
        ),
    )

    receipts = []
    output_dirs = []
    for label, script in (
        ("k_only", K_ONLY_SCRIPT),
        ("span_none", SCRIPT),
    ):
        output_dir = tmp_path / label
        command = [
            sys.executable,
            str(script),
            "--dataset",
            DATASET,
            "--chrom",
            CHROM,
            "--site-catalog",
            str(site_catalog),
            "--molecule-calls",
            str(molecule_calls),
            "--output-dir",
            str(output_dir),
            "--expected-target-components",
            "1",
            "--expected-target-sites",
            "9",
        ]
        completed = subprocess.run(
            command,
            cwd=REPO_ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert completed.returncode == 0, (
            f"command failed: {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
        receipts.append(
            json.loads((output_dir / "receipt.json").read_text())
        )
        output_dirs.append(output_dir)

    old_receipt, new_receipt = receipts
    old_dir, new_dir = output_dirs
    assert new_receipt["parameters"]["max_block_span_bp"] is None
    assert (
        new_receipt["semantic_result_sha256"]
        == old_receipt["semantic_result_sha256"]
    )
    assert new_receipt["status_counts"] == old_receipt["status_counts"]
    for key, value in old_receipt["counts"].items():
        assert new_receipt["counts"][key] == value

    old_components = _read_tsv_gz(
        old_dir / f"{DATASET}.{CHROM}.legacy_components.tsv.gz"
    )
    new_components = _read_tsv_gz(
        new_dir / f"{DATASET}.{CHROM}.legacy_components.tsv.gz"
    )
    parity_fields = (
        "legacy_component_id",
        "pre_cap_k",
        "new_block_count",
        "raw_cut_indices",
        "equal_pattern_cut_indices",
        "log1p_cut_indices",
        "raw_retained_molecule_weight",
        "raw_lost_molecule_weight",
        "status",
    )
    assert [
        tuple(row[field] for field in parity_fields)
        for row in new_components
    ] == [
        tuple(row[field] for field in parity_fields)
        for row in old_components
    ]
