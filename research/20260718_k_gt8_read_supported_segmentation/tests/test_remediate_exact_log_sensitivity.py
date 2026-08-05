import csv
import gzip
import hashlib
import importlib.util
import itertools
import json
import random
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT = (
    REPO_ROOT
    / "research/20260718_k_gt8_read_supported_segmentation"
    / "scripts/remediate_exact_log_sensitivity.py"
)
SPEC = importlib.util.spec_from_file_location("exact_log_remediation", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)

Constraint = MODULE.Constraint
evaluate_cuts = MODULE.evaluate_cuts
semantic_sha256 = MODULE.semantic_sha256
solve_ordered_partition = MODULE.solve_ordered_partition


def _best_brute_force(positions, constraints, max_block_size, mode):
    n_sites = len(positions)
    candidates = []
    for n_cuts in range(n_sites):
        for cuts in itertools.combinations(range(1, n_sites), n_cuts):
            bounds = (0,) + cuts + (n_sites,)
            if any(
                right - left > max_block_size
                for left, right in zip(bounds, bounds[1:])
            ):
                continue
            values = evaluate_cuts(
                positions,
                constraints,
                cuts,
                max_block_size=max_block_size,
            )
            objective = {
                "raw": values["raw"],
                "equal": values["patterns"],
                "exact_product": values["product"],
            }[mode]
            candidates.append(
                (
                    objective,
                    values["patterns"],
                    -values["block_count"],
                    values["gap_sum"],
                    tuple(-cut for cut in cuts),
                    cuts,
                )
            )
    assert candidates
    return max(candidates)[-1]


def test_red_team_equal_raw_pattern_product_uses_later_tie_break():
    # Two mutually exclusive retained pattern sets:
    # left weights 1+5+5 = 11, factors 2*6*6 = 72;
    # right weights 2+2+7 = 11, factors 3*3*8 = 72.
    # With k=9 and K=5, cut 4 retains the first set and cut 5 retains
    # the second.  Every primary and exact-product criterion ties, so the
    # declared final lexicographic tie-break must select cut 4.
    positions = tuple(range(1, 10))
    constraints = (
        Constraint("left_1", 4, 8, 1),
        Constraint("left_5a", 4, 8, 5),
        Constraint("left_5b", 4, 8, 5),
        Constraint("right_2a", 0, 4, 2),
        Constraint("right_2b", 0, 4, 2),
        Constraint("right_7", 0, 4, 7),
    )

    raw = solve_ordered_partition(
        positions, constraints, max_block_size=5, mode="raw"
    )
    equal = solve_ordered_partition(
        positions, constraints, max_block_size=5, mode="equal"
    )
    exact = solve_ordered_partition(
        positions, constraints, max_block_size=5, mode="exact_product"
    )

    assert raw.cuts == equal.cuts == exact.cuts == (4,)
    assert raw.objective == 11
    assert exact.retained_patterns == 3
    assert exact.objective == 72
    assert evaluate_cuts(
        positions, constraints, (5,), max_block_size=5
    ) == {
        "raw": 11,
        "patterns": 3,
        "product": 72,
        "block_count": 2,
        "gap_sum": 1,
    }


@pytest.mark.parametrize("mode", ("raw", "equal", "exact_product"))
def test_dynamic_program_matches_brute_force_oracle(mode):
    random_generator = random.Random(20260718)
    for trial in range(250):
        n_sites = random_generator.randint(2, 10)
        max_block_size = random_generator.randint(1, min(5, n_sites))
        gaps = [random_generator.randint(1, 30) for _ in range(n_sites)]
        positions = []
        cursor = 0
        for gap in gaps:
            cursor += gap
            positions.append(cursor)
        constraints = []
        for ordinal in range(random_generator.randint(0, 14)):
            lo = random_generator.randrange(n_sites)
            hi = random_generator.randrange(lo, n_sites)
            constraints.append(
                Constraint(
                    f"T{trial:03d}_{ordinal:03d}",
                    lo,
                    hi,
                    random_generator.randint(1, 30),
                )
            )
        result = solve_ordered_partition(
            positions,
            constraints,
            max_block_size=max_block_size,
            mode=mode,
        )
        expected = _best_brute_force(
            positions, constraints, max_block_size, mode
        )
        assert result.cuts == expected


def _write_tsv_gz(path, fields, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _identity(path):
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }


def _write_sidecar(path):
    Path(f"{path}.sha256").write_text(
        f"{_sha256(path)}  {path.name}\n", encoding="utf-8"
    )


def _synthetic_partition(tmp_path, legacy_log_cut="8"):
    source = tmp_path / "source"
    source.mkdir()
    prefix = "SYNTHETIC.chr1"
    components = source / f"{prefix}.legacy_components.tsv.gz"
    constraints = source / f"{prefix}.cut_constraints.tsv.gz"
    membership = source / f"{prefix}.site_membership.tsv.gz"
    receipt = source / "receipt.json"

    # K=8 production analogue of the 9-site/K=5 red-team case: only cuts
    # 7 and 8 can retain either complete 8-site interval, and no partition
    # can retain both mutually exclusive intervals.
    positions = tuple(range(1, 16))
    component_id = "SYNTHETIC:chr1:legacy_gap_50000:1:1-15"
    component_row = {
        "dataset": "SYNTHETIC",
        "chrom": "chr1",
        "legacy_component_id": component_id,
        "start1": "1",
        "end1": "15",
        "pre_cap_k": "15",
        "exact_pattern_count": "6",
        "raw_total_molecule_weight": "22",
        "raw_retained_molecule_weight": "11",
        "raw_lost_molecule_weight": "11",
        "retained_exact_pattern_count": "3",
        "lost_exact_pattern_count": "3",
        "weight_stable": "false" if legacy_log_cut != "7" else "true",
        "raw_cut_indices": "7",
        "equal_pattern_cut_indices": "7",
        "log1p_cut_indices": legacy_log_cut,
        "positions_sha256": semantic_sha256(positions),
    }
    _write_tsv_gz(components, tuple(component_row), [component_row])

    membership_rows = [
        {
            "dataset": "SYNTHETIC",
            "chrom": "chr1",
            "legacy_component_id": component_id,
            "site_index": index,
            "component_local_index": index,
            "pos1": position,
        }
        for index, position in enumerate(positions)
    ]
    _write_tsv_gz(
        membership, tuple(membership_rows[0]), membership_rows
    )

    pattern_specs = (
        ("left_1", (8, 15), 1),
        ("left_5a", (8, 15), 5),
        ("left_5b", (8, 15), 5),
        ("right_2a", (1, 8), 2),
        ("right_2b", (1, 8), 2),
        ("right_7", (1, 8), 7),
    )
    constraint_rows = [
        {
            "dataset": "SYNTHETIC",
            "chrom": "chr1",
            "legacy_component_id": component_id,
            "constraint_id": constraint_id,
            "positions": ",".join(map(str, pattern_positions)),
            "n_fixed_ra": "2",
            "span_sites": "8",
            "molecule_weight": str(weight),
        }
        for constraint_id, pattern_positions, weight in pattern_specs
    ]
    _write_tsv_gz(
        constraints, tuple(constraint_rows[0]), constraint_rows
    )

    receipt_value = {
        "schema_name": "intersubmod.k_gt8_read_supported_segmentation",
        "schema_version": "0.1.0",
        "all_pass": True,
        "parameters": {
            "max_block_size": 8,
            "tie_break": [
                "retained_weight_desc",
                "retained_pattern_count_desc",
                "block_count_asc",
                "cut_gap_sum_desc",
                "cut_indices_lexicographic_asc",
            ],
        },
        "outputs": {
            "legacy_components": _identity(components),
            "cut_constraints": _identity(constraints),
            "site_membership": _identity(membership),
        },
    }
    receipt.write_text(
        json.dumps(receipt_value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_sidecar(receipt)
    return components, constraints, membership, receipt


def test_cli_is_fail_closed_and_remediates_only_log_sensitivity(tmp_path):
    components, constraints, _, _ = _synthetic_partition(tmp_path)
    output = tmp_path / "output"
    command = [
        sys.executable,
        str(SCRIPT),
        "--components",
        str(components),
        "--cut-constraints",
        str(constraints),
        "--output-dir",
        str(output),
        "--expected-components",
        "1",
    ]
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert completed.returncode == 0, completed.stderr
    receipt = json.loads((output / "receipt.json").read_text())
    summary = json.loads((output / "summary.json").read_text())
    with gzip.open(
        output / "exact_log_sensitivity.tsv.gz",
        "rt",
        encoding="utf-8",
        newline="",
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert receipt["all_pass"] is True
    assert all(receipt["checks"].values())
    assert receipt["parameters"]["primary_raw_partition_changed"] is False
    assert summary["interpretation"]["primary_raw_partition_changed"] is False
    assert len(rows) == 1
    row = rows[0]
    assert row["raw_cut_indices"] == "7"
    assert row["equal_pattern_cut_indices"] == "7"
    assert row["legacy_log1p_cut_indices"] == "8"
    assert row["exact_product_cut_indices"] == "7"
    assert row["legacy_log_matches_exact"] == "false"
    assert row["legacy_log_exact_product_relation"] == "EQ"
    assert row["corrected_weight_stable"] == "true"
    assert row["correction_changed_stability"] == "true"
    assert row["remediation_class"] == "LEGACY_LOG_SUBOPTIMAL_LEX_TIEBREAK"

    second = subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert second.returncode == 2
    assert "output directory must not exist" in second.stderr


def test_cli_rejects_tampered_coordinate_witness_before_output(tmp_path):
    components, constraints, membership, _ = _synthetic_partition(tmp_path)
    with gzip.open(membership, "ab") as handle:
        handle.write(b"tamper")
    output = tmp_path / "must_not_exist"
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--components",
            str(components),
            "--cut-constraints",
            str(constraints),
            "--output-dir",
            str(output),
        ],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert completed.returncode == 2
    assert "source output identity mismatch" in completed.stderr
    assert not output.exists()


def test_full_mode_requires_authenticated_success_before_output(tmp_path):
    full_root = tmp_path / "unfinished_full_root"
    full_root.mkdir()
    output = tmp_path / "must_not_exist"
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--full-root",
            str(full_root),
            "--output-dir",
            str(output),
            "--expected-components",
            "408",
        ],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert completed.returncode == 2
    assert "_SUCCESS" in completed.stderr
    assert not output.exists()


REAL_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260718_k_gt8_read_supported_segmentation"
)


@pytest.mark.parametrize(
    ("chrom", "partition", "expected_components", "expected_sites"),
    (
        (
            "chr22",
            REAL_ROOT / "probes/HCC1395_chr22/partition_v2",
            6,
            98,
        ),
        (
            "chr6",
            REAL_ROOT
            / "full/HCC1395_chr1_22_v1/chromosomes/chr6/partition",
            83,
            25_657,
        ),
    ),
)
def test_real_hcc1395_probe(chrom, partition, expected_components, expected_sites):
    if not partition.is_dir():
        pytest.skip(f"real probe is unavailable: {partition}")
    prefix = f"HCC1395.{chrom}"
    rows, counts, source = MODULE.process_bundle(
        MODULE.SourceBundle(
            partition / f"{prefix}.legacy_components.tsv.gz",
            partition / f"{prefix}.cut_constraints.tsv.gz",
            partition / f"{prefix}.site_membership.tsv.gz",
            partition / "receipt.json",
        )
    )
    assert source["chrom"] == chrom
    assert len(rows) == counts["components"] == expected_components
    assert counts["sites"] == expected_sites
    assert all(row["raw_recomputed_matches_source"] == "true" for row in rows)
    assert all(row["equal_recomputed_matches_source"] == "true" for row in rows)
