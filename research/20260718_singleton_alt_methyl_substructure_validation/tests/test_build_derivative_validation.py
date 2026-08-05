from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_derivative_validation.py"
SPEC = importlib.util.spec_from_file_location("singleton_sidecar_builder", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)

PACKAGE_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "package_portable_report.py"
PACKAGE_SPEC = importlib.util.spec_from_file_location("singleton_sidecar_packager", PACKAGE_SCRIPT)
assert PACKAGE_SPEC and PACKAGE_SPEC.loader
PACKAGE_MODULE = importlib.util.module_from_spec(PACKAGE_SPEC)
sys.modules[PACKAGE_SPEC.name] = PACKAGE_MODULE
PACKAGE_SPEC.loader.exec_module(PACKAGE_MODULE)

COMPAT_DELIVERY_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "deliver_portable_artifact_scrollbar_compat.mjs"
)
QA_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "qa_portable_report.py"


def synthetic_sites() -> pd.DataFrame:
    return pd.DataFrame(
        [
            ("A", "chr1", 100, "A", "C"),
            ("A", "chr1", 50_100, "G", "T"),
            ("A", "chr1", 100_100, "C", "G"),
            ("A", "chr1", 150_101, "T", "A"),
            ("A", "chr2", 10, "A", "G"),
            ("B", "chr1", 100, "A", "T"),
        ],
        columns=["dataset", "chrom", "pos", "ref", "alt"],
    )


def test_positional_singleton_uses_transitive_gap_and_boundary() -> None:
    observed = MODULE.recompute_positional_singleton_keys(synthetic_sites())
    assert observed == {
        ("A", "chr1", 150_101, "T", "A"),
        ("A", "chr2", 10, "A", "G"),
        ("B", "chr1", 100, "A", "T"),
    }


def test_nullable_bool_parser_preserves_expected_missing_values() -> None:
    observed = MODULE.parse_nullable_bool_series(
        pd.Series(["true", "FALSE", None, np.nan]),
        "nullable",
    )
    assert observed.tolist()[:2] == [True, False]
    assert observed.isna().tolist() == [False, False, True, True]


def test_nullable_bool_parser_rejects_unknown_values() -> None:
    with pytest.raises(MODULE.ValidationError, match="unexpected boolean values"):
        MODULE.parse_nullable_bool_series(pd.Series(["true", "unknown"]), "nullable")


@pytest.mark.parametrize("value", ["chr1-22", "chr11-12", "range 1-20", "v1-2foo"])
def test_source_label_guard_allows_non_label_substrings(value: str) -> None:
    assert not MODULE.contains_source_cluster_label(value)


@pytest.mark.parametrize("value", ["1-1", "label=1-2", "source group 1-1."])
def test_source_label_guard_rejects_standalone_source_labels(value: str) -> None:
    assert MODULE.contains_source_cluster_label(value)


def test_sqlite_snapshot_roundtrip_executes_recorded_query() -> None:
    rows = [
        {"label": "first", "count": 3, "ratio": 0.5, "flag": True, "missing": None},
        {"label": "second", "count": 0, "ratio": 0.0, "flag": False, "missing": None},
    ]
    connection = MODULE.sqlite3.connect(":memory:")
    try:
        observed, sql = MODULE.roundtrip_snapshot_dataset(connection, "focused_test", rows)
        rerun = pd.read_sql_query(sql, connection)
    finally:
        connection.close()
    assert observed == rows
    assert len(rerun) == 2
    assert sql == (
        'SELECT "label", "count", "ratio", "flag", "missing" '
        'FROM "artifact_focused_test" ORDER BY "_row_order" ASC'
    )


def test_portable_packager_records_explicit_bounded_timeouts() -> None:
    args = PACKAGE_MODULE.parse_args(
        [
            "--artifact",
            "artifact.json",
            "--html",
            "report.html",
            "--receipt",
            "receipt.json",
            "--ready-timeout-ms",
            "15000",
            "--action-timeout-ms",
            "5000",
            "--timeout-ms",
            "30000",
        ]
    )
    assert args.ready_timeout_ms == 15_000
    assert args.action_timeout_ms == 5_000
    assert args.timeout_ms == 30_000


def test_compat_delivery_reuses_canonical_plugin_runtime() -> None:
    source = COMPAT_DELIVERY_SCRIPT.read_text(encoding="utf-8")
    assert "buildPortableArtifact" in source
    assert "extractPortableChartSvgs" in source
    assert "verifyPortableArtifact" in source
    assert "analytics-top-bar" in source
    assert "<html" not in source.lower()


def test_browser_qa_counts_only_enhanced_reader_visual_roots() -> None:
    source = QA_SCRIPT.read_text(encoding="utf-8")
    assert '#data-analytics-portable-reader section[data-artifact-kind="chart"]' in source
    assert '#data-analytics-portable-reader section[data-artifact-kind="table"]' in source
    assert '"report-html-frame" not in owner_class.split()' in source


def test_summarize_group_keeps_all_m2_states_disjoint() -> None:
    frame = pd.DataFrame(
        {
            "_m1_evaluable": [True, True, True, False],
            "_stable": [True, True, False, False],
            "m2_status": ["PASS", "NOT_EVALUABLE", "NOT_RUN", "NOT_RUN"],
        }
    )
    summary = MODULE.summarize_group(frame)
    assert summary["sites"] == 4
    assert summary["m1_evaluable"] == 3
    assert summary["m1_flagged"] == 2
    assert summary["m2_pass"] == 1
    assert summary["m2_not_evaluable"] == 1
    assert summary["m2_not_run"] == 2
    assert summary["m2_fail"] == 0
    assert summary["m2_determinate"] == 1
    assert summary["m2_pass_pct_all"] == 0.25
    assert summary["m2_pass_pct_determinate"] == 1.0


def test_core_distance_gate_requires_finite_symmetric_zero_diagonal() -> None:
    valid = pd.DataFrame(
        {
            "read_id": ["0", "1", "2"],
            "0": [0.0, 0.2, np.nan],
            "1": [0.2, 0.0, 0.3],
            "2": [np.nan, 0.3, 0.0],
        }
    )
    core, checks = MODULE.validate_square_distance(valid, ["0", "1"])
    assert core.shape == (2, 2)
    assert checks["core_finite"] is True
    assert checks["raw_matrix_nonfinite_cells"] == 2

    asymmetric = valid.copy()
    asymmetric.loc[1, "0"] = 0.4
    with pytest.raises(MODULE.ValidationError, match="asymmetric"):
        MODULE.validate_square_distance(asymmetric, ["0", "1"])

    nonfinite_core = valid.copy()
    nonfinite_core.loc[0, "1"] = np.nan
    with pytest.raises(MODULE.ValidationError, match="non-finite"):
        MODULE.validate_square_distance(nonfinite_core, ["0", "1"])

    nonzero_diagonal = valid.copy()
    nonzero_diagonal.loc[0, "0"] = 0.1
    with pytest.raises(MODULE.ValidationError, match="diagonal"):
        MODULE.validate_square_distance(nonzero_diagonal, ["0", "1"])


def test_methylation_gate_checks_read_and_cpg_identity() -> None:
    reads = pd.DataFrame({"read_id": ["0", "1"], "read_name": ["a", "b"]})
    methyl = pd.DataFrame({"read_id": ["0", "1"], "100": [0.1, np.nan], "200": [0.9, 0.2]})
    cpg = pd.DataFrame({"position": ["100", "200"]})
    core, checks = MODULE.validate_methylation(methyl, reads, cpg, ["0", "1"])
    assert core.shape == (2, 2)
    assert checks["core_called_cells"] == 3
    assert checks["core_missing_cells"] == 1

    wrong_cpg = pd.DataFrame({"position": ["100", "201"]})
    with pytest.raises(MODULE.ValidationError, match="column identity"):
        MODULE.validate_methylation(methyl, reads, wrong_cpg, ["0", "1"])


def test_medoid_selection_is_deterministic_and_capped() -> None:
    core_ids = [str(value) for value in range(36)]
    group_by_id = {
        read_id: ("Group A" if int(read_id) < 18 else "Group B")
        for read_id in core_ids
    }
    coordinates = np.array([int(value) % 18 for value in core_ids], dtype=float)
    distance = np.abs(coordinates[:, None] - coordinates[None, :])
    selected_first, audit_first = MODULE.select_medoid_nearest(core_ids, group_by_id, distance)
    selected_second, audit_second = MODULE.select_medoid_nearest(core_ids, group_by_id, distance)
    assert selected_first == selected_second
    assert audit_first == audit_second
    assert len(selected_first) == 30
    assert sum(group_by_id[value] == "Group A" for value in selected_first) == 15
    assert sum(group_by_id[value] == "Group B" for value in selected_first) == 15
    assert {row["group"] for row in audit_first} == {"Group A", "Group B"}


def test_medoid_selection_rejects_hp_like_group_names() -> None:
    with pytest.raises(MODULE.ValidationError, match="Group A/B"):
        MODULE.select_medoid_nearest(
            ["0", "1"],
            {"0": "1-1", "1": "1-2"},
            np.zeros((2, 2)),
        )


def test_no_clobber_output_directory_contract() -> None:
    output = MagicMock(spec=Path)
    output.resolve.return_value = MODULE.TOPIC_ROOT / "results" / "existing"
    output.parent.mkdir = MagicMock()
    output.mkdir.side_effect = FileExistsError("existing")
    with pytest.raises(FileExistsError):
        MODULE.ensure_output_dir(output)


def test_artifact_source_is_reader_safe_and_claim_guardrails_are_literal() -> None:
    source = Path(MODULE.SCRIPT if hasattr(MODULE, "SCRIPT") else SCRIPT).read_text(encoding="utf-8")
    for phrase in (
        "observed operational M2-PASS yield",
        "不是 biological/subclone prevalence",
        "PS 不在 M2 measured axes",
        "confirmed cellular subclone=0",
        "linear ancestry=0",
        "不是真陰性",
        "matched-normal",
        "tumor-REF",
    ):
        assert phrase in source
    assert '"type": "html"' in source
    assert "data:image/png;base64" in source
    assert "Group A/B" in source


def test_expected_fixed_denominators_are_locked() -> None:
    assert MODULE.EXPECTED_COUNTS == {
        "all_dataset_sites": 469_849,
        "singleton_sites": 50_432,
        "m1_evaluable": 48_347,
        "m1_flagged": 5_961,
        "m2_not_run": 44_471,
        "m2_not_evaluable": 5_913,
        "m2_fail": 18,
        "m2_pass": 30,
    }
    assert MODULE.EXPECTED_HCC1395["sites"] == 8_279
    assert MODULE.EXPECTED_HCC1395["m1_flagged"] == 734
    assert MODULE.EXPECTED_HCC1395["m2_pass"] == 2
