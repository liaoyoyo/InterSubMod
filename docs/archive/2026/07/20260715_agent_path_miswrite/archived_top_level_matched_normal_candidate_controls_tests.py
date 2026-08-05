from __future__ import annotations

import importlib.util
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "analyze_matched_normal_candidate_controls.py"
SPEC = importlib.util.spec_from_file_location("matched_normal_candidate_controls", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def make_read(
    read_id: str,
    *,
    tumor: bool,
    allele: str,
    start: int,
    hp: str = "poison-normal-hp",
) -> dict[str, str]:
    return {
        "read_id": read_id,
        "read_name": f"read-{read_id}",
        "chr": "chr1",
        "start": str(start),
        "end": str(start + 100),
        "mapq": "60",
        "hp": hp,
        "alt_support": allele,
        "is_tumor": "1" if tumor else "0",
        "strand": "+" if start % 2 == 0 else "-",
    }


def primary_core(row: dict[str, str], label: str) -> dict[str, object]:
    return {
        "read_name": row["read_name"],
        "chrom": row["chr"],
        "start": int(row["start"]),
        "end": int(row["end"]),
        "mapq": int(row["mapq"]),
        "strand": row["strand"],
        "label": label,
    }


def candidate() -> dict[str, object]:
    return {"sample": "S1", "chrom": "chr1", "pos": 123, "ref": "A", "alt": "G"}


def matrices(reads: dict[str, dict[str, str]]) -> tuple[list[str], np.ndarray, np.ndarray]:
    identifiers = list(reads)
    return (
        identifiers,
        np.zeros((len(identifiers), len(identifiers)), dtype=float),
        np.zeros((len(identifiers), 4), dtype=float),
    )


def fake_subset_factory(calls: list[list[str]]):
    def fake_subset(
        ids: list[str],
        _distance_ids: list[str],
        _distance: np.ndarray,
        _methylation_ids: list[str],
        _methylation: np.ndarray,
        _seed: int,
    ) -> dict[str, object]:
        calls.append(list(ids))
        evaluable = len(ids) >= 2 * MODULE.F.MIN_SIZE
        return {
            "status": "evaluable" if evaluable else "insufficient_reads",
            "evaluable": evaluable,
            "n_raw": len(ids),
            "n_matrix": len(ids) if evaluable else 0,
            "n_after_peel": len(ids) if evaluable else 0,
            "stable_multigroup": bool(evaluable and ids and ids[0].startswith("tr")),
            "group_sizes": {"1": len(ids)} if evaluable else {},
        }

    return fake_subset


def test_normal_methyl_background_receives_only_focal_ref_reads_and_hp_is_unused(
    monkeypatch,
) -> None:
    reads: dict[str, dict[str, str]] = {}
    for index in range(10):
        reads[f"nr{index}"] = make_read(
            f"nr{index}", tumor=False, allele="REF", start=100 + index
        )
    reads["na0"] = make_read("na0", tumor=False, allele="ALT", start=120)
    reads["nu0"] = make_read("nu0", tumor=False, allele="UNKNOWN", start=121)
    reads["nu1"] = make_read("nu1", tumor=False, allele=".", start=122)
    for index in range(6):
        reads[f"tr{index}"] = make_read(
            f"tr{index}", tumor=True, allele="REF", start=200 + index
        )
    reads["ta0"] = make_read("ta0", tumor=True, allele="ALT", start=220)
    reads["ta1"] = make_read("ta1", tumor=True, allele="ALT", start=221)
    assignment = {
        "core_reads": [
            primary_core(reads["ta0"], "group-1"),
            primary_core(reads["ta1"], "group-2"),
        ]
    }
    identifiers, distance, methylation = matrices(reads)
    calls: list[list[str]] = []
    monkeypatch.setattr(MODULE, "analyze_subset", fake_subset_factory(calls))

    first = MODULE.analyze_site_from_data(
        candidate(), reads, identifiers, distance, identifiers, methylation, assignment
    )
    changed_hp = deepcopy(reads)
    for read_id in [read_id for read_id in changed_hp if read_id.startswith("n")]:
        changed_hp[read_id]["hp"] = f"changed-{read_id}"
    calls_after_change: list[list[str]] = []
    monkeypatch.setattr(MODULE, "analyze_subset", fake_subset_factory(calls_after_change))
    second = MODULE.analyze_site_from_data(
        candidate(), changed_hp, identifiers, distance, identifiers, methylation, assignment
    )

    assert first == second
    assert calls[0] == [f"nr{index}" for index in range(10)]
    assert calls[1] == [f"tr{index}" for index in range(6)]
    assert calls_after_change == calls
    assert first["normal_called_reads"] == 11
    assert first["normal_ref_reads"] == 10
    assert first["normal_alt_reads"] == 1
    assert first["normal_unknown_reads"] == 2
    assert first["normal_focal_callability_gate"] is True
    assert first["normal_control_evaluable"] is True
    assert first["normal_ref_methyl_n_raw"] == 10
    assert first["normal_genetic_alt_absence_gate"] is False
    assert first["normal_ref_methyl_nonreplication_gate"] is True
    assert first["normal_control_status"] == "REF_METHYL_PARTITION_NOT_REPRODUCED"
    assert first["normal_hp_used"] is False
    assert first["normal_hp_policy"] == "PROHIBITED_NOT_USED"


def test_all_unknown_normal_reads_never_pass_callability_or_evaluability(monkeypatch) -> None:
    reads = {
        f"nu{index}": make_read(
            f"nu{index}", tumor=False, allele="UNKNOWN", start=100 + index
        )
        for index in range(12)
    }
    reads["ta0"] = make_read("ta0", tumor=True, allele="ALT", start=300)
    assignment = {"core_reads": [primary_core(reads["ta0"], "group-1")]}
    identifiers, distance, methylation = matrices(reads)
    calls: list[list[str]] = []
    monkeypatch.setattr(MODULE, "analyze_subset", fake_subset_factory(calls))

    row = MODULE.analyze_site_from_data(
        candidate(), reads, identifiers, distance, identifiers, methylation, assignment
    )
    assert calls[0] == []
    assert row["normal_called_reads"] == 0
    assert row["normal_ref_reads"] == 0
    assert row["normal_alt_reads"] == 0
    assert row["normal_unknown_reads"] == 12
    assert row["normal_focal_callability_gate"] is False
    assert row["normal_control_evaluable"] is False
    assert row["normal_ref_methyl_nonreplication_gate"] is None
    assert row["normal_control_status"] == "NOT_EVALUABLE"
    assert row["normal_control_not_evaluable_reason"] == "NORMAL_CALLED_READS_LT_10"


def test_callability_does_not_convert_insufficient_ref_methyl_depth_into_negative(
    monkeypatch,
) -> None:
    reads: dict[str, dict[str, str]] = {}
    for index in range(5):
        reads[f"nr{index}"] = make_read(
            f"nr{index}", tumor=False, allele="REF", start=100 + index
        )
        reads[f"na{index}"] = make_read(
            f"na{index}", tumor=False, allele="ALT", start=120 + index
        )
    reads["ta0"] = make_read("ta0", tumor=True, allele="ALT", start=300)
    assignment = {"core_reads": [primary_core(reads["ta0"], "group-1")]}
    identifiers, distance, methylation = matrices(reads)
    monkeypatch.setattr(MODULE, "analyze_subset", fake_subset_factory([]))

    row = MODULE.analyze_site_from_data(
        candidate(), reads, identifiers, distance, identifiers, methylation, assignment
    )
    assert row["normal_called_reads"] == 10
    assert row["normal_ref_reads"] == 5
    assert row["normal_focal_callability_gate"] is True
    assert row["normal_ref_methyl_evaluable"] is False
    assert row["normal_control_evaluable"] is False
    assert row["normal_ref_methyl_nonreplication_gate"] is None
    assert row["normal_control_not_evaluable_reason"] == (
        "NORMAL_REF_METHYL_INSUFFICIENT_READS"
    )


def test_summary_keeps_not_evaluable_separate_from_negative(monkeypatch) -> None:
    reads = {
        f"nu{index}": make_read(
            f"nu{index}", tumor=False, allele="UNKNOWN", start=100 + index
        )
        for index in range(12)
    }
    reads["ta0"] = make_read("ta0", tumor=True, allele="ALT", start=300)
    assignment = {"core_reads": [primary_core(reads["ta0"], "group-1")]}
    identifiers, distance, methylation = matrices(reads)
    monkeypatch.setattr(MODULE, "analyze_subset", fake_subset_factory([]))
    row = MODULE.analyze_site_from_data(
        candidate(), reads, identifiers, distance, identifiers, methylation, assignment
    )
    summary = MODULE.summarize([row])
    assert summary["n_normal_control_evaluable"] == 0
    assert summary["n_normal_control_not_evaluable"] == 1
    assert summary["n_normal_ref_methyl_nonreplicating"] == 0


def test_normal_negative_guardrail_explicitly_preserves_asm_and_cell_state() -> None:
    assert "allele-specific methylation (ASM)" in MODULE.GUARDRAIL
    assert "cell-state-dependent methylation" in MODULE.GUARDRAIL
    assert "does not prove" in MODULE.GUARDRAIL
