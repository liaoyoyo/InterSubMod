from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "summarize_exact_ps_topology.py"
)
SPEC = importlib.util.spec_from_file_location("summarize_exact_ps_topology", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_morphology_uses_rooted_graph_depth_and_branching() -> None:
    cases = {
        "no_branch_no_multilevel": [["ROOT", "A"]],
        "sister_only": [["ROOT", "A"], ["ROOT", "B"]],
        "direct_only": [["ROOT", "A"], ["A", "B"]],
        "sister_and_direct": [["ROOT", "A"], ["ROOT", "B"], ["A", "C"]],
    }
    for expected, edges in cases.items():
        assert MODULE.morphology_name(*MODULE.shape_features(edges)) == expected


def test_analysis_completeness_is_independent_of_display_truncation() -> None:
    unit = {
        "capped": False,
        "analysis_candidate_set_complete": True,
        "verification_status": "full_pass",
        "verification_complete": True,
        "verify_pass": True,
        "analysis_trees_generated": 44,
        "n_trees": 44,
        "n_trees_stored": 32,
        "n_distinct_shapes_exact": 4,
    }
    assert MODULE.is_complete(unit) is True
    assert MODULE.structural_class(unit) == "T>1|Topo>1"


def test_canonical_shape_ignores_node_labels_and_sibling_order() -> None:
    left = [["ROOT", "A"], ["ROOT", "B"], ["A", "C"]]
    right = [["X", "Z"], ["ROOT", "Y"], ["ROOT", "X"]]
    assert MODULE.canonical_shape(left) == MODULE.canonical_shape(right)
