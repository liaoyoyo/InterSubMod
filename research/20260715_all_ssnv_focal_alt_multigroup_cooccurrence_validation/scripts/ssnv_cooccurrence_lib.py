#!/usr/bin/env python3
"""Compatibility bridge to the repository-canonical sSNV co-occurrence library.

The formal producer and all topic-local tests must execute one implementation.
The canonical source of truth is InterSubMod/scripts/ssnv_cooccurrence_lib.py.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


_CANONICAL_PATH = (
    Path(__file__).resolve().parents[3] / "scripts" / "ssnv_cooccurrence_lib.py"
)
_MODULE_NAME = "_intersubmod_canonical_ssnv_cooccurrence_lib"

if _MODULE_NAME in sys.modules:
    _CANONICAL = sys.modules[_MODULE_NAME]
else:
    _SPEC = importlib.util.spec_from_file_location(_MODULE_NAME, _CANONICAL_PATH)
    if _SPEC is None or _SPEC.loader is None:
        raise ImportError(f"Cannot load canonical co-occurrence library: {_CANONICAL_PATH}")
    _CANONICAL = importlib.util.module_from_spec(_SPEC)
    sys.modules[_MODULE_NAME] = _CANONICAL
    _SPEC.loader.exec_module(_CANONICAL)

__all__ = (
    "Variant",
    "call_snv",
    "call_snvs",
    "cramer_v",
    "contingency_table",
    "group_allele_table",
    "methyl_group_marker_association",
    "fisher_freeman_halton_kx2",
    "categorical_association_table",
    "conditional_label_permutation",
    "four_state_summary",
    "benjamini_hochberg",
    "benjamini_yekutieli",
    "spatially_separated_positions",
    "count_spatially_separated_positions",
    "count_independent_positions",
)

for _PUBLIC_NAME in __all__:
    globals()[_PUBLIC_NAME] = getattr(_CANONICAL, _PUBLIC_NAME)

CANONICAL_SOURCE = _CANONICAL_PATH
