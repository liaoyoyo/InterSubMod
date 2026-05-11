"""Import-friendly alias of 00_feature_registry.py.

The file ``00_feature_registry.py`` starts with a digit, so Python cannot
import it as a module (``import 00_feature_registry`` is a SyntaxError).
This shim re-exports its public API so downstream Phase C agents can do:

    from research.feature_layered_observation.scripts.feature_registry import (
        get_feature_metadata, list_features_by_group, load_registry,
        get_group_summary, features_needing_source_lookup, GROUP_NAMES,
    )

The ``00_`` CLI script remains the authoritative source and entry point; this
module loads it via importlib to keep a single code path.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

_SPEC_PATH = Path(__file__).with_name("00_feature_registry.py")
_spec = importlib.util.spec_from_file_location("_fr_impl", _SPEC_PATH)
if _spec is None or _spec.loader is None:  # pragma: no cover
    raise ImportError(f"Cannot load registry impl from {_SPEC_PATH}")
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

# Re-export public API
GROUP_NAMES = _mod.GROUP_NAMES
REGISTRY_PATH = _mod.REGISTRY_PATH
load_registry = _mod.load_registry
get_feature_metadata = _mod.get_feature_metadata
list_features_by_group = _mod.list_features_by_group
get_group_summary = _mod.get_group_summary
features_needing_source_lookup = _mod.features_needing_source_lookup

__all__ = [
    "GROUP_NAMES",
    "REGISTRY_PATH",
    "load_registry",
    "get_feature_metadata",
    "list_features_by_group",
    "get_group_summary",
    "features_needing_source_lookup",
]
