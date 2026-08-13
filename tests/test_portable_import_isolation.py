#!/usr/bin/env python3
"""Prevent the clean test suite from importing Python modules from another clone."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_historical_consumers_resolve_only_from_current_checkout() -> None:
    probe = r"""
import importlib
import json
import sys
from pathlib import Path

repo = Path(sys.argv[1]).resolve()
module_names = (
    "scripts.analysis.build_observation_O06_verification_class",
    "scripts.analysis.build_observation_O10_read_level_methyl",
    "scripts.analysis.build_observation_O15b_loh_zone_metrics_cross_sample",
    "scripts.analysis.build_phase1a_head_to_head_baseline_table",
    "scripts.analysis.fn_verdict_readback_audit",
    "scripts.analysis.tp_fp_structure_label_association",
    "scripts.validation.metrics.compute_metrics",
)
for name in module_names:
    importlib.import_module(name)

outside_modules = {}
for name, module in sorted(sys.modules.items()):
    if not name.startswith("scripts."):
        continue
    filename = getattr(module, "__file__", None)
    if filename is None:
        continue
    resolved = Path(filename).resolve()
    if not resolved.is_relative_to(repo):
        outside_modules[name] = str(resolved)

foreign_repo_entries = []
for entry in sys.path:
    if not entry:
        continue
    candidate = Path(entry)
    if not candidate.is_absolute():
        continue
    resolved = candidate.resolve()
    if resolved == repo or resolved.is_relative_to(repo):
        continue
    if resolved.name == "InterSubMod" or (
        resolved.name == "analysis" and resolved.parent.name == "scripts"
    ):
        foreign_repo_entries.append(str(resolved))

print(json.dumps({
    "outside_modules": outside_modules,
    "foreign_repo_entries": sorted(set(foreign_repo_entries)),
}, sort_keys=True))
raise SystemExit(1 if outside_modules or foreign_repo_entries else 0)
"""
    completed = subprocess.run(
        [sys.executable, "-c", probe, str(REPO_ROOT)],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    payload = json.loads(completed.stdout)
    assert completed.returncode == 0, completed.stderr or payload
    assert payload == {"foreign_repo_entries": [], "outside_modules": {}}
