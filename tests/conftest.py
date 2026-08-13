"""Repository-wide pytest guards for clean-checkout provenance."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def _outside_scripts_modules() -> dict[str, str]:
    outside: dict[str, str] = {}
    for name, module in sorted(sys.modules.items()):
        if not name.startswith("scripts."):
            continue
        filename = getattr(module, "__file__", None)
        if filename is None:
            continue
        resolved = Path(filename).resolve()
        if not resolved.is_relative_to(REPO_ROOT):
            outside[name] = str(resolved)
    return outside


def _foreign_repo_sys_path_entries() -> list[str]:
    foreign: set[str] = set()
    for entry in sys.path:
        if not entry:
            continue
        candidate = Path(entry)
        if not candidate.is_absolute():
            continue
        resolved = candidate.resolve()
        if resolved == REPO_ROOT or resolved.is_relative_to(REPO_ROOT):
            continue
        if resolved.name == "InterSubMod" or (
            resolved.name == "analysis" and resolved.parent.name == "scripts"
        ):
            foreign.add(str(resolved))
    return sorted(foreign)


def pytest_sessionfinish(session: pytest.Session, exitstatus: int) -> None:
    """Fail if tests imported code from a different InterSubMod checkout."""

    outside = _outside_scripts_modules()
    foreign = _foreign_repo_sys_path_entries()
    if not outside and not foreign:
        return

    reporter = session.config.pluginmanager.get_plugin("terminalreporter")
    if reporter is not None:
        reporter.write_sep("=", "CLEAN CHECKOUT IMPORT ISOLATION FAILED", red=True)
        for name, path in outside.items():
            reporter.write_line(f"outside module: {name} -> {path}", red=True)
        for path in foreign:
            reporter.write_line(f"foreign sys.path entry: {path}", red=True)
    session.exitstatus = pytest.ExitCode.TESTS_FAILED
