#!/usr/bin/env python3
"""Behavioral tests for scripts/fill_report.py (anti-fabrication Layer A).

Audit 2026-06-10 finding PYTOOL-1: a present-but-null (or NaN) data value was
rendered as the literal string 'None' / 'nan' instead of being refused, so a
"no value computed" cell silently became fabricated-looking text in a report.
The contract is: every {{key}} must resolve to a REAL value; null/NaN counts as
"no source" and must trigger the same REFUSE path as a missing key.

Pure-stdlib (no pytest): run with `python3 scripts/tests/test_fill_report.py`.
Exercises the real CLI via subprocess (no mocks).
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

FILL = Path(__file__).resolve().parents[1] / "fill_report.py"

_passed = 0
_failed = 0


def _run(template_text, data_obj, extra_args=None, allow_nan=False):
    """Render template_text with data_obj; return (returncode, out_text_or_None, stderr)."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        tmpl = d / "t.tmpl"
        data = d / "d.json"
        out = d / "out.txt"
        tmpl.write_text(template_text, encoding="utf-8")
        data.write_text(json.dumps(data_obj, allow_nan=allow_nan), encoding="utf-8")
        cmd = [sys.executable, str(FILL), str(tmpl), str(data), "-o", str(out)]
        if extra_args:
            cmd += extra_args
        p = subprocess.run(cmd, capture_output=True, text=True)
        out_text = out.read_text(encoding="utf-8") if out.exists() else None
        return p.returncode, out_text, p.stderr


def check(name, cond, detail=""):
    global _passed, _failed
    if cond:
        _passed += 1
        print(f"  PASS  {name}")
    else:
        _failed += 1
        print(f"  FAIL  {name}  {detail}")


def test_null_value_is_refused():
    rc, out, err = _run("AUC={{auc}}", {"auc": None})
    check("null value -> exit 1 (refuse)", rc == 1, f"rc={rc}")
    check("null value -> no output written", out is None, f"out={out!r}")
    check("null value -> never renders literal 'None'", out is None or "None" not in out,
          f"out={out!r}")


def test_nan_value_is_refused():
    rc, out, err = _run("AUC={{auc}}", {"auc": float("nan")}, allow_nan=True)
    check("NaN value -> exit 1 (refuse)", rc == 1, f"rc={rc}")
    check("NaN value -> never renders literal 'nan'", out is None or "nan" not in out,
          f"out={out!r}")


def test_valid_value_still_renders():
    rc, out, err = _run("AUC={{auc}}", {"auc": 0.85})
    check("valid value -> exit 0", rc == 0, f"rc={rc} err={err}")
    check("valid value -> rendered", out is not None and "0.85" in out, f"out={out!r}")


def test_zero_and_false_still_render():
    # 0 and False are legitimate computed values, NOT missing — must render.
    rc, out, err = _run("n={{n}} flag={{flag}}", {"n": 0, "flag": False})
    check("zero/false -> exit 0 (not treated as missing)", rc == 0, f"rc={rc} err={err}")
    check("zero renders as 0", out is not None and "n=0" in out, f"out={out!r}")
    check("false renders", out is not None and "flag=False" in out, f"out={out!r}")


def test_null_with_allow_missing_is_marked():
    rc, out, err = _run("AUC={{auc}}", {"auc": None}, extra_args=["--allow-missing"])
    check("null + --allow-missing -> exit 0 (draft)", rc == 0, f"rc={rc} err={err}")
    check("null + --allow-missing -> {{MISSING:auc}} marker, not 'None'",
          out is not None and "MISSING:auc" in out and "None" not in out, f"out={out!r}")


if __name__ == "__main__":
    print(f"[test_fill_report] target: {FILL}")
    test_null_value_is_refused()
    test_nan_value_is_refused()
    test_valid_value_still_renders()
    test_zero_and_false_still_render()
    test_null_with_allow_missing_is_marked()
    print(f"\n[test_fill_report] {_passed} passed, {_failed} failed")
    sys.exit(1 if _failed else 0)
