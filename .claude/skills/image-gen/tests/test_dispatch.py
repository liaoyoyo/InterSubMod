"""Test dispatch routes by backend field to correct subprocess."""
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from dispatch import route_one, plan_run  # noqa

FIXTURES = Path(__file__).parent / "fixtures"


def test_plan_run_classifies_correctly():
    plan = plan_run(FIXTURES, output_dir=Path("/tmp/out"))
    by_backend = {p["backend"]: 0 for p in plan}
    for p in plan:
        by_backend[p["backend"]] = by_backend.get(p["backend"], 0) + 1
    assert by_backend.get("ai", 0) == 2, f"Expected 2 AI, got {by_backend}"
    assert by_backend.get("cairo", 0) == 2, f"Expected 2 cairo, got {by_backend}"


def test_route_one_cairo_invokes_cairo_render():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "icon.png"
        route_one(FIXTURES / "valid_icon.yaml", out, dry_run=False)
        assert out.exists() and out.stat().st_size > 500, f"icon PNG missing/too small: {out}"


def test_route_one_ai_dry_run_does_not_call_codex():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "concept.png"
        with patch("dispatch.subprocess.run") as mock_run:
            route_one(FIXTURES / "valid_concept.yaml", out, dry_run=True)
            mock_run.assert_not_called()


def main():
    tests = [test_plan_run_classifies_correctly, test_route_one_cairo_invokes_cairo_render, test_route_one_ai_dry_run_does_not_call_codex]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS: {t.__name__}")
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
            failed += 1
        except Exception as e:
            print(f"ERROR: {t.__name__}: {type(e).__name__}: {e}")
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
