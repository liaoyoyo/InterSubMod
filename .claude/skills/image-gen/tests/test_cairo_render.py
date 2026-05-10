"""Test cairo_render produces valid PNG for flow_diagram and icon."""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from cairo_render import render  # noqa

FIXTURES = Path(__file__).parent / "fixtures"


def is_png(path: Path) -> bool:
    return path.exists() and path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


def test_flow_diagram_renders():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "flow.png"
        render(FIXTURES / "valid_flow.yaml", out)
        assert is_png(out), f"Output is not a valid PNG: {out}"
        assert out.stat().st_size > 1000, f"PNG suspiciously small: {out.stat().st_size} bytes"


def test_icon_renders():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "icon.png"
        render(FIXTURES / "valid_icon.yaml", out)
        assert is_png(out), f"Output is not a valid PNG: {out}"


def test_rejects_non_cairo_type():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "bad.png"
        try:
            render(FIXTURES / "valid_concept.yaml", out)
            assert False, "Should reject backend=ai for cairo renderer"
        except ValueError as e:
            assert "backend" in str(e).lower(), f"Wrong error msg: {e}"


def main():
    tests = [test_flow_diagram_renders, test_icon_renders, test_rejects_non_cairo_type]
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
