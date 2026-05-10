"""Test interactivity_detect: returns L2 (default) or L3 (interactive)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from interactivity_detect import detect_level  # noqa

FIXTURES = Path(__file__).parent / "fixtures"


def test_short_report_is_l2():
    p = FIXTURES / "short_report.md"
    assert detect_level(p) == "l2", "Default should be L2"


def test_interactive_marker_is_l3():
    p = FIXTURES / "interactive_marker.md"
    assert detect_level(p) == "l3", "<!-- interactive: tab --> should trigger L3"


def test_explicit_l1_override():
    # Tests that --force-l1 override works
    import tempfile
    text = "<!-- html-preview: force-l1 -->\n# x"
    with tempfile.NamedTemporaryFile(suffix=".md", mode="w", delete=False) as f:
        f.write(text)
        f.flush()
        p = Path(f.name)
    try:
        assert detect_level(p) == "l1"
    finally:
        p.unlink()


def main():
    tests = [test_short_report_is_l2, test_interactive_marker_is_l3, test_explicit_l1_override]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS: {t.__name__}")
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
