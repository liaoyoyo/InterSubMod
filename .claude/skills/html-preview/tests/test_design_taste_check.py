"""Test design_taste_check static detection (regex-based, no Claude call)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from design_taste_check import detect_taboos_static  # noqa


def test_clean_html_passes():
    html = """
    <style>body { color: #141413; background: #FAF9F5; }</style>
    <h1>Title</h1>
    """
    result = detect_taboos_static(html)
    assert result["n_taboos"] == 0
    assert result["verdict"] == "clean"


def test_gradient_overuse_flagged():
    html = """
    <style>
    .a { background: linear-gradient(red, blue); }
    .b { background: linear-gradient(green, yellow); }
    </style>
    """
    result = detect_taboos_static(html)
    assert result["taboos_detected"]["gradient_overuse"] is True
    assert result["n_taboos"] >= 1


def test_glass_morphism_flagged():
    html = '<style>.x { backdrop-filter: blur(10px); }</style>'
    result = detect_taboos_static(html)
    assert result["taboos_detected"]["glass_morphism"] is True


def test_emoji_header_flagged():
    html = '<h1>🚀 Welcome 🎉</h1>'
    result = detect_taboos_static(html)
    assert result["taboos_detected"]["emoji_decorated_headers"] is True


def test_drop_shadow_text_flagged():
    html = '<style>h1 { text-shadow: 2px 2px 4px black; }</style>'
    result = detect_taboos_static(html)
    assert result["taboos_detected"]["drop_shadow_text"] is True


def test_minor_issues_verdict_one_taboo():
    html = '<style>.x { backdrop-filter: blur(10px); }</style>'
    result = detect_taboos_static(html)
    assert result["verdict"] == "minor_issues"


def test_redesign_verdict_two_taboos():
    html = """
    <style>
    .x { backdrop-filter: blur(10px); }
    h1 { text-shadow: 2px 2px 4px black; }
    </style>
    """
    result = detect_taboos_static(html)
    assert result["verdict"] == "redesign"
    assert result["n_taboos"] >= 2


def main():
    tests = [test_clean_html_passes, test_gradient_overuse_flagged, test_glass_morphism_flagged,
             test_emoji_header_flagged, test_drop_shadow_text_flagged,
             test_minor_issues_verdict_one_taboo, test_redesign_verdict_two_taboos]
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
