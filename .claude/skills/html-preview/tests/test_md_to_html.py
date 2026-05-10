"""Test md_to_html: markdown -> HTML body."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from md_to_html import convert, parse_frontmatter  # noqa

FIXTURES = Path(__file__).parent / "fixtures"


def test_parse_frontmatter():
    md = (FIXTURES / "short_report.md").read_text()
    fm, body = parse_frontmatter(md)
    assert fm.get("title") == "Short Test Report", f"Bad fm: {fm}"
    assert fm.get("date") == "2026-05-11", f"Bad date: {fm.get('date')}"
    assert "# Short Test Report" in body, "Body missing H1"


def test_convert_basic():
    md = (FIXTURES / "short_report.md").read_text()
    html = convert(md)
    assert "<h1>Short Test Report</h1>" in html
    assert "<table>" in html, "Table not rendered"
    assert "<code>" in html, "Inline code not rendered"
    assert "<pre>" in html, "Code block not rendered"
    assert "<strong>bold</strong>" in html, "Bold not rendered"
    assert '<a href="https://example.com">link to docs</a>' in html


def test_convert_strips_frontmatter():
    md = (FIXTURES / "short_report.md").read_text()
    html = convert(md)
    # Frontmatter (date, title in YAML form) should NOT appear in body
    assert "---" not in html, "YAML separator leaked into HTML"


def test_convert_handles_tables():
    md = "| A | B |\n|---|---|\n| 1 | 2 |\n"
    html = convert(md)
    assert "<table>" in html
    assert "<th>A</th>" in html
    assert "<td>1</td>" in html


def test_convert_handles_details():
    # GitHub-flavored details survive
    md = "<details><summary>Title</summary>\n\nbody\n\n</details>"
    html = convert(md)
    assert "<details>" in html
    assert "<summary>Title</summary>" in html


def main():
    tests = [
        test_parse_frontmatter, test_convert_basic, test_convert_strips_frontmatter,
        test_convert_handles_tables, test_convert_handles_details,
    ]
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
