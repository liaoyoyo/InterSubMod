"""Test dispatch end-to-end: md → companion HTML (single-file or topic-folder)."""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from dispatch import build_preview  # noqa

FIXTURES = Path(__file__).parent / "fixtures"


def test_short_report_produces_single_file():
    with tempfile.TemporaryDirectory() as tmp:
        # Copy short_report.md into tmp
        src = FIXTURES / "short_report.md"
        md = Path(tmp) / "short_report.md"
        md.write_text(src.read_text())

        result = build_preview(md)
        assert result["mode"] == "single_file", f"Got mode: {result['mode']}"
        out = result["output_path"]
        assert out.suffix == ".html"
        assert out.name == "short_report.preview.html"
        assert out.exists()
        html = out.read_text()
        assert "<title>Short Test Report</title>" in html
        assert "#D97757" in html, "design tokens missing"
        assert "<table>" in html, "table from md missing"


def test_large_report_produces_topic_folder():
    with tempfile.TemporaryDirectory() as tmp:
        src = FIXTURES / "large_report.md"
        md = Path(tmp) / "large_report.md"
        md.write_text(src.read_text())

        result = build_preview(md)
        assert result["mode"] == "topic_folder", f"Got mode: {result['mode']}"
        folder = result["output_path"]
        assert folder.is_dir()
        assert (folder / "index.html").exists()
        assert (folder / "README.md").exists()
        assert (folder / "figures").is_dir()
        assert (folder / "prompts").is_dir()


def test_taste_check_clean_passes():
    with tempfile.TemporaryDirectory() as tmp:
        src = FIXTURES / "short_report.md"
        md = Path(tmp) / "short_report.md"
        md.write_text(src.read_text())
        result = build_preview(md)
        assert result["taste_check"]["verdict"] == "clean", \
            f"Default template should be clean, got {result['taste_check']}"


def main():
    tests = [test_short_report_produces_single_file, test_large_report_produces_topic_folder,
             test_taste_check_clean_passes]
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
            import traceback
            traceback.print_exc()
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
