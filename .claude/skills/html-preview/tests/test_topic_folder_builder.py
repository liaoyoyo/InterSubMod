"""Test topic_folder_builder creates folder + README."""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from topic_folder_builder import build_folder  # noqa


def test_creates_folder_with_readme():
    with tempfile.TemporaryDirectory() as tmp:
        md = Path(tmp) / "test_report.md"
        md.write_text("---\ntitle: Test\n---\n\n# Test\n")
        folder = build_folder(md, chapters=[
            {"num": "01", "filename": "ch01_intro.html", "title": "Intro"},
            {"num": "02", "filename": "ch02_methods.html", "title": "Methods"},
        ], figures=[])
        assert folder == md.parent / "test_report"
        assert folder.exists() and folder.is_dir()
        readme = folder / "README.md"
        assert readme.exists()
        readme_text = readme.read_text()
        assert "test_report" in readme_text
        assert "ch01_intro.html" in readme_text
        assert "Methods" in readme_text


def test_creates_subdirs_prompts_figures():
    with tempfile.TemporaryDirectory() as tmp:
        md = Path(tmp) / "x.md"
        md.write_text("# x")
        folder = build_folder(md, chapters=[], figures=[])
        assert (folder / "prompts").is_dir()
        assert (folder / "figures").is_dir()


def test_includes_figures_in_readme():
    with tempfile.TemporaryDirectory() as tmp:
        md = Path(tmp) / "x.md"
        md.write_text("# x")
        folder = build_folder(md, chapters=[], figures=[
            {"filename": "fig1.png", "status": "✓ generated"},
            {"filename": "fig2.png", "status": "pending"},
        ])
        readme_text = (folder / "README.md").read_text()
        assert "fig1.png" in readme_text
        assert "pending" in readme_text


def main():
    tests = [test_creates_folder_with_readme, test_creates_subdirs_prompts_figures, test_includes_figures_in_readme]
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
