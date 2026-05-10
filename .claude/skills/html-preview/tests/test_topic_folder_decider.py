"""Test topic_folder_decider applies D17 logic."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from topic_folder_decider import should_use_topic_folder  # noqa

FIXTURES = Path(__file__).parent / "fixtures"


def test_short_report_single_file():
    short = FIXTURES / "short_report.md"
    assert should_use_topic_folder(short) is False, "Short report (~30 lines, 0 figs) should be single-file"


def test_medium_report_single_file():
    medium = FIXTURES / "medium_report.md"
    # 180 lines, 3 figures — should be SINGLE FILE per D17 (need >=200 lines OR >=5 figs)
    assert should_use_topic_folder(medium) is False, \
        "Medium (180 lines, 3 figs) should be single-file"


def test_large_report_topic_folder():
    large = FIXTURES / "large_report.md"
    # 250 lines, 5 figures — exceeds both thresholds
    assert should_use_topic_folder(large) is True, \
        "Large (250 lines, 5 figs) should be topic-folder"


def test_lines_threshold_only():
    # Hypothetical: 220 lines, 0 figures
    text = "---\ntitle: x\n---\n\n# x\n\n" + "\n".join(f"line {i}" for i in range(220))
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".md", mode="w", delete=False) as f:
        f.write(text)
        f.flush()
        path = Path(f.name)
    try:
        assert should_use_topic_folder(path) is True, ">200 lines should trigger"
    finally:
        path.unlink()


def test_figures_threshold_only():
    # Hypothetical: 30 lines, 6 figures
    text = "---\ntitle: x\n---\n\n# x\n\n" + "\n".join(f"![](f{i}.png)" for i in range(6))
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".md", mode="w", delete=False) as f:
        f.write(text)
        f.flush()
        path = Path(f.name)
    try:
        assert should_use_topic_folder(path) is True, ">=5 figures should trigger"
    finally:
        path.unlink()


def main():
    tests = [test_short_report_single_file, test_medium_report_single_file,
             test_large_report_topic_folder, test_lines_threshold_only, test_figures_threshold_only]
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
