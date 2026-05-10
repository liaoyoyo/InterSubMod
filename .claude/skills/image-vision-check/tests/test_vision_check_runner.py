"""Test runner: enumerate, validate JSON, aggregate."""
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from vision_check_runner import (  # noqa
    enumerate_images,
    validate_entry,
    aggregate,
)

SCHEMA_PATH = Path(__file__).parent.parent / "schemas" / "quality.schema.json"


def make_test_png(path: Path):
    from PIL import Image
    Image.new("RGB", (32, 32), (200, 100, 50)).save(path, "PNG")


def test_enumerate_images_finds_pngs():
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        make_test_png(d / "fig1.png")
        make_test_png(d / "fig2.png")
        (d / "not_a_png.txt").write_text("ignored")
        out = enumerate_images(d)
        names = sorted(p.name for p in out)
        assert names == ["fig1.png", "fig2.png"], f"Wrong: {names}"


def test_validate_entry_pass():
    entry = {
        "image": "fig1.png",
        "version": 1,
        "type": "concept_diagram",
        "backend": "ai",
        "checks": {
            "prompt_fidelity": True,
            "readability": True,
            "focal_clarity": True,
            "design_taboos_avoided": True,
            "print_friendly": True,
            "cross_figure_consistency": True,
        },
        "score": "6/6",
        "verdict": "pass",
        "suggestions": [],
        "checked_at": "2026-05-10T14:00:00Z",
    }
    assert validate_entry(entry, SCHEMA_PATH) == [], "Should be valid"


def test_validate_entry_fail_bad_score():
    entry = {
        "image": "fig1.png",
        "version": 1,
        "checks": {
            "prompt_fidelity": True,
            "readability": True,
            "focal_clarity": True,
            "design_taboos_avoided": True,
            "print_friendly": True,
            "cross_figure_consistency": True,
        },
        "score": "7/6",  # invalid
        "verdict": "pass",
    }
    errors = validate_entry(entry, SCHEMA_PATH)
    assert errors, "Should report invalid score format"


def test_aggregate_summary():
    entries = [
        {"verdict": "pass"},
        {"verdict": "pass"},
        {"verdict": "partial"},
        {"verdict": "fail"},
    ]
    s = aggregate(entries)
    assert s["n_pass"] == 2 and s["n_partial"] == 1 and s["n_fail"] == 1
    assert s["overall"] == "mixed"
    s2 = aggregate([{"verdict": "pass"}, {"verdict": "pass"}])
    assert s2["overall"] == "all_pass"


def main():
    tests = [test_enumerate_images_finds_pngs, test_validate_entry_pass, test_validate_entry_fail_bad_score, test_aggregate_summary]
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
