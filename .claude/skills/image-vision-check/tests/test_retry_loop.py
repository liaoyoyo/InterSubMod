"""Test apply_suggestions_to_yaml appends suggestions to constraints field."""
import sys
import tempfile
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from vision_check_runner import apply_suggestions_to_yaml  # noqa


def test_appends_to_constraints():
    with tempfile.TemporaryDirectory() as tmp:
        yaml_path = Path(tmp) / "p.yaml"
        yaml_path.write_text(yaml.safe_dump({
            "type": "concept_diagram",
            "subject": "x",
            "constraints": ["original"],
            "output": {"size": "1024x1024", "quality": "high"},
            "backend": "ai",
        }))
        suggestions = ["constraints: add 'no drop shadow'"]
        apply_suggestions_to_yaml(yaml_path, suggestions)
        new = yaml.safe_load(yaml_path.read_text())
        assert "no drop shadow" in new["constraints"], f"Constraint not added: {new['constraints']}"
        assert "original" in new["constraints"], "Original constraint dropped"


def test_handles_missing_constraints_field():
    with tempfile.TemporaryDirectory() as tmp:
        yaml_path = Path(tmp) / "p.yaml"
        yaml_path.write_text(yaml.safe_dump({
            "type": "icon",
            "subject": "x",
            "output": {"size": "256x256", "quality": "n/a"},
            "backend": "cairo",
        }))
        apply_suggestions_to_yaml(yaml_path, ["constraints: add 'high contrast glyph'"])
        new = yaml.safe_load(yaml_path.read_text())
        assert "constraints" in new and "high contrast glyph" in new["constraints"]


def main():
    failed = 0
    for t in [test_appends_to_constraints, test_handles_missing_constraints_field]:
        try:
            t()
            print(f"PASS: {t.__name__}")
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
