"""Validate prompt YAML files conform to expected schema."""
import sys
from pathlib import Path
import yaml

SKILL_DIR = Path(__file__).parent.parent
TEMPLATES = SKILL_DIR / "templates"
FIXTURES = SKILL_DIR / "tests" / "fixtures"

REQUIRED_FIELDS = {"type", "subject", "constraints", "output"}
VALID_TYPES = {"concept_diagram", "flow_diagram", "data_mockup", "icon"}

def validate(path: Path) -> list[str]:
    errors = []
    data = yaml.safe_load(path.read_text())
    if not isinstance(data, dict):
        return [f"{path}: top-level not a dict"]
    missing = REQUIRED_FIELDS - data.keys()
    if missing:
        errors.append(f"{path}: missing fields {missing}")
    if data.get("type") not in VALID_TYPES:
        errors.append(f"{path}: type='{data.get('type')}' not in {VALID_TYPES}")
    out = data.get("output", {})
    if "size" not in out or "quality" not in out:
        errors.append(f"{path}: output.{{size,quality}} missing")
    return errors

def main():
    files = list(TEMPLATES.glob("*.yaml")) + list(FIXTURES.glob("valid_*.yaml"))
    if not files:
        print(f"FAIL: no yaml found under {TEMPLATES} or {FIXTURES}")
        return 1
    all_errors = []
    for f in files:
        all_errors.extend(validate(f))
    if all_errors:
        for e in all_errors:
            print(f"FAIL: {e}")
        return 1
    print(f"PASS: {len(files)} yaml files valid")
    return 0

if __name__ == "__main__":
    sys.exit(main())
