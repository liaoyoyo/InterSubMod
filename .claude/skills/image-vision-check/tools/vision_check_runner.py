"""vision_check_runner — enumerate images, validate Claude-produced quality JSON, aggregate.

Note on architecture:
This runner does NOT itself call Claude vision. The pattern is:
1. Runner enumerates PNGs and writes `<figures_dir>/_check_instructions.json`.
2. Calling Claude agent reads the instructions, uses its own `Read` tool on each
   PNG, and writes per-image JSON entries (matching `schemas/quality.schema.json`)
   into `<figures_dir>/_pending_entries/<basename>.json`.
3. Runner is re-invoked with `--aggregate` to validate + merge into final
   `<figures_dir>/quality.json`.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

# Optional jsonschema; we degrade gracefully if missing.
try:
    import jsonschema as _js  # type: ignore
    _HAS_JSONSCHEMA = True
except ImportError:
    _HAS_JSONSCHEMA = False


def enumerate_images(figures_dir: Path) -> list[Path]:
    return sorted(figures_dir.glob("*.png"))


def validate_entry(entry: dict, schema_path: Path) -> list[str]:
    """Return list of error strings; empty list means valid."""
    if not _HAS_JSONSCHEMA:
        # Manual minimal check
        required = {"image", "version", "checks", "score", "verdict"}
        missing = required - entry.keys()
        if missing:
            return [f"missing required fields: {missing}"]
        score = entry.get("score", "")
        if not (
            isinstance(score, str)
            and len(score) == 3
            and score[1] == "/"
            and score[0].isdigit()
            and score[2] == "6"
            and 0 <= int(score[0]) <= 6
        ):
            return [f"invalid score: {score}"]
        return []
    schema = json.loads(Path(schema_path).read_text())
    validator = _js.Draft202012Validator(schema)
    return [
        f"{e.message} at /{'/'.join(map(str, e.path))}"
        for e in validator.iter_errors(entry)
    ]


def aggregate(entries: list[dict]) -> dict:
    n_pass = sum(1 for e in entries if e.get("verdict") == "pass")
    n_partial = sum(1 for e in entries if e.get("verdict") == "partial")
    n_fail = sum(1 for e in entries if e.get("verdict") == "fail")
    n_human = sum(1 for e in entries if e.get("verdict") == "needs_human")
    overall = "all_pass" if (n_partial + n_fail + n_human) == 0 else "mixed"
    return {
        "n_pass": n_pass,
        "n_partial": n_partial,
        "n_fail": n_fail,
        "n_needs_human": n_human,
        "overall": overall,
    }


def write_instructions(figures_dir: Path, prompts_dir: Path | None) -> Path:
    images = enumerate_images(figures_dir)
    instructions = {
        "figures_dir": str(figures_dir.resolve()),
        "prompts_dir": str(prompts_dir.resolve()) if prompts_dir else None,
        "schema_path": str(
            (Path(__file__).parent.parent / "schemas" / "quality.schema.json").resolve()
        ),
        "master_prompt": str(
            (Path(__file__).parent.parent / "prompts" / "vision_review_master.md").resolve()
        ),
        "checklists_dir": str(
            (Path(__file__).parent.parent / "checklists").resolve()
        ),
        "pending_dir": str((figures_dir / "_pending_entries").resolve()),
        "images": [
            {
                "png": str(img.resolve()),
                "yaml": str((prompts_dir / (img.stem + ".yaml")).resolve())
                if prompts_dir and (prompts_dir / (img.stem + ".yaml")).exists()
                else None,
            }
            for img in images
        ],
    }
    out = figures_dir / "_check_instructions.json"
    out.write_text(json.dumps(instructions, indent=2, ensure_ascii=False))
    (figures_dir / "_pending_entries").mkdir(exist_ok=True)
    return out


def aggregate_pending(figures_dir: Path) -> Path:
    schema_path = Path(__file__).parent.parent / "schemas" / "quality.schema.json"
    pending = figures_dir / "_pending_entries"
    if not pending.exists():
        print(f"No pending entries at {pending}", file=sys.stderr)
        sys.exit(2)

    entries: list[dict] = []
    errors: list[str] = []
    for j in sorted(pending.glob("*.json")):
        try:
            entry = json.loads(j.read_text())
        except json.JSONDecodeError as e:
            errors.append(f"{j.name}: {e}")
            continue
        errs = validate_entry(entry, schema_path)
        if errs:
            errors.append(f"{j.name}: {'; '.join(errs)}")
            continue
        entries.append(entry)

    if errors:
        print("Validation errors:")
        for e in errors:
            print(f"  {e}")

    summary = aggregate(entries)
    quality = {
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "figures_dir": str(figures_dir.resolve()),
        "summary": summary,
        "entries": entries,
    }
    out = figures_dir / "quality.json"
    out.write_text(json.dumps(quality, indent=2, ensure_ascii=False))
    print(f"[vision_check] aggregated {len(entries)} entries → {out}")
    print(f"[vision_check] summary: {summary}")
    return out


def apply_suggestions_to_yaml(yaml_path: Path, suggestions: list[str]) -> None:
    """Apply textual suggestions from quality.json to a prompt YAML.

    Currently supports the 'constraints: add ...' suggestion form.
    Other suggestion forms are logged but not auto-applied (D8 max 1 retry,
    keep simple).
    """
    import yaml as _yaml
    data = _yaml.safe_load(Path(yaml_path).read_text())
    if not isinstance(data, dict):
        return
    constraints = list(data.get("constraints") or [])
    for s in suggestions:
        prefix = "constraints: add"
        if s.startswith(prefix):
            tail = s[len(prefix):].strip()
            tail = tail.strip("'\"`")
            if tail and tail not in constraints:
                constraints.append(tail)
    data["constraints"] = constraints
    Path(yaml_path).write_text(_yaml.safe_dump(data, allow_unicode=True, sort_keys=False))


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="vision_check_runner")
    parser.add_argument("figures_dir", type=Path)
    parser.add_argument(
        "--prompts",
        type=Path,
        default=None,
        help="paired prompts dir (defaults to <figures_dir>/../prompts/)",
    )
    parser.add_argument(
        "--aggregate",
        action="store_true",
        help="aggregate pending JSON entries into final quality.json",
    )
    args = parser.parse_args(argv[1:])

    if not args.figures_dir.is_dir():
        print(f"figures_dir not found: {args.figures_dir}", file=sys.stderr)
        return 2

    if args.aggregate:
        aggregate_pending(args.figures_dir)
        return 0

    prompts_dir = args.prompts
    if prompts_dir is None:
        candidate = args.figures_dir.parent / "prompts"
        prompts_dir = candidate if candidate.is_dir() else None

    instructions = write_instructions(args.figures_dir, prompts_dir)
    print(f"[vision_check] instructions written → {instructions}")
    print(f"[vision_check] next: Claude agent reads instructions, processes each image,")
    print(f"[vision_check] writes to _pending_entries/<basename>.json,")
    print(f"[vision_check] then re-run this script with --aggregate")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
