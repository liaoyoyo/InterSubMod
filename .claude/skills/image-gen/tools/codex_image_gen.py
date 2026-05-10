"""Wrap codex CLI exec for AI-track image generation.

Reads prompt YAML, constructs codex exec command with $imagegen, runs it,
expects output PNG to land in --image-dir.

Backend: codex CLI v0.125+ (ChatGPT OAuth).
Model: gpt-image-2 (codex default for $imagegen as of 2026-04-21).
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import List

import yaml


def extract_prompt_text(yaml_path: Path) -> str:
    """Read prompt YAML and render to natural-language string for codex.

    Format: subject + key_elements (if any) + labels + style hints + constraints.
    Always appends '$imagegen' to trigger codex CLI's image gen skill.
    """
    data = yaml.safe_load(yaml_path.read_text())
    parts: list[str] = []

    parts.append(data.get("subject", ""))

    elems = data.get("key_elements") or []
    if elems:
        parts.append("Include: " + "; ".join(elems) + ".")

    labels = data.get("labels") or {}
    for pos, lab in labels.items():
        if lab:
            parts.append(f"Label at {pos}: '{lab}'.")

    chart = data.get("chart_type")
    if chart:
        ax = data.get("axes", {})
        parts.append(f"Chart type: {chart}. X axis: {ax.get('x_label','')}, Y axis: {ax.get('y_label','')}.")
        parts.append(f"Data trend hint: {data.get('data_hint','')}")

    style = data.get("style", {})
    pal = style.get("palette") or []
    if pal:
        parts.append(f"Palette: {', '.join(pal)}. Background: {style.get('background','white')}.")

    cons = data.get("constraints") or []
    if cons:
        parts.append("Constraints: " + "; ".join(cons) + ".")

    out = data.get("output", {})
    parts.append(f"Output size: {out.get('size','1024x1024')}, quality: {out.get('quality','high')}.")

    parts.append("$imagegen")
    return " ".join(p for p in parts if p)


def build_command(yaml_path: Path, output_dir: Path) -> List[str]:
    """Construct codex exec CLI command list."""
    prompt_text = extract_prompt_text(yaml_path)
    return [
        "codex", "exec",
        "--image-dir", str(output_dir),
        "--ask-for-approval", "never",
        prompt_text,
    ]


def run(yaml_path: Path, output_dir: Path, *, dry_run: bool = False) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = build_command(yaml_path, output_dir)
    if dry_run:
        print("DRY RUN:", " ".join(repr(c) for c in cmd))
        return 0
    print(f"[codex_image_gen] Running for {yaml_path.name} → {output_dir}/")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print("STDOUT:", proc.stdout, file=sys.stderr)
        print("STDERR:", proc.stderr, file=sys.stderr)
        print(f"[codex_image_gen] codex exec exited {proc.returncode}", file=sys.stderr)
    return proc.returncode


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print("Usage: codex_image_gen.py <prompt.yaml> <output_dir> [--dry-run]")
        return 2
    yaml_path = Path(argv[1])
    out_dir = Path(argv[2])
    dry = "--dry-run" in argv[3:]
    return run(yaml_path, out_dir, dry_run=dry)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
