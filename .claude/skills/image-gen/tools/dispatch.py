"""Main dispatch entry for image-gen skill.

Reads prompts from <prompts_dir>/*.yaml, routes by `backend` field:
  - ai    → codex_image_gen.run + codex_output_collector.collect_latest
  - cairo → cairo_render.render
Then pil_postprocess.resize_to_target standardizes output sizes.

Cost preview before run.
Idempotent: skips if target PNG exists (unless --force).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import yaml

# Add tools dir to import path so siblings work when run as script
_TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(_TOOLS))

import cairo_render  # noqa: E402
import codex_image_gen  # noqa: E402
import codex_output_collector  # noqa: E402
import pil_postprocess  # noqa: E402


def plan_run(prompts_dir: Path, output_dir: Path) -> list[dict]:
    """Return list of dicts {prompt_path, backend, target_path}."""
    plan = []
    for p in sorted(prompts_dir.glob("*.yaml")):
        data = yaml.safe_load(p.read_text())
        backend = data.get("backend", "ai")
        target = output_dir / (p.stem + ".png")
        plan.append({"prompt_path": p, "backend": backend, "target_path": target})
    return plan


def route_one(prompt_path: Path, target_path: Path, *, dry_run: bool = False) -> int:
    """Dispatch one prompt to its backend, write PNG to target_path."""
    data = yaml.safe_load(Path(prompt_path).read_text())
    backend = data.get("backend", "ai")
    target_path.parent.mkdir(parents=True, exist_ok=True)

    if backend == "cairo":
        if dry_run:
            print(f"DRY RUN [cairo]: {prompt_path.name} → {target_path}")
            return 0
        cairo_render.render(prompt_path, target_path)
        print(f"[cairo] {prompt_path.name} → {target_path}")
        return 0

    # AI backend: codex exec + collector
    if dry_run:
        cmd = codex_image_gen.build_command(prompt_path, target_path.parent)
        print(f"DRY RUN [ai]: {' '.join(repr(c) for c in cmd)}")
        return 0

    started_at = time.time()
    rc = codex_image_gen.run(prompt_path, target_path.parent)
    if rc != 0:
        print(f"[ai] codex exec failed for {prompt_path.name} (exit {rc})", file=sys.stderr)
        return rc
    # Rescue PNG if codex wrote elsewhere
    candidate_dirs = [target_path.parent] + codex_output_collector.CODEX_DEFAULT_DIRS
    moved = codex_output_collector.collect_latest(
        candidate_dirs, target_path.parent,
        target_name=target_path.name, since_epoch=started_at,
    )
    if moved is None:
        print(f"[ai] PNG not found in {candidate_dirs} for {prompt_path.name}", file=sys.stderr)
        return 1
    print(f"[ai] {prompt_path.name} → {moved}")
    return 0


def cost_preview(plan: list[dict]) -> str:
    n_ai = sum(1 for p in plan if p["backend"] == "ai")
    n_cairo = sum(1 for p in plan if p["backend"] == "cairo")
    return (
        f"[image-gen] {len(plan)} prompts: AI×{n_ai} (codex $imagegen, ~30k tokens/each, eats ChatGPT quota), "
        f"cairo×{n_cairo} (zero cost). Output dir: see plan."
    )


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="image-gen dispatch")
    parser.add_argument("prompts_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true",
                        help="re-render even if target PNG exists")
    parser.add_argument("--postprocess-size", default=None,
                        help="resize all outputs to e.g. 1024x1024")
    args = parser.parse_args(argv[1:])

    plan = plan_run(args.prompts_dir, args.output_dir)
    print(cost_preview(plan))
    if args.dry_run:
        for p in plan:
            print(f"  {p['backend']:5s}  {p['prompt_path'].name} → {p['target_path']}")
        return 0

    failed = 0
    for p in plan:
        if p["target_path"].exists() and not args.force:
            print(f"[skip] {p['target_path']} (use --force to re-render)")
            continue
        rc = route_one(p["prompt_path"], p["target_path"])
        if rc != 0:
            failed += 1

    if args.postprocess_size:
        w, h = args.postprocess_size.lower().split("x")
        size = (int(w), int(h))
        for p in plan:
            if p["target_path"].exists():
                pil_postprocess.resize_to_target(p["target_path"], p["target_path"], size)
                print(f"[pil] resized {p['target_path']} → {size}")

    if failed:
        print(f"[image-gen] {failed} of {len(plan)} prompts failed", file=sys.stderr)
        return 1
    print(f"[image-gen] {len(plan)} prompts succeeded")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
