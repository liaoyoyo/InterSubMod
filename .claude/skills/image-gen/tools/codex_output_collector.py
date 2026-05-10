"""Find newest PNG across candidate dirs (created after `since_epoch`) and move it.

Use case: codex CLI may write generated PNG to `--image-dir` OR fall back to
`~/.codex/generated_images/`. This collector scans both, picks the latest PNG
created after the codex call started, and renames it to a target filename.
"""
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Iterable, Optional


def collect_latest(
    candidate_dirs: Iterable[Path],
    target_dir: Path,
    *,
    target_name: str,
    since_epoch: float,
) -> Optional[Path]:
    """Return path to moved file, or None if no eligible PNG found.

    `since_epoch` filters PNGs created at-or-after that mtime.
    """
    target_dir.mkdir(parents=True, exist_ok=True)
    candidates: list[Path] = []
    for d in candidate_dirs:
        if not d.exists():
            continue
        candidates.extend(p for p in d.glob("*.png") if p.stat().st_mtime >= since_epoch)
    if not candidates:
        return None
    latest = max(candidates, key=lambda p: p.stat().st_mtime)
    dst = target_dir / target_name
    shutil.move(str(latest), str(dst))
    return dst


CODEX_DEFAULT_DIRS = [
    Path.home() / ".codex" / "generated_images",
]


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 4:
        print("Usage: codex_output_collector.py <target_dir> <target_name> <since_epoch> [<extra_dir>...]")
        return 2
    target_dir = Path(argv[1])
    target_name = argv[2]
    since_epoch = float(argv[3])
    extras = [Path(a) for a in argv[4:]]
    dirs = extras + CODEX_DEFAULT_DIRS
    result = collect_latest(dirs, target_dir, target_name=target_name, since_epoch=since_epoch)
    if result is None:
        print("[collector] No PNG found in any candidate dir", file=sys.stderr)
        return 1
    print(f"[collector] Moved → {result}")
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
