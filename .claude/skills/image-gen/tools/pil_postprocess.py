"""PIL post-processing: resize to standard sizes + base64 inline encoding."""
from __future__ import annotations

import base64
from pathlib import Path

from PIL import Image


def resize_to_target(src: Path, dst: Path, target_size: tuple[int, int]) -> None:
    """Resize src PNG to target_size. Preserves aspect ratio via thumbnail + padding (white bg)."""
    img = Image.open(src)
    img.thumbnail(target_size, Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", target_size, (255, 255, 255))
    paste_x = (target_size[0] - img.size[0]) // 2
    paste_y = (target_size[1] - img.size[1]) // 2
    if img.mode in ("RGBA", "LA"):
        canvas.paste(img, (paste_x, paste_y), img.convert("RGBA").split()[-1])
    else:
        canvas.paste(img, (paste_x, paste_y))
    canvas.save(dst, "PNG", optimize=True)


def to_base64(src: Path) -> str:
    """Return data:image/png;base64,... string."""
    return "data:image/png;base64," + base64.b64encode(Path(src).read_bytes()).decode("ascii")


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 3:
        print("Usage: pil_postprocess.py resize <src.png> <dst.png> [<W>x<H>]")
        print("       pil_postprocess.py base64 <src.png>")
        return 2
    cmd = argv[1]
    if cmd == "resize":
        if len(argv) < 4:
            print("Usage: pil_postprocess.py resize <src.png> <dst.png> [<W>x<H>]")
            return 2
        size = (1024, 1024)
        if len(argv) >= 5:
            w, h = argv[4].lower().split("x")
            size = (int(w), int(h))
        resize_to_target(Path(argv[2]), Path(argv[3]), size)
        print(f"[pil] resized → {argv[3]} ({size[0]}x{size[1]})")
        return 0
    if cmd == "base64":
        print(to_base64(Path(argv[2])))
        return 0
    print(f"Unknown cmd: {cmd}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
