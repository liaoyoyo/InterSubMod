"""topic_folder_builder — build {basename}/ topic folder structure (D17).

Creates: foo/ + foo/prompts/ + foo/figures/ + foo/README.md
README is rendered from components/readme_template.md.jinja.
"""
from __future__ import annotations

import datetime
from pathlib import Path
from typing import Iterable

from jinja2 import Template

_SKILL_DIR = Path(__file__).resolve().parent.parent
_README_TEMPLATE_PATH = _SKILL_DIR / "components" / "readme_template.md.jinja"
SKILL_VERSION = "2.0.0"


def build_folder(md_path: Path, *, chapters: Iterable[dict], figures: Iterable[dict]) -> Path:
    """Build topic folder next to md_path. Returns folder Path.

    chapters: list of {"num": "01", "filename": "ch01_x.html", "title": "..."}
    figures: list of {"filename": "fig1.png", "status": "..."}
    """
    basename = md_path.stem
    folder = md_path.parent / basename
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "prompts").mkdir(exist_ok=True)
    (folder / "figures").mkdir(exist_ok=True)
    (folder / "figures" / ".gitkeep").touch(exist_ok=True)
    (folder / "prompts" / ".gitkeep").touch(exist_ok=True)

    tmpl_text = _README_TEMPLATE_PATH.read_text()
    rendered = Template(tmpl_text).render(
        basename=basename,
        skill_version=SKILL_VERSION,
        generated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        chapters=list(chapters),
        figures=list(figures),
    )
    (folder / "README.md").write_text(rendered)
    return folder


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 2:
        print("Usage: topic_folder_builder.py <md_file>")
        return 2
    md = Path(argv[1])
    folder = build_folder(md, chapters=[], figures=[])
    print(f"[topic_folder_builder] built {folder}")
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
