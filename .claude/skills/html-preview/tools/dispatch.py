"""dispatch — main entry for html-preview skill.

Workflow:
  1. Detect topic-folder vs single-file (D17)
  2. Detect L2 vs L3 interactivity level (D15 default L2)
  3. Convert .md body → HTML (md_to_html)
  4. Render base_lN.html.jinja with design tokens inline
  5. Inline images (PNG → base64)
  6. Static taste audit (D14/D18)
  7. Write output (single .preview.html OR topic folder index.html)
"""
from __future__ import annotations

import datetime
import sys
from pathlib import Path

from jinja2 import Template

_TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(_TOOLS))

import design_taste_check  # noqa: E402
import inline_assets  # noqa: E402
import interactivity_detect  # noqa: E402
import md_to_html  # noqa: E402
import topic_folder_builder  # noqa: E402
import topic_folder_decider  # noqa: E402

_SKILL_DIR = _TOOLS.parent
_TEMPLATES_DIR = _SKILL_DIR / "templates"
_DESIGN_TOKENS_PATH = _TEMPLATES_DIR / "design_tokens.css"
SKILL_VERSION = "2.0.0"


def _render(md_path: Path, level: str) -> str:
    """Render md to full HTML string at given level."""
    md_text = md_path.read_text()
    fm, _ = md_to_html.parse_frontmatter(md_text)
    body_html = md_to_html.convert(md_text)

    # Inline any local <img src="...">
    body_html = inline_assets.inline_images(body_html, base_dir=md_path.parent)

    template_path = _TEMPLATES_DIR / f"base_{level}.html.jinja"
    if not template_path.exists():
        # Fallback to L2 if requested level doesn't exist yet
        template_path = _TEMPLATES_DIR / "base_l2.html.jinja"
    tmpl_text = template_path.read_text()
    tokens_css = _DESIGN_TOKENS_PATH.read_text()

    rendered = Template(tmpl_text).render(
        title=fm.get("title", md_path.stem),
        design_tokens_css=tokens_css,
        source_md_relpath=md_path.name,
        skill_version=SKILL_VERSION,
        generated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        body_html=body_html,
    )
    return rendered


def build_preview(md_path: Path) -> dict:
    """Build companion HTML for md_path. Returns:
        {
          "mode": "single_file" | "topic_folder",
          "output_path": Path (file for single_file, folder for topic_folder),
          "level": "l1"|"l2"|"l3",
          "taste_check": {...}
        }
    """
    md_path = Path(md_path).resolve()
    use_folder = topic_folder_decider.should_use_topic_folder(md_path)
    level = interactivity_detect.detect_level(md_path)
    html = _render(md_path, level)
    taste = design_taste_check.detect_taboos_static(html)

    if use_folder:
        folder = topic_folder_builder.build_folder(
            md_path, chapters=[], figures=[],
        )
        index_path = folder / "index.html"
        index_path.write_text(html)
        return {
            "mode": "topic_folder",
            "output_path": folder,
            "level": level,
            "taste_check": taste,
        }
    else:
        out = md_path.parent / (md_path.stem + ".preview.html")
        out.write_text(html)
        return {
            "mode": "single_file",
            "output_path": out,
            "level": level,
            "taste_check": taste,
        }


def main(argv: list[str]) -> int:
    import argparse
    parser = argparse.ArgumentParser(prog="html-preview dispatch")
    parser.add_argument("md_path", type=Path)
    parser.add_argument("--rebuild", action="store_true",
                        help="overwrite existing companion (default: skip if exists)")
    parser.add_argument("--polished", action="store_true",
                        help="(D19, opt-in) use frontend-design plugin (Phase 2.5 not implemented yet)")
    args = parser.parse_args(argv[1:])

    if not args.md_path.exists():
        print(f"Source not found: {args.md_path}", file=sys.stderr)
        return 2
    if args.polished:
        print("[html-preview] --polished requires frontend-design plugin invocation; not yet wired in Phase 2 base. "
              "Falling back to MVP.css default.", file=sys.stderr)

    result = build_preview(args.md_path)
    print(f"[html-preview] mode={result['mode']} level={result['level']} → {result['output_path']}")
    print(f"[html-preview] taste_check: verdict={result['taste_check']['verdict']} "
          f"n_taboos={result['taste_check']['n_taboos']}")
    if result["taste_check"]["suggestions"]:
        for s in result["taste_check"]["suggestions"]:
            print(f"  - {s}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
