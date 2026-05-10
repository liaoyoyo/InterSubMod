"""Program-track image rendering via pycairo.

Supports: flow_diagram (nodes + edges), icon (shape + glyph).
Deterministic, no AI randomness, no API quota. Fast.
"""
from __future__ import annotations

import math
from pathlib import Path

import cairo
import yaml


def _parse_size(size_str: str) -> tuple[int, int]:
    w, h = size_str.lower().split("x")
    return int(w), int(h)


def _parse_hex(c: str) -> tuple[float, float, float, float]:
    """Hex color #RRGGBB or transparent → RGBA tuple in 0-1 range."""
    if c == "transparent":
        return (0.0, 0.0, 0.0, 0.0)
    c = c.lstrip("#")
    r, g, b = int(c[0:2], 16), int(c[2:4], 16), int(c[4:6], 16)
    return (r / 255, g / 255, b / 255, 1.0)


def _draw_flow(ctx: cairo.Context, data: dict, w: int, h: int) -> None:
    style = data.get("style", {})
    palette = style.get("palette") or ["#1f2937", "#3b82f6"]
    text_rgba = _parse_hex(palette[0])
    accent_rgba = _parse_hex(palette[1])
    bg_rgba = _parse_hex(style.get("background", "white") if style.get("background") != "white" else "#ffffff")
    font_size = float(style.get("font_size", 14))
    font_family = style.get("font_family", "DejaVu Sans")

    # Background
    ctx.set_source_rgba(*bg_rgba)
    ctx.rectangle(0, 0, w, h)
    ctx.fill()

    ctx.select_font_face(font_family, cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.set_font_size(font_size)

    nodes = data.get("nodes", [])
    node_by_id = {n["id"]: n for n in nodes}

    # Draw edges first (under nodes)
    for e in data.get("edges", []):
        a = node_by_id.get(e["from"])
        b = node_by_id.get(e["to"])
        if not a or not b:
            continue
        ctx.set_source_rgba(*accent_rgba)
        ctx.set_line_width(2.0)
        ctx.move_to(a["x"], a["y"])
        ctx.line_to(b["x"], b["y"])
        ctx.stroke()
        # Arrow head
        angle = math.atan2(b["y"] - a["y"], b["x"] - a["x"])
        ah = 10.0
        ctx.move_to(b["x"], b["y"])
        ctx.line_to(b["x"] - ah * math.cos(angle - math.pi / 7),
                    b["y"] - ah * math.sin(angle - math.pi / 7))
        ctx.line_to(b["x"] - ah * math.cos(angle + math.pi / 7),
                    b["y"] - ah * math.sin(angle + math.pi / 7))
        ctx.close_path()
        ctx.fill()

    # Draw nodes
    for n in nodes:
        label = n.get("label", "")
        ext = ctx.text_extents(label)
        pad_x, pad_y = 16.0, 12.0
        box_w = ext.width + 2 * pad_x
        box_h = ext.height + 2 * pad_y
        x0 = n["x"] - box_w / 2
        y0 = n["y"] - box_h / 2
        # Node bg
        ctx.set_source_rgba(*bg_rgba)
        ctx.rectangle(x0, y0, box_w, box_h)
        ctx.fill()
        # Node border
        ctx.set_source_rgba(*text_rgba)
        ctx.set_line_width(1.5)
        ctx.rectangle(x0, y0, box_w, box_h)
        ctx.stroke()
        # Label
        ctx.move_to(n["x"] - ext.width / 2 - ext.x_bearing,
                    n["y"] + ext.height / 2 - (ext.height + ext.y_bearing))
        ctx.show_text(label)


def _draw_icon(ctx: cairo.Context, data: dict, w: int, h: int) -> None:
    style = data.get("style", {})
    fill_rgba = _parse_hex(data.get("fill_color", "#10b981"))
    glyph_rgba = _parse_hex(data.get("glyph_color", "#ffffff"))
    bg = style.get("background", "transparent")
    bg_rgba = _parse_hex(bg) if bg != "white" else (1, 1, 1, 1)
    pad = float(style.get("padding_px", 24))
    shape = data.get("shape", "circle")
    glyph = data.get("glyph", "?")

    # Background
    ctx.set_source_rgba(*bg_rgba)
    ctx.rectangle(0, 0, w, h)
    ctx.fill()

    # Shape
    ctx.set_source_rgba(*fill_rgba)
    cx, cy = w / 2, h / 2
    radius = (min(w, h) - 2 * pad) / 2
    if shape == "circle":
        ctx.arc(cx, cy, radius, 0, 2 * math.pi)
        ctx.fill()
    elif shape in ("square", "rounded_square"):
        x0, y0 = pad, pad
        side = min(w, h) - 2 * pad
        if shape == "rounded_square":
            r = side * 0.15
            ctx.move_to(x0 + r, y0)
            ctx.line_to(x0 + side - r, y0)
            ctx.arc(x0 + side - r, y0 + r, r, -math.pi / 2, 0)
            ctx.line_to(x0 + side, y0 + side - r)
            ctx.arc(x0 + side - r, y0 + side - r, r, 0, math.pi / 2)
            ctx.line_to(x0 + r, y0 + side)
            ctx.arc(x0 + r, y0 + side - r, r, math.pi / 2, math.pi)
            ctx.line_to(x0, y0 + r)
            ctx.arc(x0 + r, y0 + r, r, math.pi, 3 * math.pi / 2)
            ctx.close_path()
            ctx.fill()
        else:
            ctx.rectangle(x0, y0, side, side)
            ctx.fill()
    else:
        # Fallback: circle
        ctx.arc(cx, cy, radius, 0, 2 * math.pi)
        ctx.fill()

    # Glyph
    ctx.set_source_rgba(*glyph_rgba)
    font_size = radius * 1.3
    ctx.select_font_face("DejaVu Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    ctx.set_font_size(font_size)
    ext = ctx.text_extents(glyph)
    ctx.move_to(cx - ext.width / 2 - ext.x_bearing,
                cy + ext.height / 2 - (ext.height + ext.y_bearing))
    ctx.show_text(glyph)


def render(yaml_path: Path, output: Path) -> None:
    """Render YAML prompt to PNG.

    Raises ValueError if backend is not 'cairo'.
    """
    data = yaml.safe_load(Path(yaml_path).read_text())
    if data.get("backend") != "cairo":
        raise ValueError(
            f"cairo_render only handles backend='cairo'; got '{data.get('backend')}' "
            f"in {yaml_path}"
        )
    out = data.get("output", {})
    w, h = _parse_size(out.get("size", "1024x768"))
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, w, h)
    ctx = cairo.Context(surface)

    t = data.get("type")
    if t == "flow_diagram":
        _draw_flow(ctx, data, w, h)
    elif t == "icon":
        _draw_icon(ctx, data, w, h)
    else:
        raise ValueError(f"cairo_render does not support type='{t}'")

    surface.write_to_png(str(output))


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 3:
        print("Usage: cairo_render.py <prompt.yaml> <output.png>")
        return 2
    render(Path(argv[1]), Path(argv[2]))
    print(f"[cairo] rendered → {argv[2]}")
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
