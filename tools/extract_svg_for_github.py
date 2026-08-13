#!/usr/bin/env python3
"""Extract selected inline SVG figures from the public explanation pages.

The default source and output directories are resolved from this script's Git
checkout, never from the caller's current directory or a hard-coded machine
path.  ``--source-dir`` and ``--output-dir`` accept either repo-relative or
absolute paths.  Rendering happens in an automatically cleaned temporary
directory under the selected output directory and each result is replaced
atomically only after SVG and PNG validation succeeds.

Examples:

    python3 tools/extract_svg_for_github.py --names architecture-overview
    python3 tools/extract_svg_for_github.py \
        --source-dir docs/explain --output-dir docs/images \
        --names architecture-overview methylation-circularity
    python3 tools/extract_svg_for_github.py --verify-only --names funnel-7samples
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE_DIR = REPO_ROOT / "docs" / "explain"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "docs" / "images"
BACKGROUND = "#FAF9F5"
SVG_RE = re.compile(r"<svg\b.*?</svg>", re.DOTALL | re.IGNORECASE)
CSS_COMMENT_RE = re.compile(r"/\*.*?\*/", re.DOTALL)
CSS_URL_RE = re.compile(
    r"url\s*\(\s*(?:(?P<quote>['\"])(?P<quoted>.*?)"
    r"(?P=quote)|(?P<unquoted>[^)]*?))\s*\)",
    re.IGNORECASE | re.DOTALL,
)
ALLOWED_RENDER_REQUEST_SCHEMES = {"about", "blob", "data", "file"}


@dataclass(frozen=True)
class Figure:
    source_file: str
    svg_index: int
    output_name: str
    caption: str


FIGURES = (
    Figure(
        "11_system-map-overview.standalone.html",
        1,
        "architecture-overview",
        "五層系統全景：資料 → 上游工具 → 兩支 C++ → Python → HTML",
    ),
    Figure(
        "11_system-map-overview.standalone.html",
        2,
        "methylation-circularity",
        "為什麼甲基化不能用來重建譜系（cis-ASM 循環）",
    ),
    Figure(
        "11_system-map-overview.standalone.html",
        3,
        "funnel-7samples",
        "7 technical datasets／6 biological IDs funnel：候選位點依可排序性與圖形狀態分層",
    ),
    Figure(
        "12_intersubmod-io.standalone.html",
        1,
        "ism-internal-pipeline",
        "InterSubMod 內部八階段處理鏈",
    ),
    Figure(
        "13_longlineage-io.standalone.html",
        1,
        "longlineage-funnel",
        "LongLineage 位點流失漏斗（目前不形成已確認細胞譜系）",
    ),
    Figure(
        "14_upstream-data.standalone.html",
        1,
        "upstream-toolchain",
        "上游前處理鏈與 sidecar 串流設計",
    ),
    Figure(
        "15_python-html-layer.standalone.html",
        1,
        "workstation-refuse-design",
        "工作站生成器的拒絕渲染設計（防止規格所列必填欄位被靜默省略）",
    ),
    Figure(
        "16_how-to-run.standalone.html",
        1,
        "howto-six-steps",
        "操作六步驟與各自的驗收點",
    ),
)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_from_repo(value: str | Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = REPO_ROOT / path
    return path.resolve()


def display_path(path: Path) -> str:
    try:
        return path.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path)


def parse_view_box(svg: str) -> tuple[float, float, float, float]:
    root = ET.fromstring(svg)
    raw = root.attrib.get("viewBox")
    if raw is None:
        raise ValueError("SVG root has no viewBox; aspect ratio cannot be verified")
    values = [float(value) for value in re.split(r"[\s,]+", raw.strip()) if value]
    if len(values) != 4 or values[2] <= 0 or values[3] <= 0:
        raise ValueError(f"invalid SVG viewBox: {raw!r}")
    return values[0], values[1], values[2], values[3]


def validate_css_resource_references(css: str, location: str) -> None:
    """Allow local SVG fragment references but reject imported resources."""
    without_comments = CSS_COMMENT_RE.sub("", css)
    if re.search(r"@\s*import\b", without_comments, flags=re.IGNORECASE):
        raise ValueError(f"unsafe CSS @import in {location}")

    url_openings = list(re.finditer(r"url\s*\(", without_comments, flags=re.IGNORECASE))
    url_matches = list(CSS_URL_RE.finditer(without_comments))
    if len(url_matches) != len(url_openings):
        raise ValueError(f"malformed or unsupported CSS url() in {location}")

    for match in url_matches:
        value = (match.group("quoted") if match.group("quote") else match.group("unquoted") or "").strip()
        if re.fullmatch(r"#[A-Za-z_][A-Za-z0-9_.:-]*", value) is None:
            raise ValueError(
                f"unsafe non-fragment CSS url() in {location}: {value!r}"
            )


def validate_svg_security(svg: str) -> None:
    """Reject active or remotely loaded SVG content before browser rendering."""
    root = ET.fromstring(svg)
    for element in root.iter():
        local_name = element.tag.rsplit("}", 1)[-1].lower()
        if local_name in {"script", "foreignobject"}:
            raise ValueError(f"unsafe SVG element: {local_name}")
        if local_name == "style":
            validate_css_resource_references(
                "".join(element.itertext()), "<style> element"
            )
        for raw_name, raw_value in element.attrib.items():
            attr_name = raw_name.rsplit("}", 1)[-1].lower()
            value = raw_value.strip().lower()
            if attr_name.startswith("on"):
                raise ValueError(f"unsafe SVG event attribute: {attr_name}")
            if attr_name == "href" and value and not value.startswith("#"):
                raise ValueError(f"unsafe non-fragment SVG href: {raw_value!r}")
            if attr_name == "style" or re.search(
                r"url\s*\(", raw_value, flags=re.IGNORECASE
            ):
                validate_css_resource_references(
                    raw_value, f"{local_name}[{attr_name}]"
                )


def add_background(svg: str) -> str:
    """Insert an opaque background while retaining the source viewBox."""
    x, y, width, height = parse_view_box(svg)
    rect = (
        f'<rect x="{x:g}" y="{y:g}" width="{width:g}" height="{height:g}" '
        f'fill="{BACKGROUND}"/>'
    )
    for anchor in ("</desc>", "</title>"):
        if anchor in svg:
            offset = svg.index(anchor) + len(anchor)
            return f"{svg[:offset]}\n  {rect}{svg[offset:]}"
    offset = svg.index(">") + 1
    return f"{svg[:offset]}\n  {rect}{svg[offset:]}"


def extract_inline_svg(source_path: Path, svg_index: int) -> tuple[str, bytes]:
    source_bytes = source_path.read_bytes()
    blocks = SVG_RE.findall(source_bytes.decode("utf-8"))
    if len(blocks) < svg_index:
        raise ValueError(
            f"{display_path(source_path)} contains {len(blocks)} SVG blocks; "
            f"cannot select #{svg_index}"
        )
    return blocks[svg_index - 1], source_bytes


def build_svg_document(figure: Figure, source_path: Path) -> tuple[str, dict[str, str]]:
    inline_svg, source_bytes = extract_inline_svg(source_path, figure.svg_index)
    validate_svg_security(inline_svg)
    svg = add_background(inline_svg)
    opening_tag = svg.split(">", 1)[0]
    if "xmlns=" not in opening_tag:
        svg = svg.replace("<svg", '<svg xmlns="http://www.w3.org/2000/svg"', 1)
    source_sha256 = sha256_bytes(source_bytes)
    inline_sha256 = sha256_bytes(inline_svg.encode("utf-8"))
    source_label = display_path(source_path)
    document = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f"<!-- {figure.caption}\n"
        f"     source: {source_label} (svg #{figure.svg_index})\n"
        f"     source_document_sha256: {source_sha256}\n"
        f"     inline_svg_sha256: {inline_sha256}\n"
        "     generated_by: tools/extract_svg_for_github.py; do not edit directly. -->\n"
        f"{svg}\n"
    )
    ET.fromstring(document)
    return document, {
        "source": source_label,
        "source_document_sha256": source_sha256,
        "inline_svg_sha256": inline_sha256,
    }


def render_png(svg_document: str, png_path: Path) -> tuple[int, int]:
    """Render only the SVG element at the exact viewBox aspect ratio."""
    from PIL import Image
    from playwright.sync_api import sync_playwright

    _, _, view_width, view_height = parse_view_box(svg_document)
    pixel_width = min(max(round(view_width), 700), 1100)
    pixel_height = round(pixel_width * view_height / view_width)
    inline_svg = re.sub(r"<\?xml[^>]*\?>\s*", "", svg_document, count=1)
    markup = (
        '<!doctype html><meta charset="utf-8">'
        "<style>html,body{margin:0;padding:0;background:"
        f"{BACKGROUND};width:{pixel_width}px;height:{pixel_height}px;overflow:hidden}}"
        f"svg{{display:block;width:{pixel_width}px!important;height:{pixel_height}px!important}}</style>"
        f"{inline_svg}"
    )
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(
                viewport={"width": pixel_width, "height": pixel_height},
                device_scale_factor=1,
            )
            disallowed_requests: list[str] = []

            def record_request(request: object) -> None:
                url = str(getattr(request, "url", ""))
                scheme = url.split(":", 1)[0].lower() if ":" in url else ""
                if scheme not in ALLOWED_RENDER_REQUEST_SCHEMES:
                    disallowed_requests.append(url)

            def block_disallowed_route(route: object) -> None:
                request = getattr(route, "request")
                url = str(getattr(request, "url", ""))
                scheme = url.split(":", 1)[0].lower() if ":" in url else ""
                if scheme not in ALLOWED_RENDER_REQUEST_SCHEMES:
                    getattr(route, "abort")("blockedbyclient")
                else:
                    getattr(route, "continue_")()

            page.on("request", record_request)
            page.route("**/*", block_disallowed_route)
            page.set_content(markup, wait_until="load")
            page.evaluate("document.fonts.ready")
            if disallowed_requests:
                attempted = ", ".join(sorted(set(disallowed_requests))[:5])
                raise ValueError(
                    "SVG rendering attempted a disallowed external request: "
                    f"{attempted}"
                )
            locator = page.locator("svg")
            box = locator.bounding_box()
            if box is None or round(box["width"]) != pixel_width or round(box["height"]) != pixel_height:
                raise ValueError(f"rendered SVG bounds do not preserve viewBox ratio: {box}")
            locator.screenshot(path=str(png_path), animations="disabled")
        finally:
            browser.close()
    with Image.open(png_path) as image:
        image.verify()
    with Image.open(png_path) as image:
        if image.format != "PNG" or image.size != (pixel_width, pixel_height):
            raise ValueError(
                f"invalid PNG output: format={image.format}, size={image.size}, "
                f"expected={(pixel_width, pixel_height)}"
            )
    return pixel_width, pixel_height


def selected_figures(names: Iterable[str] | None) -> list[Figure]:
    by_name = {figure.output_name: figure for figure in FIGURES}
    if not names:
        return list(FIGURES)
    unknown = sorted(set(names) - set(by_name))
    if unknown:
        raise ValueError(f"unknown figure name(s): {', '.join(unknown)}")
    return [by_name[name] for name in names]


def process_figure(
    figure: Figure,
    source_dir: Path,
    output_dir: Path,
    temp_dir: Path,
    verify_only: bool,
) -> dict[str, object]:
    source_path = source_dir / figure.source_file
    svg_path = output_dir / f"{figure.output_name}.svg"
    png_path = output_dir / f"{figure.output_name}.png"
    if not source_path.is_file():
        raise FileNotFoundError(f"source HTML not found: {source_path}")

    svg_document, provenance = build_svg_document(figure, source_path)
    temporary_svg = temp_dir / svg_path.name
    temporary_png = temp_dir / png_path.name
    temporary_svg.write_text(svg_document, encoding="utf-8")
    dimensions = render_png(svg_document, temporary_png)
    expected_svg_sha256 = sha256_file(temporary_svg)
    expected_png_sha256 = sha256_file(temporary_png)

    if verify_only:
        if not svg_path.is_file() or not png_path.is_file():
            raise FileNotFoundError(f"generated pair missing: {svg_path}, {png_path}")
        actual_svg_sha256 = sha256_file(svg_path)
        actual_png_sha256 = sha256_file(png_path)
        if actual_svg_sha256 != expected_svg_sha256 or actual_png_sha256 != expected_png_sha256:
            raise ValueError(
                f"generated pair is stale for {figure.output_name}: "
                f"svg={actual_svg_sha256 == expected_svg_sha256}, "
                f"png={actual_png_sha256 == expected_png_sha256}"
            )
    else:
        os.replace(temporary_svg, svg_path)
        os.replace(temporary_png, png_path)
        actual_svg_sha256 = sha256_file(svg_path)
        actual_png_sha256 = sha256_file(png_path)

    return {
        "name": figure.output_name,
        **provenance,
        "svg_path": display_path(svg_path),
        "svg_sha256": actual_svg_sha256,
        "png_path": display_path(png_path),
        "png_sha256": actual_png_sha256,
        "png_dimensions": list(dimensions),
        "status": "VERIFIED" if verify_only else "GENERATED_AND_VERIFIED",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-dir",
        default=str(DEFAULT_SOURCE_DIR),
        help="source HTML directory; relative paths resolve from the current checkout",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="SVG/PNG output directory; relative paths resolve from the current checkout",
    )
    parser.add_argument(
        "--names",
        nargs="+",
        choices=[figure.output_name for figure in FIGURES],
        help="one or more output names (default: all configured figures)",
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="regenerate in a temporary directory and compare hashes without changing outputs",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source_dir = resolve_from_repo(args.source_dir)
    output_dir = resolve_from_repo(args.output_dir)
    if not source_dir.is_dir():
        raise FileNotFoundError(f"source directory not found: {source_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    figures = selected_figures(args.names)
    records: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix=".extract-svg-", dir=output_dir) as temporary:
        temp_dir = Path(temporary)
        for figure in figures:
            records.append(
                process_figure(figure, source_dir, output_dir, temp_dir, args.verify_only)
            )
    print(
        json.dumps(
            {
                "source_dir": display_path(source_dir),
                "output_dir": display_path(output_dir),
                "mode": "VERIFY_ONLY" if args.verify_only else "GENERATE",
                "count": len(records),
                "figures": records,
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
