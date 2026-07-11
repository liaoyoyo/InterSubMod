#!/usr/bin/env python3
"""Embed the precompiled Recharts report runtime into a self-contained HTML artifact."""

from __future__ import annotations

import argparse
import base64
import gzip
import json
import re
import textwrap
from pathlib import Path

RUNTIME_MARKER = "<!-- DATA_ANALYTICS_HTML_REPORT_RUNTIME -->"
RUNTIME_ASSET = "html-report-runtime.html"
RUNTIME_SOURCE_ID = "data-analytics-html-report-runtime-source"
LOCAL_RUNTIME_SOURCE = "html-report-runtime-source.v0.2.6.txt"
LOCAL_RUNTIME_STYLES = "html-report-runtime-styles.v0.2.6.css"
CHART_RUNTIME_STYLES = """
[data-recharts-chart] [data-recharts-live] { display: none !important; }
[data-recharts-chart][data-recharts-ready="true"] [data-recharts-fallback] { display: none !important; }
[data-recharts-chart][data-recharts-ready="true"] [data-recharts-live] { display: block !important; }
[data-recharts-chart] [data-recharts-live] svg {
  display: block !important;
  height: 100% !important;
  max-width: 100% !important;
  min-width: 0 !important;
  width: 100% !important;
}
""".strip()
SOURCE_TOOLTIP_SCRIPT = "html-report-source-tooltips.v0.2.6.js"


def plugin_root() -> Path:
    return Path(__file__).resolve().parent


def read_runtime_html(runtime_html: Path | None) -> str:
    if runtime_html is not None:
        return runtime_html.read_text(encoding="utf-8")

    local_source = plugin_root() / LOCAL_RUNTIME_SOURCE
    local_styles = plugin_root() / LOCAL_RUNTIME_STYLES
    if local_source.is_file() and local_styles.is_file():
        encoded = local_source.read_text(encoding="utf-8").strip()
        scripts = gzip.decompress(base64.b64decode(encoded)).decode("utf-8")
        styles = local_styles.read_text(encoding="utf-8")
        return f'<style>{styles}</style>\n<script type="module">{scripts}</script>\n'

    asset_dir = plugin_root() / "assets"
    parts = sorted(asset_dir.glob(f"{RUNTIME_ASSET}.gz.b64.part*"))
    if not parts:
        raise ValueError(
            "HTML report runtime bundle parts are missing; use the packaged plugin or run npm run build."
        )
    encoded = "".join(part.read_text(encoding="utf-8").strip() for part in parts)
    return gzip.decompress(base64.b64decode(encoded)).decode("utf-8")


def read_source_tooltip_script() -> str:
    return (plugin_root() / SOURCE_TOOLTIP_SCRIPT).read_text(encoding="utf-8")


def runtime_fragments(runtime_html: str) -> tuple[str, str]:
    if "Redirecting to the local widget source" in runtime_html:
        raise ValueError("HTML report runtime asset is still a development redirect; run npm run build.")

    styles = "\n".join(
        match.group(1)
        for match in re.finditer(r"<style[^>]*>(.*?)</style>", runtime_html, flags=re.DOTALL | re.IGNORECASE)
    )
    scripts = "\n".join(
        match.group(1)
        for match in re.finditer(
            r"<script[^>]*type=[\"']module[\"'][^>]*>(.*?)</script>",
            runtime_html,
            flags=re.DOTALL | re.IGNORECASE,
        )
    )
    if not scripts:
        raise ValueError("HTML report runtime bundle does not contain an inline module script.")
    return styles, scripts


def safe_json_script_text(payload: object) -> str:
    text = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    return (
        text.replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("\u2028", "\\u2028")
        .replace("\u2029", "\\u2029")
    )


def compressed_runtime_source(scripts: str) -> str:
    """Keep the large runtime inert until the readable fallback has painted.

    ChatGPT's HTML source view pays a high startup cost when every report
    contains the raw minified React/Recharts bundle. Gzip keeps the exact
    precompiled runtime self-contained while avoiding giant JavaScript lines
    in the authored HTML.
    """

    encoded = base64.b64encode(gzip.compress(scripts.encode("utf-8"), mtime=0)).decode("ascii")
    return textwrap.fill(encoded, width=76)


def runtime_loader_script() -> str:
    """Return a tiny replayable loader for the compressed runtime.

    ChatGPT can serialize a rendered desktop preview and later reopen that
    DOM at mobile width. Keep the compact source and loader in the document
    so the serialized preview can rerun the runtime and repair stale SVG
    geometry; the large decompressed module is still removed after it runs.
    """

    return f"""(() => {{
  const source = document.getElementById("{RUNTIME_SOURCE_ID}");
  if (!(source instanceof HTMLTemplateElement) || typeof DecompressionStream !== "function") {{
    return;
  }}

  function sizeIsStale(current, serialized) {{
    if (!Number.isFinite(current) || current <= 0 || !Number.isFinite(serialized) || serialized <= 0) {{
      return true;
    }}
    return Math.abs(current - serialized) > Math.max(4, Math.max(current, serialized) * 0.08);
  }}

  function shouldBoot() {{
    for (const host of document.querySelectorAll("[data-recharts-chart]")) {{
      if (host.getAttribute("data-recharts-ready") !== "true") return true;
      const mount = host.querySelector("[data-recharts-live]");
      if (!(mount instanceof HTMLElement)) return true;
      const svg = mount.querySelector("svg");
      if (!svg) {{
        if (mount.querySelector(".leaderboard, .heatmap-grid-panel")) continue;
        return true;
      }}
      const viewBox = String(svg.getAttribute("viewBox") ?? "").trim().split(/[\\s,]+/).map(Number);
      if (viewBox.length !== 4) return true;
      const rect = mount.getBoundingClientRect();
      if (sizeIsStale(rect.width, viewBox[2]) || sizeIsStale(rect.height, viewBox[3])) return true;
    }}
    return false;
  }}

  async function boot() {{
    try {{
      const encoded = source.content.textContent.replace(/\\s/g, "");
      const binary = atob(encoded);
      const bytes = new Uint8Array(binary.length);
      for (let index = 0; index < binary.length; index += 1) {{
        bytes[index] = binary.charCodeAt(index);
      }}
      const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
      const runtimeSource = await new Response(stream).text();
      const runtime = document.createElement("script");
      runtime.dataset.dataAnalyticsHtmlReportRuntime = "true";
      runtime.textContent = runtimeSource;
      document.head.append(runtime);
      runtime.remove();
    }} catch {{}}
  }}

  let started = false;
  function start() {{
    if (started) return;
    started = true;
    if (!shouldBoot()) return;
    void boot();
  }}

  function schedule() {{
    if (typeof window.requestIdleCallback === "function") {{
      window.requestIdleCallback(start, {{ timeout: 1500 }});
    }}
    window.setTimeout(start, 1500);
  }}

  if (document.readyState === "complete") {{
    schedule();
  }} else {{
    window.addEventListener("load", schedule, {{ once: true }});
  }}
}})();"""


def chart_ids(payload: object) -> list[str]:
    if not isinstance(payload, dict):
        raise ValueError("Payload must be a JSON object.")
    charts = payload.get("charts")
    if not isinstance(charts, list):
        raise ValueError("Payload must contain a charts array.")
    ids: list[str] = []
    for index, chart in enumerate(charts):
        if not isinstance(chart, dict) or not str(chart.get("id", "")).strip():
            raise ValueError(f"charts[{index}] must have a non-empty id.")
        chart_id = str(chart["id"]).strip()
        if chart_id in ids:
            raise ValueError(f"Duplicate chart id: {chart_id}")
        ids.append(chart_id)
    return ids


def validate_shell(shell: str, ids: list[str]) -> None:
    if shell.count(RUNTIME_MARKER) != 1:
        raise ValueError(f"HTML shell must contain exactly one {RUNTIME_MARKER} marker.")
    marker_start = shell.index(RUNTIME_MARKER)
    host_pattern = re.compile(
        r"\bdata-recharts-chart\s*=\s*([\"'])(.*?)\1",
        flags=re.IGNORECASE,
    )
    hosts = [
        (match.start(), match.group(2))
        for match in host_pattern.finditer(shell[:marker_start])
    ]
    if host_pattern.search(shell[marker_start:]):
        raise ValueError("HTML shell chart hosts must appear before the runtime marker.")

    starts: list[tuple[int, str]] = []
    for chart_id in ids:
        matching_hosts = [start for start, host_id in hosts if host_id == chart_id]
        if not matching_hosts:
            raise ValueError(f'HTML shell is missing data-recharts-chart="{chart_id}".')
        if len(matching_hosts) != 1:
            raise ValueError(f"HTML shell must contain exactly one chart host for {chart_id}.")
        starts.append((matching_hosts[0], chart_id))

    unexpected_ids = sorted({host_id for _, host_id in hosts if host_id not in ids})
    if unexpected_ids:
        raise ValueError(
            "HTML shell has chart hosts without payload entries: " + ", ".join(unexpected_ids)
        )

    starts.sort()
    for index, (start, chart_id) in enumerate(starts):
        end = starts[index + 1][0] if index + 1 < len(starts) else marker_start
        scope = shell[start:end]
        if len(re.findall(r"\bdata-recharts-live\b", scope, flags=re.IGNORECASE)) != 1:
            raise ValueError(f"Chart {chart_id} needs its own data-recharts-live mount.")
        if len(re.findall(r"\bdata-recharts-fallback\b", scope, flags=re.IGNORECASE)) != 1:
            raise ValueError(f"Chart {chart_id} needs its own data-recharts-fallback content.")
        fallback = re.search(
            r"<[^>]*\bdata-recharts-fallback\b[^>]*>(.*?)(?=<[^>]*\bdata-recharts-live\b)",
            scope,
            flags=re.DOTALL | re.IGNORECASE,
        )
        if not fallback:
            raise ValueError(f"Chart {chart_id} needs its own data-recharts-fallback content.")
        fallback_body = re.sub(r"<!--.*?-->", "", fallback.group(1), flags=re.DOTALL)
        if not re.search(
            r"<(?:path|line|polyline|polygon|rect|circle|ellipse|text|table|tr|td|th)\b",
            fallback_body,
            flags=re.IGNORECASE,
        ):
            raise ValueError(f"Chart {chart_id} needs a readable SVG or table fallback.")


def embed(shell: str, payload: object, runtime_html: str) -> str:
    ids = chart_ids(payload)
    validate_shell(shell, ids)
    styles, scripts = runtime_fragments(runtime_html)
    embedded_styles = "\n".join(part for part in (styles, CHART_RUNTIME_STYLES) if part)
    fragments = [
        '<script data-data-analytics-html-report-source-tooltips="true">'
        + read_source_tooltip_script().replace("</script", "<\\/script")
        + "</script>",
        f'<script id="analytics-payload" type="application/json">{safe_json_script_text(payload)}</script>',
        f'<style data-data-analytics-html-report-runtime="true">{embedded_styles}</style>',
    ]
    fragments.append(
        f'<template id="{RUNTIME_SOURCE_ID}" '
        f'data-data-analytics-html-report-runtime-source="true">\n'
        f"{compressed_runtime_source(scripts)}\n"
        "</template>"
    )
    fragments.append(
        '<script data-data-analytics-html-report-runtime-loader="true">'
        f"{runtime_loader_script()}"
        "</script>"
    )
    return shell.replace(RUNTIME_MARKER, "\n".join(fragments))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Authored HTML shell with runtime marker.")
    parser.add_argument("--payload", required=True, type=Path, help="JSON payload with a charts array.")
    parser.add_argument("--output", required=True, type=Path, help="Self-contained HTML output path.")
    parser.add_argument(
        "--runtime-html",
        type=Path,
        help="Optional uncompressed runtime HTML for local tests; packaged runs use bundled parts.",
    )
    args = parser.parse_args()

    shell = args.input.read_text(encoding="utf-8")
    payload = json.loads(args.payload.read_text(encoding="utf-8"))
    output = embed(shell, payload, read_runtime_html(args.runtime_html))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8")
    print(f"Wrote {args.output} with {len(chart_ids(payload))} Recharts chart(s).")


if __name__ == "__main__":
    main()
