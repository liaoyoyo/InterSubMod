#!/usr/bin/env python3
"""Security regression tests for the public SVG/PNG extractor."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools/extract_svg_for_github.py"
SPEC = importlib.util.spec_from_file_location("extract_svg_for_github_for_test", SCRIPT)
if SPEC is None or SPEC.loader is None:  # pragma: no cover - import infrastructure failure
    raise RuntimeError(f"cannot load SVG extractor: {SCRIPT}")
EXTRACTOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = EXTRACTOR
SPEC.loader.exec_module(EXTRACTOR)


class ExtractSvgForGithubSecurityTest(unittest.TestCase):
    def test_fragment_urls_are_allowed(self) -> None:
        svg = """
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 10">
          <defs>
            <linearGradient id="gradient"><stop offset="0" stop-color="#123456"/></linearGradient>
            <marker id="arrow"><path d="M0,0 L4,2 L0,4 z"/></marker>
          </defs>
          <style>.filled { fill: url('#gradient'); }</style>
          <rect class="filled" width="10" height="10" style="stroke:url(#gradient)"/>
          <path d="M10,5 L20,5" marker-end="url(#arrow)"/>
        </svg>
        """
        EXTRACTOR.validate_svg_security(svg)

    def test_css_import_is_rejected(self) -> None:
        for css in (
            "@import url(https://example.invalid/theme.css);",
            "@import 'https://example.invalid/theme.css';",
            "@im/* comment */port url(https://example.invalid/theme.css);",
        ):
            with self.subTest(css=css):
                svg = (
                    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10">'
                    f"<style>{css}</style><rect width=\"10\" height=\"10\"/></svg>"
                )
                with self.assertRaisesRegex(ValueError, "CSS @import"):
                    EXTRACTOR.validate_svg_security(svg)

    def test_non_fragment_css_urls_are_rejected_everywhere(self) -> None:
        cases = {
            "style_element": "<style>.x{fill:url(https://example.invalid/pixel)}</style><rect class='x'/>",
            "style_attribute": "<rect style='fill:url(data:image/png;base64,AAAA)'/>",
            "presentation_attribute": "<rect fill='url(file:///etc/passwd)'/>",
        }
        for location, body in cases.items():
            with self.subTest(location=location):
                svg = (
                    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10">'
                    f"{body}</svg>"
                )
                with self.assertRaisesRegex(ValueError, "non-fragment CSS url"):
                    EXTRACTOR.validate_svg_security(svg)

    def test_renderer_fails_if_static_validation_is_bypassed(self) -> None:
        svg = """<?xml version="1.0"?>
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10">
          <style>@import url(https://example.invalid/theme.css);</style>
          <rect width="10" height="10" fill="#123456"/>
        </svg>
        """
        with tempfile.TemporaryDirectory(prefix="ism-svg-request-test-") as temporary:
            output = Path(temporary) / "should-not-exist.png"
            with self.assertRaisesRegex(ValueError, "disallowed external request"):
                EXTRACTOR.render_png(svg, output)
            self.assertFalse(output.exists())

    def test_renderer_uses_default_chromium_sandbox_and_renders_local_svg(self) -> None:
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn("--no-sandbox", source)
        svg = """<?xml version="1.0"?>
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 5">
          <defs><linearGradient id="g"><stop stop-color="#123456"/></linearGradient></defs>
          <rect width="10" height="5" fill="url(#g)"/>
        </svg>
        """
        with tempfile.TemporaryDirectory(prefix="ism-svg-sandbox-test-") as temporary:
            output = Path(temporary) / "local.png"
            self.assertEqual(EXTRACTOR.render_png(svg, output), (700, 350))
            self.assertTrue(output.is_file())


if __name__ == "__main__":
    unittest.main()
