"""Test inline_assets: replace <img src="..."> with data: URIs."""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from inline_assets import inline_images  # noqa
from PIL import Image


def make_test_png(path: Path):
    Image.new("RGB", (16, 16), (50, 100, 200)).save(path, "PNG")


def test_inlines_existing_relative_png():
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        png = d / "fig.png"
        make_test_png(png)
        html = '<img src="fig.png" alt="x">'
        out = inline_images(html, base_dir=d)
        assert "data:image/png;base64," in out, f"Not inlined: {out}"
        assert 'src="fig.png"' not in out, "Original src not replaced"
        assert 'alt="x"' in out, "Other attributes preserved"


def test_keeps_external_url_unchanged():
    with tempfile.TemporaryDirectory() as tmp:
        html = '<img src="https://example.com/img.png">'
        out = inline_images(html, base_dir=Path(tmp))
        assert out == html, "External URL should be untouched"


def test_keeps_data_uri_unchanged():
    with tempfile.TemporaryDirectory() as tmp:
        html = '<img src="data:image/png;base64,abc">'
        out = inline_images(html, base_dir=Path(tmp))
        assert out == html, "data: URI should be untouched"


def test_missing_file_keeps_src_with_warning_comment():
    with tempfile.TemporaryDirectory() as tmp:
        html = '<img src="missing.png">'
        out = inline_images(html, base_dir=Path(tmp))
        assert 'src="missing.png"' in out, "Should keep src on miss"
        assert "<!--" in out and "missing" in out.lower(), "Should warn via comment"


def main():
    tests = [test_inlines_existing_relative_png, test_keeps_external_url_unchanged,
             test_keeps_data_uri_unchanged, test_missing_file_keeps_src_with_warning_comment]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS: {t.__name__}")
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
