"""Test pil_postprocess.resize_to_target and to_base64."""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from pil_postprocess import resize_to_target, to_base64  # noqa
from PIL import Image


def make_test_png(path: Path, w: int, h: int) -> None:
    Image.new("RGB", (w, h), (200, 200, 200)).save(path, "PNG")


def test_resize_preserves_aspect_or_pads():
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "in.png"
        dst = Path(tmp) / "out.png"
        make_test_png(src, 800, 400)  # 2:1
        resize_to_target(src, dst, target_size=(1024, 1024))
        img = Image.open(dst)
        assert img.size == (1024, 1024), f"Expected 1024x1024, got {img.size}"


def test_to_base64_starts_with_data_uri():
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "in.png"
        make_test_png(src, 100, 100)
        b64 = to_base64(src)
        assert b64.startswith("data:image/png;base64,"), f"Bad prefix: {b64[:40]}"
        assert len(b64) > 100, "Base64 too short"


def main():
    tests = [test_resize_preserves_aspect_or_pads, test_to_base64_starts_with_data_uri]
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
