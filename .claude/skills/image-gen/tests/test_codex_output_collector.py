"""Test that collector finds latest PNG and moves to target."""
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from codex_output_collector import collect_latest  # noqa


def test_finds_and_moves_latest_png():
    with tempfile.TemporaryDirectory() as src_str, tempfile.TemporaryDirectory() as dst_str:
        src = Path(src_str)
        dst = Path(dst_str)
        # Create 2 PNGs with different mtimes
        old = src / "old.png"
        old.write_bytes(b"\x89PNG\r\n\x1a\n" + b"old")
        time.sleep(0.05)
        new = src / "new.png"
        new.write_bytes(b"\x89PNG\r\n\x1a\n" + b"new")

        target_name = "moved_fig1.png"
        result = collect_latest([src], dst, target_name=target_name, since_epoch=0.0)

        assert result is not None, "collect_latest returned None"
        assert result.name == target_name, f"Wrong name: {result.name}"
        assert (dst / target_name).read_bytes().endswith(b"new"), "Wrong PNG content"


def test_returns_none_when_no_pngs():
    with tempfile.TemporaryDirectory() as src_str, tempfile.TemporaryDirectory() as dst_str:
        src = Path(src_str)
        dst = Path(dst_str)
        result = collect_latest([src], dst, target_name="x.png", since_epoch=0.0)
        assert result is None, "Should return None when no PNGs found"


def main():
    tests = [test_finds_and_moves_latest_png, test_returns_none_when_no_pngs]
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
