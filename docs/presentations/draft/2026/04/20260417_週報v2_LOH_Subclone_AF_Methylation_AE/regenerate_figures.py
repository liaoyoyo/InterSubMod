"""Regenerate all 5 key figures for Weekly Report v2.

Usage:
    python3 regenerate_figures.py

Output: figures/01_*.png through 05_*.png at 300 DPI.
"""
import os
import sys
from figure_specs import FIGURE_FUNCS


def main():
    print("=" * 60)
    print("Regenerating figures for Weekly Report v2")
    print("=" * 60)
    for tag, func in FIGURE_FUNCS:
        print(f"  Building figure {tag} ({func.__name__}) ...", end=" ")
        try:
            path = func()
            size_kb = os.path.getsize(path) / 1024
            print(f"OK ({size_kb:.0f} KB)")
            print(f"    → {os.path.relpath(path)}")
        except Exception as e:
            print(f"FAIL: {e}")
            return 1
    print("=" * 60)
    print(f"All {len(FIGURE_FUNCS)} figures generated at 300 DPI.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
