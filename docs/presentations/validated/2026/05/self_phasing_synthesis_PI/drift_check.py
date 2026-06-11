#!/usr/bin/env python3
"""
drift_check.py — Detect manual edits to preview/*.html vs last build snapshot.

Workflow:
    1. After build_html.py, run: python drift_check.py --snapshot
       → copies preview/*.html → .build_snapshot/
    2. Manually edit preview/*.html in browser/editor (style tweaks).
    3. Run: python drift_check.py --check
       → reports which slide drifted, which field (canvas_html/speaker/tier3),
         and where in build_html.py to write back.
    4. Apply Edit patches to build_html.py based on report.
    5. Re-run build_html.py + drift_check.py --check → expect 0 drift.

Usage:
    python drift_check.py --snapshot          # take baseline (run after build)
    python drift_check.py --check             # detect drift
    python drift_check.py --check --verbose   # show raw diff lines
    python drift_check.py --check --slide 09b_hp33_mechanism  # check one slide
"""
import argparse
import difflib
import re
import shutil
from pathlib import Path

HERE = Path(__file__).parent
PREVIEW = HERE / "preview"
SNAPSHOT = HERE / ".build_snapshot"
BUILD_PY = HERE / "build_html.py"


def normalize_html(text):
    """Normalize HTML for diff: strip whitespace, unify escape, drop comments."""
    text = re.sub(r"<!--.*?-->", "", text, flags=re.DOTALL)  # drop comments
    text = text.replace("&gt;", ">").replace("&lt;", "<")
    text = text.replace("&amp;", "&").replace("&#39;", "'").replace("&quot;", '"')
    lines = []
    for line in text.splitlines():
        line = re.sub(r"\s+", " ", line).strip()
        if line:
            lines.append(line)
    return lines


def classify_drift_field(line):
    """Guess which add() field the drift line belongs to."""
    if 'class="slide-title' in line or "<h1" in line:
        return "title"
    if 'class="en-subtitle"' in line:
        return "en"
    if 'class="speak-text"' in line or "[標題]" in line or "[結論]" in line:
        return "speaker"
    if 'class="tier3"' in line or "[ORAL-OPTIONAL]" in line:
        return "tier3"
    if 'class="slide-canvas"' in line or "</article>" in line:
        return "canvas_html (boundary)"
    return "canvas_html"


def find_py_line(slide_id):
    """Locate the `add(id="{slide_id}"` line in build_html.py."""
    py_text = BUILD_PY.read_text(encoding="utf-8")
    for i, line in enumerate(py_text.splitlines(), 1):
        if f'add(id="{slide_id}"' in line:
            return i
    return -1


def snapshot_cmd():
    """Copy preview/*.html → .build_snapshot/ as baseline."""
    SNAPSHOT.mkdir(exist_ok=True)
    n = 0
    for html in PREVIEW.glob("slide_*.html"):
        shutil.copy2(html, SNAPSHOT / html.name)
        n += 1
    # also snapshot index.html
    idx = PREVIEW / "index.html"
    if idx.exists():
        shutil.copy2(idx, SNAPSHOT / "index.html")
        n += 1
    print(f"[snapshot] {n} files → .build_snapshot/")


def check_cmd(verbose=False, only_slide=None):
    """Detect drift between preview/ and .build_snapshot/."""
    if not SNAPSHOT.exists():
        print("[error] .build_snapshot/ not found. Run --snapshot first.")
        return 1
    drift_count = 0
    for current in sorted(PREVIEW.glob("slide_*.html")):
        slide_id = current.stem.replace("slide_", "")  # e.g. 09b_hp33_mechanism
        if only_slide and slide_id != only_slide:
            continue
        baseline = SNAPSHOT / current.name
        if not baseline.exists():
            print(f"[NEW]   {current.name} — no baseline (new slide?)")
            drift_count += 1
            continue
        cur_lines = normalize_html(current.read_text(encoding="utf-8"))
        base_lines = normalize_html(baseline.read_text(encoding="utf-8"))
        if cur_lines == base_lines:
            continue
        # drift detected — classify + locate
        drift_count += 1
        diff = list(difflib.unified_diff(base_lines, cur_lines, lineterm="", n=0))
        added = [l for l in diff if l.startswith("+") and not l.startswith("+++")]
        removed = [l for l in diff if l.startswith("-") and not l.startswith("---")]
        fields = set()
        for l in added + removed:
            fields.add(classify_drift_field(l[1:]))
        py_line = find_py_line(slide_id)
        py_loc = f"build_html.py:{py_line}" if py_line > 0 else "build_html.py:?"
        print(f"[DRIFT] slide_{slide_id}.html")
        print(f"        +{len(added)} / -{len(removed)} lines (normalized)")
        print(f"        likely field(s): {', '.join(sorted(fields))}")
        print(f"        py location: {py_loc} (add(id=\"{slide_id}\"))")
        if verbose:
            for l in diff[:20]:
                print(f"        {l}")
            if len(diff) > 20:
                print(f"        ... (+{len(diff)-20} more lines)")
        print()
    if drift_count == 0:
        print("[ok] 0 drift — preview/ matches .build_snapshot/")
        return 0
    print(f"[summary] {drift_count} slide(s) drifted")
    return 1


def main():
    ap = argparse.ArgumentParser(description="Drift check for build_html.py outputs")
    ap.add_argument("--snapshot", action="store_true", help="copy preview to baseline")
    ap.add_argument("--check", action="store_true", help="detect drift")
    ap.add_argument("--verbose", action="store_true", help="show raw diff")
    ap.add_argument("--slide", default=None, help="check only one slide id")
    args = ap.parse_args()
    if args.snapshot:
        snapshot_cmd()
    elif args.check:
        return check_cmd(verbose=args.verbose, only_slide=args.slide)
    else:
        ap.print_help()


if __name__ == "__main__":
    raise SystemExit(main() or 0)
