#!/usr/bin/env python3
"""fill_report.py — Layer A: by-construction anti-fabrication report renderer.

Industry pattern: literate / reproducible reporting (Quarto / RMarkdown inline code).
The number in the report IS the computed value pulled from a data file at render
time, so a fabricated number is *physically impossible* — you cannot type a value
that is not a key in the data file.

Usage:
    python3 scripts/fill_report.py <template> <data.json> -o <out>
    python3 scripts/fill_report.py report.tmpl.html metrics.json -o report.html

Template placeholders use {{dotted.key}} resolved against data.json (nested dicts
and list indices supported, e.g. {{regions.0.auc}}).

Contract (the whole point — fail loud, never silently fabricate):
  - Every {{key}} MUST resolve to a value in data.json.
  - If ANY placeholder is unresolved -> ERROR, refuse to write output (exit 1),
    list the offending keys. This prevents shipping {{auc}} or a stale value.
  - --allow-missing renders unresolved keys as {{MISSING:key}} and exits 0 with a
    warning (for in-progress drafting only; never use for validated/PI output).

Every value written is therefore traceable to <data.json>:<key>, which is exactly
what number_provenance.py audit/gate verifies downstream.
"""
import argparse
import json
import re
import sys

PLACEHOLDER = re.compile(r"\{\{\s*([A-Za-z0-9_.\-]+)\s*\}\}")


def resolve(data, dotted):
    """Resolve a dotted path against nested dict/list. Return (found, value)."""
    cur = data
    for part in dotted.split("."):
        if isinstance(cur, dict) and part in cur:
            cur = cur[part]
        elif isinstance(cur, list):
            try:
                cur = cur[int(part)]
            except (ValueError, IndexError):
                return False, None
        else:
            return False, None
    return True, cur


def fmt(value):
    """Render a resolved value as text. Lists/dicts -> compact JSON."""
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False)
    return str(value)


def main():
    ap = argparse.ArgumentParser(description="By-construction report renderer (anti-fabrication Layer A).")
    ap.add_argument("template")
    ap.add_argument("data_json")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--allow-missing", action="store_true",
                    help="render unresolved keys as {{MISSING:key}} (drafting only; NOT for validated/PI)")
    args = ap.parse_args()

    try:
        with open(args.template, "r", encoding="utf-8") as f:
            tmpl = f.read()
        with open(args.data_json, "r", encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        print(f"[fill_report] ERROR reading inputs: {e}", file=sys.stderr)
        return 1

    missing = []

    def sub(m):
        key = m.group(1)
        found, val = resolve(data, key)
        if not found:
            missing.append(key)
            return f"{{{{MISSING:{key}}}}}"
        return fmt(val)

    rendered = PLACEHOLDER.sub(sub, tmpl)

    if missing:
        uniq = sorted(set(missing))
        print(f"[fill_report] {len(uniq)} placeholder(s) have NO source in {args.data_json}:",
              file=sys.stderr)
        for k in uniq:
            print(f"  - {{{{{k}}}}}", file=sys.stderr)
        if not args.allow_missing:
            print("[fill_report] REFUSING to write — every number must trace to the data file.\n"
                  "  Fix the data file or the placeholder. (--allow-missing for drafts only.)",
                  file=sys.stderr)
            return 1
        print("[fill_report] --allow-missing set: writing with {{MISSING:..}} markers (DRAFT only).",
              file=sys.stderr)

    try:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(rendered)
    except OSError as e:
        print(f"[fill_report] ERROR writing {args.out}: {e}", file=sys.stderr)
        return 1

    n_filled = len(PLACEHOLDER.findall(tmpl)) - len(missing)
    print(f"[fill_report] OK -> {args.out} ({n_filled} value(s) injected from {args.data_json})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
