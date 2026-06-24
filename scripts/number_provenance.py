#!/usr/bin/env python3
"""number_provenance.py — Layer B (write-time gate) + Layer C (task-end audit table).

Anti-fabrication mechanism for the 2026-06-01 "fabricated metric in HTML" postmortem.
The bug: an AI types a metric (AUC=0.985) from *expectation*, not from a file it just
read. Pure text rules failed twice (postmortem §9) -> only a mechanical check works.

This script extracts "metric-shaped" numbers from a report and checks whether each one
literally appears in a bounded set of source files (declared data_sources / sibling
_assets/ / same-dir .json|.tsv|.txt|.csv). A number with no source is a fabrication
candidate.

Subcommands
-----------
gate   : (hook) read Claude Code PreToolUse JSON on stdin. TIERED:
           - validated/ or pi_reports/ path  -> STRICT: unsourced metric => exit 2 (block)
           - any other path                  -> ADVISORY: print reminder, exit 0
         FAIL-OPEN: any error / no content / no sources => exit 0 (a broken hook must
         never block all writes; this is NOT the `|| exit 0` neutering bug — genuine
         detections in strict paths still exit 2).
         Override: a `<!-- provenance-verified: reason -->` marker in the content
         (or `data_sources:` frontmatter listing every value's source) => exit 0.

audit  : (CLI) python3 number_provenance.py audit <report> [--sources P ...]
         Print a markdown provenance table: metric value | found in (file:line) | status.
         For task-end recording of "the basis for each number" (Layer C).

Design for the 真實 × 消耗 × 錯誤率 balance:
  真實  : catches hand-typed numbers absent from any source file.
  消耗  : deterministic substring scan, flat-rate ~0 token cost; bounded file set.
  錯誤率: only "metric-shaped" tokens (>=2-decimal floats, AUC=/p=/%/Δ/n=...) are
          checked; years/dates/section refs skipped; substring match tolerates trailing
          precision (report 0.866 found inside source 0.8662). Derived numbers computed
          in-head are the known false-negative hole -> covered by Layer A (fill_report.py).
"""
import argparse
import json
import os
import re
import sys

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod"

# ---- metric-shaped number extraction ----------------------------------------
# P1: floats with >=2 decimals  -> 0.985, 44.89, 0.8662  (years/section ints excluded)
_P1 = re.compile(r"(?<![\w.])\d+\.\d{2,}(?![\w])")
# P2: tagged metrics            -> AUC=0.5, p<1e-3, OR=0.194, n=40, Δβ=-0.122, median=0.97
_P2 = re.compile(
    r"(?:AUC|OR|RR|HR|F1|R2|R²|p|P|r|n|N|median|mean|Δβ|Δ|delta)\s*[=<>:]\s*"
    r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
)
# P3: percentages               -> 44.89%, 5%
_P3 = re.compile(r"(\d+(?:\.\d+)?)\s*%")

# tokens to ignore even if they match (low information / structural)
_DATE = re.compile(r"\b\d{4}-\d{2}-\d{2}\b")
_YEAR = re.compile(r"^(?:19|20)\d{2}$")

# range sanity: metrics physically bounded to [0,1] — out-of-range = definitely wrong/fabricated.
# Zero-false-positive subset only: % skipped (increases can exceed 100%), r skipped (ambiguous radius/count).
# Word-boundary before name prevents "step=3"->p, "fp=5"/"hp=1"->p false matches. 2026-06-25 架構稽核 §13 第4道.
_RANGE = re.compile(
    r"(?<![A-Za-z0-9])(AUC|Jaccard|Cramér'?s?\s*V|CramérV|Cramer\s*V|p-value|p)\s*[=<>:]\s*"
    r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
)


def range_sanity(text):
    """Return [(snippet, metric_label, value)] for metrics whose value is impossible:
    AUC / p / Jaccard / Cramér's V must be in [0,1]. Never throws (used inside fail-open gate)."""
    bad = []
    try:
        for m in _RANGE.finditer(_DATE.sub(" ", text)):
            name, raw = m.group(1), m.group(2)
            try:
                val = float(raw)
            except ValueError:
                continue
            if 0.0 <= val <= 1.0:
                continue
            n = name.replace("'", "").replace(" ", "").lower()
            label = ("AUC" if n.startswith("auc") else "Jaccard" if n.startswith("jaccard")
                     else "Cramér's V" if "v" in n else "p-value")
            bad.append((m.group(0).strip(), label, val))
    except Exception:
        return []
    return bad

SOURCE_EXT = (".json", ".tsv", ".txt", ".csv", ".t.txt", ".log", ".jsonl")
MAX_SOURCE_FILES = 80
MAX_SCAN_BYTES = 20 * 1024 * 1024  # per-file cap


def extract_metrics(text):
    """Return a sorted unique list of metric-shaped numeric token strings."""
    text_nodate = _DATE.sub(" ", text)  # drop dates so 2026-06-01 doesn't leak ints
    toks = set()
    toks.update(_P1.findall(text_nodate))
    toks.update(_P2.findall(text_nodate))
    toks.update(_P3.findall(text_nodate))
    out = set()
    for t in toks:
        t = t.strip()
        if not t:
            continue
        # drop bare years and tiny integers (1, 2, 3 ...) — not "metric-shaped"
        if _YEAR.match(t):
            continue
        if "." not in t and "e" not in t.lower():
            # pure integer: only keep if reasonably large (sample sizes etc.)
            try:
                if abs(int(t)) < 10:
                    continue
            except ValueError:
                continue
        out.add(t)
    return sorted(out)


# ---- source gathering + scanning --------------------------------------------
def _read_frontmatter(text):
    """Extract YAML-ish frontmatter block (--- ... ---) or HTML comment header."""
    fm = {}
    m = re.match(r"^---\n(.*?)\n---\n", text, re.DOTALL)
    block = m.group(1) if m else ""
    for line in block.splitlines():
        if ":" in line:
            k, _, v = line.partition(":")
            fm[k.strip()] = v.strip()
    return fm, text


def gather_sources(report_path, declared):
    """Bounded source file set (respects search-discipline §12 — no repo-root recursion)."""
    files = []
    seen = set()

    def add(p):
        ap = os.path.abspath(p)
        if ap in seen or not os.path.isfile(ap):
            return
        try:
            if os.path.getsize(ap) > 200 * 1024 * 1024:  # skip multi-GB BAM-like
                return
        except OSError:
            return
        seen.add(ap)
        files.append(ap)

    for p in declared:
        add(p if os.path.isabs(p) else os.path.join(ROOT, p))

    if report_path:
        d = os.path.dirname(os.path.abspath(report_path))
        for cand_dir in (d, os.path.join(d, "_assets"),
                         os.path.join(d, os.path.splitext(os.path.basename(report_path))[0] + "_assets")):
            if os.path.isdir(cand_dir):
                for name in sorted(os.listdir(cand_dir)):
                    if name.lower().endswith(SOURCE_EXT):
                        add(os.path.join(cand_dir, name))
                    if len(files) >= MAX_SOURCE_FILES:
                        break
    return files[:MAX_SOURCE_FILES]


def find_token(token, sources):
    """Return 'file:line' where token appears as a substring, else None."""
    for path in sources:
        try:
            scanned = 0
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                for i, line in enumerate(f, 1):
                    scanned += len(line)
                    if scanned > MAX_SCAN_BYTES:
                        break
                    if token in line:
                        rel = os.path.relpath(path, ROOT)
                        return f"{rel}:{i}"
        except OSError:
            continue
    return None


# ---- gate (Layer B) ----------------------------------------------------------
def is_strict_path(p):
    return bool(p) and ("/docs/reports/validated/" in p or "/docs/reports/pi_reports/" in p
                        or "/docs/experiments/validated/" in p or "/pi_reports/" in p)


def cmd_gate():
    """PreToolUse hook entry. FAIL-OPEN on any error (exit 0)."""
    try:
        raw = sys.stdin.read()
        data = json.loads(raw) if raw.strip() else {}
    except (ValueError, OSError):
        return 0
    tool = data.get("tool_name", "")
    if tool not in ("Edit", "Write"):
        return 0
    ti = data.get("tool_input", {})
    fpath = ti.get("file_path", "") or ""
    if not (fpath.endswith(".md") or fpath.endswith(".html")):
        return 0
    content = ti.get("content") or ti.get("new_string") or ""
    if not content:
        return 0

    # explicit human/AI override -> trust
    if "provenance-verified:" in content:
        return 0

    fm, _ = _read_frontmatter(content)
    declared = []
    if fm.get("data_sources"):
        declared = [s.strip().strip("[]\"' ") for s in fm["data_sources"].split(",") if s.strip()]

    metrics = extract_metrics(content)
    rng = range_sanity(content)  # range-impossible values (AUC/p/Jaccard/V ∉ [0,1]) = definitely wrong
    if not metrics and not rng:
        return 0

    sources = gather_sources(fpath, declared)
    unsourced = [m for m in metrics if find_token(m, sources) is None] if metrics else []

    strict = is_strict_path(fpath)
    if not unsourced and not rng:
        return 0

    base = os.path.basename(fpath)
    rng_block = ""
    if rng:
        rng_block = ("\n  🔴 RANGE-IMPOSSIBLE（定義上不可能，必為錯/捏造）: "
                     + "; ".join(f"{s}（{lbl} 須∈[0,1]）" for s, lbl, _v in rng[:6]))
    if strict:
        unsourced_block = ""
        if unsourced:
            unsourced_block = (f"\n  {len(unsourced)} metric(s) found NO source file: "
                               f"{', '.join(unsourced[:12])}{' ...' if len(unsourced) > 12 else ''}")
        msg = (f"[number_provenance] 🔴 BLOCK Write to validated/PI path: {base}"
               f"{unsourced_block}{rng_block}\n"
               f"  → 數字須溯源且範圍合理。Fix: 寫檔→讀回→貼真值；或 frontmatter `data_sources:`；或 fill_report.py。\n"
               f"  → Override (genuinely sourced elsewhere): <!-- provenance-verified: <來源> -->")
        print(msg, file=sys.stderr)
        return 2
    # advisory tier
    unsourced_block = ""
    if unsourced:
        unsourced_block = (f" {len(unsourced)} metric(s) not found in nearby sources: "
                           f"{', '.join(unsourced[:8])}{' ...' if len(unsourced) > 8 else ''}")
    msg = (f"[number_provenance] ⚠ {base}:{unsourced_block}{rng_block}\n"
           f"  → 確認每數字 Read 自檔(非預期) + 範圍合理。Postmortem 2026-06-01.")
    print(msg)
    return 0


# ---- audit (Layer C) ---------------------------------------------------------
def cmd_audit(report, extra_sources):
    try:
        with open(report, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
    except OSError as e:
        print(f"ERROR reading {report}: {e}", file=sys.stderr)
        return 1
    fm, _ = _read_frontmatter(content)
    declared = list(extra_sources)
    if fm.get("data_sources"):
        declared += [s.strip().strip("[]\"' ") for s in fm["data_sources"].split(",") if s.strip()]

    metrics = extract_metrics(content)
    sources = gather_sources(report, declared)

    print(f"# Number Provenance Audit — {os.path.relpath(os.path.abspath(report), ROOT)}\n")
    print(f"> sources scanned ({len(sources)}): "
          + (", ".join(os.path.relpath(s, ROOT) for s in sources) if sources else "**none found**"))
    print(f"> metrics detected: {len(metrics)}\n")
    print("| metric | source (file:line) | status |")
    print("|--------|--------------------|--------|")
    n_ok = 0
    for m in metrics:
        loc = find_token(m, sources)
        if loc:
            n_ok += 1
            print(f"| `{m}` | {loc} | ✅ sourced |")
        else:
            print(f"| `{m}` | — | 🔴 NO SOURCE |")
    n_bad = len(metrics) - n_ok
    print(f"\n**{n_ok}/{len(metrics)} sourced; {n_bad} unsourced.** "
          + ("✅ all numbers traceable." if n_bad == 0
             else "🔴 unsourced numbers must be verified or removed before validated/PI."))
    rng = range_sanity(content)
    if rng:
        print("\n**🔴 RANGE-IMPOSSIBLE（必為錯，AUC/p/Jaccard/V 須 ∈ [0,1]）:**")
        for s, lbl, v in rng:
            print(f"- `{s}` — {lbl} 實為 {v}，超出 [0,1]")
    else:
        print("\n✅ range sanity: no impossible AUC/p/Jaccard/V values.")
    return 0  # audit is a reporting tool — always exits 0, never blocks


def main():
    ap = argparse.ArgumentParser(description="Number provenance: gate (hook) + audit (table).")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("gate")
    pa = sub.add_parser("audit")
    pa.add_argument("report")
    pa.add_argument("--sources", nargs="*", default=[])
    args = ap.parse_args()
    if args.cmd == "gate":
        return cmd_gate()
    if args.cmd == "audit":
        return cmd_audit(args.report, args.sources)
    return 0


if __name__ == "__main__":
    sys.exit(main())
