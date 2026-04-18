#!/usr/bin/env python3
"""query_ledger.py - Query the evidence ledger for experiment history.

Search and filter the evidence_ledger.jsonl by verdict, hypothesis, sample,
date range, or free-text search. Outputs a formatted summary table.

Usage:
    # Show all entries
    python3 scripts/validation/ledger/query_ledger.py

    # Filter by verdict
    python3 scripts/validation/ledger/query_ledger.py --verdict positive

    # Filter by hypothesis keyword
    python3 scripts/validation/ledger/query_ledger.py --search "QS"

    # Filter by sample and date range
    python3 scripts/validation/ledger/query_ledger.py --sample HCC1395 --after 2026-04-01

    # Filter only validation_experiment entries
    python3 scripts/validation/ledger/query_ledger.py --type validation_experiment

    # JSON output for programmatic use
    python3 scripts/validation/ledger/query_ledger.py --verdict positive --json
"""

import argparse
import json
import os
import sys
from datetime import datetime


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
DEFAULT_LEDGER = os.path.join(PROJECT_ROOT, "research", "autoresearch", "evidence_ledger.jsonl")


def load_ledger(ledger_path):
    """Load all entries from JSONL ledger."""
    entries = []
    if not os.path.isfile(ledger_path):
        return entries

    with open(ledger_path, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                entries.append(json.loads(line))
            except json.JSONDecodeError as e:
                print(f"[WARN] Line {line_num}: invalid JSON: {e}", file=sys.stderr)

    return entries


def filter_entries(entries, verdict=None, search=None, sample=None,
                   hypothesis_id=None, entry_type=None, track=None,
                   after=None, before=None):
    """Filter ledger entries by criteria."""
    filtered = []

    for entry in entries:
        # Verdict filter
        if verdict and entry.get("verdict", "") != verdict:
            continue

        # Type filter
        if entry_type and entry.get("type", "") != entry_type:
            continue

        # Track filter
        if track:
            entry_track = entry.get("pipeline_track", "")
            entry_mode = entry.get("validation_mode", "")
            entry_cmode = entry.get("mode", "")
            if track == "to" and not any("to" in str(v).lower() for v in [entry_track, entry_mode, entry_cmode]):
                continue
            if track == "paired" and any("to" in str(v).lower() for v in [entry_track, entry_mode, entry_cmode]):
                continue

        # Sample filter
        if sample:
            datasets = entry.get("datasets_tested", [])
            result_f1 = entry.get("result_f1", {})
            if sample not in datasets and sample not in result_f1:
                continue

        # Hypothesis ID filter
        if hypothesis_id and entry.get("hypothesis_id", "") != hypothesis_id:
            continue

        # Free-text search (hypothesis, key_observations, cycle_id)
        if search:
            search_lower = search.lower()
            searchable = " ".join([
                str(entry.get("hypothesis", "")),
                str(entry.get("key_observations", "")),
                str(entry.get("cycle_id", "")),
                str(entry.get("hypothesis_id", "")),
                str(entry.get("human_notes", "")),
            ]).lower()
            if search_lower not in searchable:
                continue

        # Date filters
        ts = entry.get("timestamp", "")
        if after and ts:
            try:
                entry_date = ts[:10]  # YYYY-MM-DD portion
                if entry_date < after:
                    continue
            except (ValueError, IndexError):
                pass

        if before and ts:
            try:
                entry_date = ts[:10]
                if entry_date > before:
                    continue
            except (ValueError, IndexError):
                pass

        filtered.append(entry)

    return filtered


def format_table(entries, verbose=False):
    """Format entries as a readable table."""
    if not entries:
        return "No matching entries found."

    lines = []
    lines.append(f"Found {len(entries)} entries:\n")

    # Summary table
    lines.append(f"{'#':>3} | {'Date':10} | {'ID':15} | {'HypID':8} | {'Verdict':16} | {'Track':7} | {'Hypothesis':50}")
    lines.append("-" * 120)

    for i, entry in enumerate(entries, 1):
        ts = entry.get("timestamp", "")[:10] if entry.get("timestamp") else "N/A"
        cycle_id = entry.get("cycle_id", "N/A")[:15]
        hyp_id = entry.get("hypothesis_id", "")[:8]
        verdict = entry.get("verdict", "N/A")[:16]
        track = entry.get("pipeline_track", entry.get("validation_mode", ""))[:7]
        hypothesis = entry.get("hypothesis", "N/A")[:50]

        lines.append(f"{i:3d} | {ts:10} | {cycle_id:15} | {hyp_id:8} | {verdict:16} | {track:7} | {hypothesis}")

    if verbose:
        lines.append("")
        lines.append("=" * 80)
        lines.append("DETAILED RESULTS")
        lines.append("=" * 80)

        for i, entry in enumerate(entries, 1):
            lines.append(f"\n--- Entry #{i}: {entry.get('cycle_id', 'N/A')} ---")
            lines.append(f"  Hypothesis: {entry.get('hypothesis', 'N/A')}")
            lines.append(f"  Verdict: {entry.get('verdict', 'N/A')}")
            lines.append(f"  Human: {entry.get('human_decision', 'N/A')}")

            # F1 results
            result_f1 = entry.get("result_f1", {})
            delta_f1 = entry.get("delta_f1", {})
            if result_f1:
                lines.append(f"  F1 Results:")
                for sample in sorted(result_f1.keys()):
                    f1 = result_f1[sample]
                    delta = delta_f1.get(sample, "")
                    delta_str = f" (delta={delta:+.6f})" if isinstance(delta, (int, float)) else ""
                    lines.append(f"    {sample}: {f1:.4f}{delta_str}")

            # Cross-sample stats
            css = entry.get("cross_sample_stats", {})
            if css:
                lines.append(f"  Cross-sample: mean_delta={css.get('mean_delta_f1', 'N/A')}, "
                             f"p={css.get('paired_t_pvalue', 'N/A')}")

            # Key observations
            obs = entry.get("key_observations", "")
            if obs:
                lines.append(f"  Observations: {obs[:200]}...")

            lines.append(f"  Notes: {entry.get('human_notes', 'N/A')}")

    return "\n".join(lines)


def format_summary_stats(entries):
    """Generate summary statistics from entries."""
    if not entries:
        return ""

    verdicts = {}
    tracks = {}
    for entry in entries:
        v = entry.get("verdict", "unknown")
        verdicts[v] = verdicts.get(v, 0) + 1
        t = entry.get("pipeline_track", entry.get("validation_mode", "unknown"))
        tracks[t] = tracks.get(t, 0) + 1

    lines = [
        "\n--- Summary ---",
        f"Total entries: {len(entries)}",
        f"Verdicts: {', '.join(f'{k}={v}' for k, v in sorted(verdicts.items()))}",
        f"Tracks: {', '.join(f'{k}={v}' for k, v in sorted(tracks.items()))}",
    ]

    # Date range
    timestamps = [e.get("timestamp", "")[:10] for e in entries if e.get("timestamp")]
    if timestamps:
        lines.append(f"Date range: {min(timestamps)} to {max(timestamps)}")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Query the evidence ledger")
    parser.add_argument("--ledger", default=DEFAULT_LEDGER,
                        help="Path to evidence_ledger.jsonl")
    parser.add_argument("--verdict", default=None,
                        help="Filter by verdict (positive, negative, neutral, etc.)")
    parser.add_argument("--search", default=None,
                        help="Free-text search in hypothesis, observations, notes")
    parser.add_argument("--sample", default=None,
                        help="Filter by sample name (e.g., HCC1395)")
    parser.add_argument("--hypothesis-id", default=None,
                        help="Filter by hypothesis ID (e.g., H011)")
    parser.add_argument("--type", default=None,
                        help="Filter by entry type (e.g., validation_experiment)")
    parser.add_argument("--track", default=None, choices=["paired", "to"],
                        help="Filter by analysis track")
    parser.add_argument("--after", default=None,
                        help="Show entries after date (YYYY-MM-DD)")
    parser.add_argument("--before", default=None,
                        help="Show entries before date (YYYY-MM-DD)")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Show detailed results")
    parser.add_argument("--json", action="store_true",
                        help="Output as JSON array")
    parser.add_argument("--stats", action="store_true",
                        help="Show summary statistics")
    parser.add_argument("--last", type=int, default=None,
                        help="Show only the last N entries")
    args = parser.parse_args()

    entries = load_ledger(args.ledger)
    if not entries:
        print(f"[ERROR] No entries found in {args.ledger}", file=sys.stderr)
        sys.exit(1)

    filtered = filter_entries(
        entries,
        verdict=args.verdict,
        search=args.search,
        sample=args.sample,
        hypothesis_id=args.hypothesis_id,
        entry_type=args.type,
        track=args.track,
        after=args.after,
        before=args.before,
    )

    # Limit to last N
    if args.last and len(filtered) > args.last:
        filtered = filtered[-args.last:]

    if args.json:
        print(json.dumps(filtered, indent=2, ensure_ascii=False))
    else:
        print(format_table(filtered, verbose=args.verbose))
        if args.stats:
            print(format_summary_stats(filtered))


if __name__ == "__main__":
    main()
