#!/bin/bash
# cache_telemetry.sh — Analyze Claude Code transcript JSONL for cache hit rate.
#
# Purpose: Quantify Prompt Cache effectiveness (Anthropic 90% cost reduction claim).
# Source: P2 audit M7 — no telemetry was the gap; this hook closes it.
# Usage:
#   bash InterSubMod/scripts/hooks/cache_telemetry.sh              # latest session
#   bash InterSubMod/scripts/hooks/cache_telemetry.sh --all        # all sessions
#   bash InterSubMod/scripts/hooks/cache_telemetry.sh <path.jsonl> # specific transcript
# Output: stdout summary

set -o pipefail

export CT_TRANSCRIPT_DIR="/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod"
export CT_ARG="${1:-latest}"

python3 << 'PYEOF'
import json
import sys
import os
import datetime
from glob import glob

arg = os.environ.get("CT_ARG", "latest")
transcript_dir = os.environ.get("CT_TRANSCRIPT_DIR", "")

if arg == "latest":
    files = sorted(glob(f"{transcript_dir}/*.jsonl"), key=os.path.getmtime, reverse=True)
    targets = files[:1] if files else []
elif arg == "--all":
    targets = sorted(glob(f"{transcript_dir}/*.jsonl"), key=os.path.getmtime)
elif arg.endswith(".jsonl"):
    targets = [arg]
else:
    targets = sorted(glob(f"{transcript_dir}/*.jsonl"), key=os.path.getmtime, reverse=True)[:1]

if not targets:
    print("No transcripts found", file=sys.stderr)
    sys.exit(1)

# Pricing (Opus 4.7 per million tokens, approximation as of 2026-05)
PRICING = {"input": 15.0, "cache_read": 1.50, "cache_creation": 18.75, "output": 75.0}

agg = {"files": 0, "turns": 0, "input": 0, "cache_read": 0, "cache_creation": 0, "output": 0, "per_file": []}

for tpath in targets:
    fs = {"file": os.path.basename(tpath), "turns": 0, "input": 0, "cache_read": 0, "cache_creation": 0, "output": 0}
    try:
        with open(tpath, encoding="utf-8") as f:
            for line in f:
                try:
                    d = json.loads(line)
                    usage = d.get("message", {}).get("usage", {})
                    if usage:
                        fs["turns"] += 1
                        fs["input"] += usage.get("input_tokens", 0)
                        fs["cache_read"] += usage.get("cache_read_input_tokens", 0)
                        fs["cache_creation"] += usage.get("cache_creation_input_tokens", 0)
                        fs["output"] += usage.get("output_tokens", 0)
                except Exception:
                    pass
    except Exception as e:
        print(f"Error reading {tpath}: {e}", file=sys.stderr)
        continue

    agg["files"] += 1
    for k in ("turns", "input", "cache_read", "cache_creation", "output"):
        agg[k] += fs[k]
    agg["per_file"].append(fs)

total_in = agg["input"] + agg["cache_read"] + agg["cache_creation"]
cache_hit_rate = agg["cache_read"] / max(1, total_in)
cache_create_rate = agg["cache_creation"] / max(1, total_in)

cost_actual = (agg["input"] * PRICING["input"] + agg["cache_read"] * PRICING["cache_read"]
               + agg["cache_creation"] * PRICING["cache_creation"] + agg["output"] * PRICING["output"]) / 1e6
cost_no_cache = (total_in * PRICING["input"] + agg["output"] * PRICING["output"]) / 1e6
savings = cost_no_cache - cost_actual
savings_pct = savings * 100 / max(0.0001, cost_no_cache)

print("# cache_telemetry — Prompt Cache effectiveness")
print(f"Date:           {datetime.datetime.now().isoformat(timespec='seconds')}")
print(f"Files scanned:  {agg['files']}")
print(f"Total turns:    {agg['turns']}")
print()
print("## Token totals")
print(f"  input (uncached):      {agg['input']:>15,}")
print(f"  cache_read:            {agg['cache_read']:>15,}")
print(f"  cache_creation:        {agg['cache_creation']:>15,}")
print(f"  output:                {agg['output']:>15,}")
print()
print("## Cache effectiveness")
print(f"  Cache hit rate:        {cache_hit_rate * 100:>6.1f}%  (cache_read / total_input)")
print(f"  Cache create rate:     {cache_create_rate * 100:>6.1f}%")
print(f"  Uncached rate:         {agg['input'] / max(1, total_in) * 100:>6.1f}%")
print()
print("## Cost analysis (Opus 4.7 approx pricing)")
print(f"  Actual cost:           ${cost_actual:>10.4f}")
print(f"  Without cache (est):   ${cost_no_cache:>10.4f}")
print(f"  Savings:               ${savings:>10.4f}  ({savings_pct:.1f}%)")
print()
if cache_hit_rate < 0.5:
    print("⚠ Cache hit rate < 50% — review static prefix (CLAUDE.md / AGENTS.md / rules)")
elif cache_hit_rate < 0.7:
    print("🟡 Cache hit rate 50-70% — moderate. Improvement possible.")
else:
    print(f"✅ Cache hit rate ≥ 70% — good static/dynamic separation. (Anthropic claim: up to 90%)")

if len(agg["per_file"]) > 1:
    print()
    print("## Top 5 sessions by cost")
    def session_cost(x):
        return (x["input"]*PRICING["input"] + x["cache_read"]*PRICING["cache_read"]
                + x["cache_creation"]*PRICING["cache_creation"] + x["output"]*PRICING["output"]) / 1e6
    for s in sorted(agg["per_file"], key=lambda x: -session_cost(x))[:5]:
        print(f"  {s['file'][:40]:<42} ${session_cost(s):>8.4f} (turns={s['turns']})")
PYEOF
