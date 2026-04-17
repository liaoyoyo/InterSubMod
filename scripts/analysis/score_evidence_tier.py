#!/usr/bin/env python3
"""Score evidence_ledger.jsonl entries by tier_flags.

Usage:
    python3 scripts/analysis/score_evidence_tier.py                # audit all entries
    python3 scripts/analysis/score_evidence_tier.py --last         # audit last entry only
    python3 scripts/analysis/score_evidence_tier.py --validate     # exit 1 if tier > confidence_cap
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
LEDGER = PROJECT_ROOT / "research/autoresearch/evidence_ledger.jsonl"

FLAG_TIER_MIN = {
    "within_group_ols": 3,
    "af_bin_stratified": 3,
    "permutation_tested": 3,
    "nested_cv": 3,
    "multi_method_agreement": 4,
    "multi_sample_consistent": 4,
    "literature_aligned": 5,
    "independent_reproduction": 5,
}


def compute_confidence_cap(flags: list[str]) -> int:
    """Return the max tier achievable given the flags present."""
    if not flags:
        return 1
    # Tier 3 requires any one tier-3 flag
    # Tier 4 requires any one tier-4 flag AND at least one tier-3 flag
    # Tier 5 requires any tier-5 flag AND tier-3 + tier-4 requirements met
    has_t3 = any(FLAG_TIER_MIN.get(f, 99) <= 3 for f in flags)
    has_t4 = any(FLAG_TIER_MIN.get(f, 99) == 4 for f in flags)
    has_t5 = any(FLAG_TIER_MIN.get(f, 99) == 5 for f in flags)

    if has_t5 and has_t4 and has_t3:
        return 5
    if has_t4 and has_t3:
        return 4
    if has_t3:
        return 3
    return 2  # multi-sample without controls


def audit(entries: list[dict], strict: bool = False) -> int:
    """Audit entries; return count of violations."""
    violations = 0
    for i, e in enumerate(entries, 1):
        hid = e.get("hypothesis_id", "unknown")
        tier = e.get("tier")
        flags = e.get("tier_flags", [])
        cap = compute_confidence_cap(flags)

        status = "OK"
        if tier is None:
            status = "MISSING tier"
        elif tier > cap:
            status = f"OVER-CLAIMED (tier={tier} > cap={cap})"
            violations += 1
        elif tier < cap:
            status = f"UNDER-CLAIMED (tier={tier} < cap={cap})"

        print(f"  [{i:3d}] {hid:30s} tier={tier} cap={cap} flags={len(flags):2d} → {status}")

    print(f"\nTotal violations: {violations}")
    if strict and violations > 0:
        return 1
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--last", action="store_true", help="Audit only the last entry")
    parser.add_argument("--validate", action="store_true",
                        help="Exit 1 if any over-claimed tier found")
    parser.add_argument("--ledger", type=Path, default=LEDGER)
    args = parser.parse_args()

    if not args.ledger.exists():
        print(f"[error] ledger not found: {args.ledger}", file=sys.stderr)
        return 2

    entries = []
    with args.ledger.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                entries.append(json.loads(line))
            except json.JSONDecodeError as e:
                print(f"[warn] skip malformed line: {e}", file=sys.stderr)

    if args.last:
        entries = entries[-1:]

    return audit(entries, strict=args.validate)


if __name__ == "__main__":
    sys.exit(main())
