#!/usr/bin/env python3
"""
Resilient Waterfall harness — Phase 5 EVALUATE gate.

Reads cycle artifacts, computes 6 risk components (precedent_similarity is
Path B and stays null), composite risk score, and applies two-stage verdict
(base + per-component override). Sweeps pitfalls_table.json against
plan.hypothesis. Writes evaluation.json conforming to evaluation.schema.json.

Usage:
    python3 _evaluator_run.py <cycle_id>

Exit codes:
    0 = approve_tier or downgrade_tier (cycle proceeds, possibly with caveat)
    1 = pending_review or exit_negative (human decision required)
    2 = error (artifacts missing, etc.)
"""

from __future__ import annotations

import json
import re
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]  # InterSubMod/
STATE_ROOT = REPO_ROOT / "state"
LEDGER = REPO_ROOT / "research" / "autoresearch" / "evidence_ledger.jsonl"
PITFALLS_PATH = Path(__file__).resolve().parent / "pitfalls_table.json"

# Composite weights (Path A: 6 components, precedent_similarity=null skipped)
WEIGHTS = {
    "multi_sample_consistency":  0.27,
    "effect_size_stability":     0.22,
    "precondition_freshness":    0.17,
    "subgroup_homogeneity":      0.17,
    "pitfall_coverage_score":    0.17,
}
# tier_support_alignment carried but not weighted into composite by default
# (it is informational; large mismatches surface via pitfall_hits)


def utcnow_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_json(path: Path) -> dict | None:
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as e:
        sys.stderr.write(f"WARN: {path} not valid JSON ({e}); skipping.\n")
        return None


def cycle_dir(cycle_id: str) -> Path:
    """Resolve cycle directory: prefer state/cycles/, fallback to state/retro_cycles/.

    Allows /run-evaluator to operate on retrospective fixtures (Day 6 Drill 1)
    without polluting active.json or live cycles/.
    """
    primary = STATE_ROOT / "cycles" / cycle_id
    if primary.is_dir():
        return primary
    retro = STATE_ROOT / "retro_cycles" / cycle_id
    if retro.is_dir():
        return retro
    return primary  # caller will hit the load-error path


def load_cycle_artifacts(cycle_id: str) -> dict[str, Any]:
    d = cycle_dir(cycle_id)
    artifacts = {
        "state":      load_json(d / "state.json"),
        "plan":       load_json(d / "plan.json"),
        "precheck":   load_json(d / "precheck.json"),
        "pilot":      load_json(d / "pilot.json"),
        "generalize": load_json(d / "generalize.json"),
    }
    if artifacts["state"] is None or artifacts["plan"] is None:
        sys.stderr.write(f"ERROR: cycle {cycle_id} missing required state.json or plan.json\n")
        sys.exit(2)
    return artifacts


# -------------------- Component calculations -------------------- #

def comp_multi_sample_consistency(generalize: dict | None) -> float | None:
    if not generalize:
        return None
    consistency = generalize.get("consistency", {})
    n_pass = consistency.get("n_samples_passed")
    n_total = consistency.get("n_samples_total")
    if n_total is None or n_total == 0:
        return None
    return max(0.0, min(1.0, n_pass / n_total))


def comp_effect_size_stability(generalize: dict | None) -> float | None:
    if not generalize:
        return None
    samples = generalize.get("samples", [])
    abs_values = [abs(s["metric_value"]) for s in samples if "metric_value" in s]
    if len(abs_values) < 2:
        return None
    mx = max(abs_values)
    if mx == 0:
        return 0.0
    return max(0.0, min(1.0, min(abs_values) / mx))


def comp_precondition_freshness(precheck: dict | None) -> float:
    if precheck is None:
        return 0.5  # neutral when missing
    verdict = precheck.get("verdict", "WARN")
    return {"PASS": 1.0, "WARN": 0.7, "BLOCKED": 0.3}.get(verdict, 0.5)


def comp_subgroup_homogeneity(pilot: dict | None, generalize: dict | None) -> float | None:
    """
    Two strategies, in order of preference:
    1. If pilot.metric_results contains subgroup breakdown (future schema), use stddev/mean.
    2. Fallback: cross-sample coefficient of variation from generalize.samples[].metric_value.
    """
    if generalize:
        samples = generalize.get("samples", [])
        values = [s["metric_value"] for s in samples if "metric_value" in s]
        if len(values) >= 3:
            mean = statistics.mean(values)
            if mean == 0:
                return 0.0
            std = statistics.stdev(values) if len(values) > 1 else 0.0
            cv = std / abs(mean)
            return max(0.0, min(1.0, 1.0 - cv))
    return None


def comp_pitfall_coverage_score(plan: dict, pilot: dict | None) -> tuple[float, list[dict]]:
    """
    Sweep pitfalls_table.json. For each entry whose trigger_keywords match
    plan.hypothesis, evaluate auto_check_rule. Returns (score, hits[]).

    score = 1 - n_unchecked_relevant / n_relevant
    if no relevant pitfalls → score = 1.0
    """
    hypothesis_text = (plan.get("hypothesis") or "").lower()
    pitfalls_data = load_json(PITFALLS_PATH)
    if not pitfalls_data:
        return 1.0, []

    pitfalls = pitfalls_data.get("pitfalls", [])
    hits = []
    n_relevant = 0
    n_unchecked = 0

    for pf in pitfalls:
        keywords = [k.lower() for k in pf.get("trigger_keywords", [])]
        if not any(k in hypothesis_text for k in keywords):
            continue

        n_relevant += 1
        # Heuristic check: is the relevant guard rail done?
        guard_passed = _evaluate_pitfall_guard(pf, pilot)
        evidence = ""
        if not guard_passed:
            n_unchecked += 1
            evidence = f"trigger_keyword matched in hypothesis; guard '{pf.get('detect_via', [])}' not pass"
            hits.append({
                "pitfall": pf["id"],
                "severity": pf.get("severity", "warn"),
                "evidence": evidence,
            })

    if n_relevant == 0:
        return 1.0, []
    return max(0.0, 1.0 - n_unchecked / n_relevant), hits


def _evaluate_pitfall_guard(pitfall: dict, pilot: dict | None) -> bool:
    """
    Conservative: require the guard artifact to be 'pass' in pilot.confound_guard.
    If pilot/confound_guard missing → assume guard NOT passed (safer).
    """
    if pilot is None:
        return False
    cg = pilot.get("confound_guard", {})

    # Map pitfall_id to confound_guard field (heuristic)
    guard_map = {
        "P-01": cg.get("af_bin_check"),
        "P-02": cg.get("within_group_ols"),
        "P-06": cg.get("n_reads_confound"),
    }
    guard_status = guard_map.get(pitfall["id"])
    if guard_status is None:
        # Pitfalls without a direct confound_guard field (P-03/P-04/P-05) → mark unchecked
        return False
    return guard_status == "pass"


def comp_tier_support_alignment(plan: dict, state: dict) -> float | None:
    """
    Compares plan-stated tier intent vs evidence_ledger stability for the
    same hypothesis_id. If ledger entry missing or no tier set, returns None.
    """
    hypothesis_id = plan.get("hypothesis_id") or state.get("hypothesis_id")
    tier_used = state.get("tier")  # current tier as integer-string or symbol
    try:
        tier_int = int(tier_used)
    except (TypeError, ValueError):
        return None  # tier=='pending' / 'NEGATIVE' → not applicable

    if hypothesis_id is None or not LEDGER.is_file():
        return None

    stability_grade = None
    for line in LEDGER.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            entry = json.loads(line)
        except json.JSONDecodeError:
            continue
        if entry.get("hypothesis_id") == hypothesis_id:
            sval = str(entry.get("stability") or "")
            m = re.search(r"\b([1-5])\b", sval)
            if m:
                stability_grade = int(m.group(1))

    if stability_grade is None:
        return None
    return min(stability_grade, tier_int) / tier_int


# -------------------- Verdict logic -------------------- #

def compute_verdict(components: dict) -> tuple[str, str, dict]:
    """Two-stage verdict: composite + per-component override."""
    contributing = [(k, v) for k, v in components.items() if v is not None and k in WEIGHTS]
    total_weight = sum(WEIGHTS[k] for k, _ in contributing)
    if total_weight == 0:
        return "approve_tier", "approve_tier", {"risk_base": 0.0, "n_low": 0, "n_critical": 0}

    risk_base = sum(WEIGHTS[k] * (1.0 - v) for k, v in contributing) / total_weight

    # Stage 1: base verdict
    if risk_base > 0.7:
        base_verdict = "pending_review"
    elif risk_base > 0.4:
        base_verdict = "downgrade_tier"
    else:
        base_verdict = "approve_tier"

    # Stage 2: per-component minima override
    in_use = [v for k, v in components.items() if v is not None and k in WEIGHTS]
    n_low = sum(1 for v in in_use if v < 0.4)
    n_critical = sum(1 for v in in_use if v < 0.2)

    if n_critical >= 1 or n_low >= 3:
        final_verdict = "pending_review"
    elif n_low >= 1:
        rank = ["approve_tier", "downgrade_tier", "pending_review", "exit_negative"]
        final_verdict = "downgrade_tier" if rank.index(base_verdict) < rank.index("downgrade_tier") else base_verdict
    else:
        final_verdict = base_verdict

    debug = {"risk_base": round(risk_base, 4), "n_low": n_low, "n_critical": n_critical}
    return base_verdict, final_verdict, debug


def tier_recommendation(state: dict, final_verdict: str) -> str:
    """Return string from {'1','2','3','4','5','NEGATIVE'} or 'pending'."""
    current = state.get("tier", "pending")
    if final_verdict == "exit_negative":
        return "NEGATIVE"
    if final_verdict == "pending_review":
        return current  # do not change; await human review
    if final_verdict == "downgrade_tier":
        try:
            return str(max(1, int(current) - 1))
        except (TypeError, ValueError):
            return current
    return current  # approve_tier: keep


# -------------------- Output rendering -------------------- #

def render_summary(cycle_id: str, payload: dict, debug: dict) -> str:
    lines = [f"[/run-evaluator] cycle_id={cycle_id}"]
    lines.append("  Components:")
    for k in WEIGHTS:
        v = payload["risk_components"].get(k)
        if v is None:
            lines.append(f"    [{k:<28}] —  (skipped, no data)")
            continue
        marker = "✓" if v >= 0.6 else ("⚠" if v >= 0.4 else "✗")
        lines.append(f"    [{k:<28}] {v:>4.2f} {marker}")
    tsa = payload["risk_components"].get("tier_support_alignment")
    if tsa is not None:
        lines.append(f"  tier_support_alignment (info): {tsa:.2f}")
    lines.append(f"  composite risk_base: {debug['risk_base']:.3f}")
    lines.append(f"  n_low (<0.4): {debug['n_low']}    n_critical (<0.2): {debug['n_critical']}")
    lines.append(f"  verdict: {payload['verdict']}")
    lines.append(f"  tier_recommendation: {payload['tier_recommendation']}")
    if payload["pitfall_hits"]:
        lines.append(f"  pitfall_hits: {len(payload['pitfall_hits'])}")
        for h in payload["pitfall_hits"]:
            lines.append(f"    - {h['pitfall']} [{h['severity']}] {h['evidence'][:80]}")
    return "\n".join(lines)


def write_evaluation(cycle_id: str, payload: dict) -> Path:
    out = cycle_dir(cycle_id) / "evaluation.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    return out


# -------------------- M2 failure attribution -------------------- #

# Component → primary category mapping for M2 attribution
_COMPONENT_CATEGORY = {
    "multi_sample_consistency": "consistency_violation",
    "effect_size_stability":    "consistency_violation",
    "precondition_freshness":   "precondition_violation",
    "subgroup_homogeneity":     "confound_bias",
    "pitfall_coverage_score":   "confound_bias",
}

# Pitfall id → category mapping
_PITFALL_CATEGORY = {
    "P-01": "confound_bias",
    "P-02": "confound_bias",
    "P-03": "dataset_schema",
    "P-04": "dataset_schema",
    "P-05": "confound_bias",
    "P-06": "confound_bias",
    "P-07": "precedent_pattern",
}


def build_failure_attribution(
    state: dict | None,
    components: dict,
    pitfall_hits: list[dict],
    final_verdict: str,
) -> dict | None:
    """M2 structured failure attribution; only set when verdict ≠ approve_tier."""
    if final_verdict == "approve_tier":
        return None

    categories: list[str] = []
    components_below: list[str] = []
    for k, v in components.items():
        if v is None or k not in WEIGHTS:
            continue
        if v < 0.4:
            components_below.append(f"{k}={v:.2f}")
            cat = _COMPONENT_CATEGORY.get(k)
            if cat:
                categories.append(cat)

    pitfall_ids = [h["pitfall"] for h in pitfall_hits]
    for pid in pitfall_ids:
        cat = _PITFALL_CATEGORY.get(pid)
        if cat:
            categories.append(cat)

    # Dedupe preserving order
    categories = list(dict.fromkeys(categories)) or ["unknown"]
    primary = categories[0]

    # Confidence: rule match (pitfall) > multi-component > single
    if pitfall_hits:
        confidence = "high"
    elif len(components_below) >= 2:
        confidence = "medium"
    else:
        confidence = "low"

    return {
        "categories": categories,
        "primary_category": primary,
        "phase_at_failure": (state or {}).get("phase", "P5_EVALUATE"),
        "skill_at_failure": "run-evaluator",
        "gate_at_failure": "P5_EVALUATOR",
        "components_below_threshold": components_below,
        "pitfalls_hit": pitfall_ids,
        "confidence": confidence,
        "diagnosed_at": utcnow_iso(),
        "auto_recovery_attempted": False,
        "auto_recovery_succeeded": None,
        "user_intervention_required": final_verdict in ("pending_review", "exit_negative"),
    }


def update_state_failure_attribution(cycle_id: str, attribution: dict | None) -> None:
    """Write/clear failure_attribution in state.json (denormalized snapshot)."""
    state_path = cycle_dir(cycle_id) / "state.json"
    if not state_path.exists():
        return
    try:
        state = json.loads(state_path.read_text())
    except json.JSONDecodeError:
        return
    if attribution is None:
        state.pop("failure_attribution", None)
    else:
        state["failure_attribution"] = attribution
    state_path.write_text(json.dumps(state, indent=2, ensure_ascii=False) + "\n")


# -------------------- Main -------------------- #

def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: _evaluator_run.py <cycle_id>\n")
        sys.exit(2)

    cycle_id = sys.argv[1]
    art = load_cycle_artifacts(cycle_id)

    components = {
        "multi_sample_consistency": comp_multi_sample_consistency(art["generalize"]),
        "effect_size_stability":    comp_effect_size_stability(art["generalize"]),
        "precondition_freshness":   comp_precondition_freshness(art["precheck"]),
        "subgroup_homogeneity":     comp_subgroup_homogeneity(art["pilot"], art["generalize"]),
    }
    coverage_score, pitfall_hits = comp_pitfall_coverage_score(art["plan"], art["pilot"])
    components["pitfall_coverage_score"] = coverage_score
    components["tier_support_alignment"] = comp_tier_support_alignment(art["plan"], art["state"])
    components["precedent_similarity"] = None  # Path A skip

    base_verdict, final_verdict, debug = compute_verdict(components)
    rec_tier = tier_recommendation(art["state"], final_verdict)

    required_followups = []
    if final_verdict == "pending_review":
        required_followups.append("Human review required: examine low/critical components and pitfall_hits before tier change.")
    if final_verdict == "downgrade_tier":
        required_followups.append("Tier reduced by 1 from current. Re-run /run-evaluator after addressing the lowest component.")
    if any(c is not None and c < 0.4 for k, c in components.items() if k in WEIGHTS):
        offenders = [k for k, c in components.items() if c is not None and c < 0.4 and k in WEIGHTS]
        required_followups.append(f"Components below 0.4 (need attention): {', '.join(offenders)}")

    # M2 failure attribution (only when verdict ≠ approve_tier)
    attribution = build_failure_attribution(art["state"], components, pitfall_hits, final_verdict)

    payload = {
        "schema_version": "1.0",
        "cycle_id": cycle_id,
        "ran_at": utcnow_iso(),
        "retraction_risk": debug["risk_base"],
        "risk_components": components,
        "pitfall_hits": pitfall_hits,
        "tier_recommendation": rec_tier,
        "required_followups": required_followups,
        "verdict": final_verdict,
    }
    if attribution is not None:
        payload["failure_attribution"] = attribution

    out_path = write_evaluation(cycle_id, payload)
    update_state_failure_attribution(cycle_id, attribution)  # mirror to state.json
    print(render_summary(cycle_id, payload, debug))
    print(f"\nWritten to: {out_path.relative_to(REPO_ROOT)}")
    print(f"Base verdict: {base_verdict} → Final verdict: {final_verdict} (override: {base_verdict != final_verdict})")
    if attribution:
        print(f"Failure attribution: {attribution['primary_category']} (confidence={attribution['confidence']}); state.json updated.")

    sys.exit(1 if final_verdict in ("pending_review", "exit_negative") else 0)


if __name__ == "__main__":
    main()
