#!/bin/bash
# pre_tier_upgrade_check.sh
# PreToolUse hook (matcher: Edit|Write) — block tier ⭐4 / ⭐5 upgrade
# unless the corresponding cycle has an evaluation.json (P5 EVALUATE
# completed) OR an explicit override comment is present.
#
# Triggers: Edit/Write of
#   - InterSubMod/state/cycles/*/state.json
#   - InterSubMod/docs/experiments/INDEX.md
#   - InterSubMod/docs/CURRENT_FOCUS.md
#
# Detection: regex on new_string for tier upgrade markers (⭐4 / ⭐5 /
# "tier": "4" / "tier": "5"), then check the relevant cycle's
# evaluation.json. If missing AND no override comment → exit 2 (block).
#
# Override syntax (in the same Edit/Write payload):
#   <!-- tier-upgrade-override: <reason> -->
#
# Override audit trail: appended to InterSubMod/state/invalidation/tier_overrides.jsonl

set -u

REPO_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
AUDIT_LOG="${REPO_ROOT}/state/invalidation/tier_overrides.jsonl"

# gate_exit: hard-block (exit 2) ONLY for the structured authoritative path
# (state/cycles/*/state.json). For prose docs (INDEX.md / CURRENT_FOCUS.md) the
# ⭐4/⭐5 marker fires constantly in narrative ("⭐4 需 COLO829"), so blocking there
# would be a false-positive nightmare — degrade to ADVISORY (exit 0 + warning).
# Cross-artifact tier consistency on prose docs is covered by /provenance-tier-audit.
gate_exit() {
    if [ "${WATCH_KIND:-}" = "state_json" ]; then
        exit 2
    fi
    echo "  ↑ ADVISORY only (prose doc); the ENFORCED exit-2 gate is on state/cycles/*/state.json." >&2
    echo "    For INDEX.md / CURRENT_FOCUS.md tier consistency, run /provenance-tier-audit." >&2
    exit 0
}

INPUT=$(cat)

FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)
NEW_STRING=$(echo "$INPUT" | jq -r '.tool_input.new_string // .tool_input.content // empty' 2>/dev/null)

[ -z "$FILE_PATH" ] && exit 0
[ -z "$NEW_STRING" ] && exit 0

# Bail unless the file is one of the 3 watched paths
case "$FILE_PATH" in
    *"/state/cycles/"*"/state.json")
        WATCH_KIND="state_json"
        ;;
    *"/docs/experiments/INDEX.md")
        WATCH_KIND="index_md"
        ;;
    *"/docs/CURRENT_FOCUS.md")
        WATCH_KIND="current_focus"
        ;;
    *)
        exit 0
        ;;
esac

# Detect a tier-upgrade pattern: ⭐4 / ⭐5 OR "tier": "4" / "tier": "5"
# (also support tier=4 / tier 4 / Tier ⭐4 with tolerant matching)
if ! echo "$NEW_STRING" | grep -qE '(⭐[45]|"tier"[[:space:]]*:[[:space:]]*"[45]"|[Tt]ier.*[45])'; then
    exit 0
fi

# Resolve cycle_id
CYCLE_ID=""
case "$WATCH_KIND" in
    state_json)
        # Path: .../state/cycles/<cycle_id>/state.json
        CYCLE_ID=$(echo "$FILE_PATH" | sed -nE 's|.*/state/cycles/([^/]+)/state\.json|\1|p')
        ;;
    *)
        # For INDEX.md / CURRENT_FOCUS.md, cycle_id may be embedded in new_string.
        # Look for canonical cycle_id pattern YYYYMMDD-HHMM-{slug}
        CYCLE_ID=$(echo "$NEW_STRING" | grep -oE '[0-9]{8}-[0-9]{4}-[a-z0-9_-]+' | head -1)
        ;;
esac

# Check for explicit override comment in new_string
if echo "$NEW_STRING" | grep -qE '<!--[[:space:]]*tier-upgrade-override:[[:space:]]*[^[:space:]]'; then
    REASON=$(echo "$NEW_STRING" | sed -nE 's|.*<!--[[:space:]]*tier-upgrade-override:[[:space:]]*([^>]*)-->.*|\1|p' | head -1 | sed 's/[[:space:]]*$//')
    TS=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    mkdir -p "$(dirname "$AUDIT_LOG")"
    jq -n -c \
        --arg ts "$TS" \
        --arg file "$FILE_PATH" \
        --arg cycle "${CYCLE_ID:-unknown}" \
        --arg kind "$WATCH_KIND" \
        --arg reason "$REASON" \
        '{timestamp: $ts, file: $file, cycle_id: $cycle, kind: $kind, reason: $reason, action: "override_allowed"}' \
        >> "$AUDIT_LOG"
    echo "[TierUpgradeHook] tier ⭐4-5 upgrade override accepted for cycle '${CYCLE_ID:-unknown}': ${REASON}" >&2
    echo "  audit logged → state/invalidation/tier_overrides.jsonl" >&2
    exit 0
fi

# No override; check if evaluation.json exists for this cycle
if [ -n "$CYCLE_ID" ]; then
    EVAL_PATH="${REPO_ROOT}/state/cycles/${CYCLE_ID}/evaluation.json"
    if [ -f "$EVAL_PATH" ]; then
        # Verify tier_recommendation matches the upgrade target
        REC=$(jq -r '.tier_recommendation // empty' "$EVAL_PATH" 2>/dev/null)
        VERDICT=$(jq -r '.verdict // empty' "$EVAL_PATH" 2>/dev/null)
        if [ "$VERDICT" = "approve_tier" ]; then
            echo "[TierUpgradeHook] tier upgrade allowed: cycle ${CYCLE_ID} has evaluation.json (verdict=approve_tier, tier_recommendation=${REC})." >&2
            exit 0
        else
            echo "[TierUpgradeHook] BLOCKED: cycle ${CYCLE_ID} evaluation.json verdict='${VERDICT}' (not approve_tier)." >&2
            RISK=$(jq -r '.retraction_risk // "n/a"' "$EVAL_PATH" 2>/dev/null)
            echo "  retraction_risk: ${RISK}" >&2
            echo "  Required: re-run /run-evaluator after addressing low components, OR add override:" >&2
            echo "  <!-- tier-upgrade-override: <reason ≤300 chars> -->" >&2
            gate_exit
        fi
    fi

    # evaluation.json missing
    echo "[TierUpgradeHook] BLOCKED: cycle ${CYCLE_ID} has no evaluation.json (P5 EVALUATE not run)." >&2
    echo "  Run: /run-evaluator ${CYCLE_ID}" >&2
    echo "  OR add explicit override:" >&2
    echo "  <!-- tier-upgrade-override: <reason ≤300 chars> -->" >&2
    gate_exit
fi

# Could not resolve cycle_id (e.g. INDEX.md upgrade row without recognizable id)
echo "[TierUpgradeHook] BLOCKED: tier ⭐4-5 marker detected in ${FILE_PATH} but no cycle_id recognized." >&2
echo "  Either include the cycle_id (YYYYMMDD-HHMM-slug) near the tier change," >&2
echo "  OR add explicit override:" >&2
echo "  <!-- tier-upgrade-override: <reason; legacy entries pre-dating harness need this> -->" >&2
gate_exit
