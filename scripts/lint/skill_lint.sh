#!/usr/bin/env bash
# skill_lint.sh — Validate SKILL.md frontmatter format for all .claude/skills/
#
# Per Anthropic best-practice article 4 (Skills): description must include
#   1. Action verb at start (directive, not descriptive noun phrase)
#   2. USE WHEN clause (positive trigger keywords)
#   3. SKIP WHEN clause (negative trigger / boundary statement)
#   4. evals.json existence (warned, not failed)
#
# Exit 0 if all required checks PASS; exit 1 if any FAIL.
# WARN counts (missing evals.json) do not affect exit code.
#
# Usage: bash scripts/lint/skill_lint.sh [--skills-dir PATH]

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SKILLS_DIR="${REPO_ROOT}/.claude/skills"

if [[ "${1:-}" == "--skills-dir" && -n "${2:-}" ]]; then
    SKILLS_DIR="$2"
fi

if [[ ! -d "${SKILLS_DIR}" ]]; then
    echo "ERROR: skills dir not found: ${SKILLS_DIR}" >&2
    exit 2
fi

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
# Extract description value from a SKILL.md frontmatter block.
# Handles both single-line `description: ...` and multi-line `description: |` form.
extract_description() {
    local file="$1"
    awk '
        BEGIN { in_fm=0; in_desc=0; desc=""; line_count=0 }
        /^---[[:space:]]*$/ {
            line_count++
            if (line_count == 1) { in_fm=1; next }
            if (line_count == 2) { exit }
        }
        in_fm == 1 && /^description:[[:space:]]*\|[[:space:]]*$/ {
            in_desc=2; next
        }
        in_fm == 1 && /^description:[[:space:]]*/ {
            sub(/^description:[[:space:]]*/, "")
            desc = $0
            in_desc = 1
            next
        }
        in_desc == 2 && /^[[:space:]]+/ {
            sub(/^[[:space:]]+/, "")
            desc = desc " " $0
            next
        }
        in_desc == 2 && /^[^[:space:]]/ {
            in_desc = 0
        }
        END { print desc }
    ' "$file"
}

# Extract name field
extract_name() {
    local file="$1"
    awk '
        BEGIN { in_fm=0; line_count=0 }
        /^---[[:space:]]*$/ {
            line_count++
            if (line_count == 1) { in_fm=1; next }
            if (line_count == 2) { exit }
        }
        in_fm == 1 && /^name:[[:space:]]*/ {
            sub(/^name:[[:space:]]*/, "")
            print
            exit
        }
    ' "$file"
}

# Check if description starts with action verb (directive form).
# Heuristic: flag if starts with these descriptive patterns:
#   - "Skill" / "skill"
#   - "X 的" (CJK noun phrase)
#   - "A/An <noun>" without verb
# We accept if first word is a verb (not "Skill"/"This"/"A"/"An" + noun pattern)
# or if the description starts with a non-Latin (CJK) action word.
is_directive_form() {
    local desc="$1"
    local first_word
    first_word=$(printf '%s' "$desc" | awk '{print $1}')
    # Strip leading markdown emphasis
    first_word="${first_word#\*\*}"
    first_word="${first_word#__}"

    # Bad patterns (descriptive-form indicators)
    case "$first_word" in
        Skill|skill|"This"|"A"|"An")
            return 1
            ;;
    esac

    # Bad: ends with "的" right after first word (CJK possessive noun phrase)
    if printf '%s' "$desc" | grep -qE '^[^[:space:]]+的[[:space:]]'; then
        return 1
    fi

    # Otherwise treat as directive
    return 0
}

# Check description contains USE WHEN trigger keyword
has_use_when() {
    local desc="$1"
    if printf '%s' "$desc" | grep -qE '(USE WHEN|觸發|TRIGGER when|Use when)'; then
        return 0
    fi
    return 1
}

# Check description contains SKIP WHEN clause
has_skip_when() {
    local desc="$1"
    if printf '%s' "$desc" | grep -qE '(SKIP WHEN|DO NOT USE WHEN|不適用|SKIP:)'; then
        return 0
    fi
    return 1
}

# -----------------------------------------------------------------------------
# Main scan
# -----------------------------------------------------------------------------
TOTAL=0
PASS_COUNT=0
FAIL_COUNT=0
WARN_COUNT=0

# Header
printf '%-32s %-18s %-10s %-12s %-6s %s\n' \
    "Skill" "Description form" "USE WHEN" "SKIP WHEN" "Evals" "Status"
printf '%-32s %-18s %-10s %-12s %-6s %s\n' \
    "--------------------------------" "------------------" "----------" "------------" "------" "------"

# Iterate over each skill directory
while IFS= read -r -d '' skill_md; do
    skill_dir="$(dirname "$skill_md")"
    skill_name="$(basename "$skill_dir")"
    TOTAL=$((TOTAL + 1))

    # Extract fields
    name=$(extract_name "$skill_md")
    desc=$(extract_description "$skill_md")

    # Run checks
    fail_reasons=()
    warn_reasons=()

    # 1. name field present
    if [[ -z "$name" ]]; then
        fail_reasons+=("missing-name")
    fi

    # 2. description present
    if [[ -z "$desc" ]]; then
        fail_reasons+=("missing-description")
        form_label="missing"
        use_label="-"
        skip_label="-"
    else
        # 3. Form check
        if is_directive_form "$desc"; then
            form_label="directive"
        else
            form_label="descriptive"
            fail_reasons+=("descriptive-form")
        fi

        # 4. USE WHEN check
        if has_use_when "$desc"; then
            use_label="OK"
        else
            use_label="MISS"
            fail_reasons+=("no-use-when")
        fi

        # 5. SKIP WHEN check
        if has_skip_when "$desc"; then
            skip_label="OK"
        else
            skip_label="MISS"
            fail_reasons+=("no-skip-when")
        fi
    fi

    # 6. evals.json (warn-only)
    if [[ -f "${skill_dir}/evals.json" ]]; then
        evals_label="OK"
    else
        evals_label="MISS"
        warn_reasons+=("no-evals")
        WARN_COUNT=$((WARN_COUNT + 1))
    fi

    # Status verdict
    if [[ ${#fail_reasons[@]} -eq 0 ]]; then
        status="PASS"
        PASS_COUNT=$((PASS_COUNT + 1))
    else
        status="FAIL [$(IFS=,; echo "${fail_reasons[*]}")]"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi

    printf '%-32s %-18s %-10s %-12s %-6s %s\n' \
        "$skill_name" "$form_label" "$use_label" "$skip_label" "$evals_label" "$status"

done < <(find "${SKILLS_DIR}" -mindepth 2 -maxdepth 2 -name SKILL.md -print0 | sort -z)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo
echo "Summary"
echo "-------"
printf '  Total skills scanned : %d\n' "$TOTAL"
printf '  PASS                 : %d\n' "$PASS_COUNT"
printf '  FAIL                 : %d\n' "$FAIL_COUNT"
printf '  WARN (no evals.json) : %d\n' "$WARN_COUNT"
echo

if [[ $FAIL_COUNT -gt 0 ]]; then
    echo "Result: FAIL — $FAIL_COUNT skill(s) need description fixes"
    exit 1
fi

echo "Result: PASS — all $TOTAL skills have valid description format"
exit 0
