#!/bin/bash
# UserPromptSubmit hook — inject reminder about .md path output format.
#
# Rule: When listing .md file paths for the user to confirm/navigate, prefix with
# "InterSubMod/..." so the user can easily locate files from the project root.
# Example:
#   BAD:  docs/experiments/in_progress/2026/04/XX.md
#   GOOD: InterSubMod/docs/experiments/in_progress/2026/04/XX.md

cat <<'EOF'
{"hookSpecificOutput":{"hookEventName":"UserPromptSubmit","additionalContext":"[OUTPUT RULE — .md path format]\nWhen listing .md file paths for user confirmation or navigation, ALWAYS prefix with 'InterSubMod/...' (project root) so the user can locate files from the repository root.\n\nBAD:  docs/experiments/in_progress/2026/04/XX.md\nBAD:  /big7_disk/liaoyoyo2001/InterSubMod/docs/...\nGOOD: InterSubMod/docs/experiments/in_progress/2026/04/XX.md\n\nApplies to: new reports, existing reports, any file reference the user needs to open. Use relative-from-repo-root format starting with 'InterSubMod/'."}}
EOF
