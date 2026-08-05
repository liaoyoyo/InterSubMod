#!/usr/bin/env bash
set -euo pipefail

# The packaged renderer verifies the report before the compatibility patch is
# inserted. Hide non-overlay scrollbars only for that initial packaging pass;
# the final HTML is re-verified with the unwrapped Chromium executable.
: "${REAL_CHROMIUM:?Set REAL_CHROMIUM to the chrome-headless-shell executable}"
exec "$REAL_CHROMIUM" --hide-scrollbars "$@"
