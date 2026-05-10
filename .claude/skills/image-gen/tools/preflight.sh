#!/usr/bin/env bash
# preflight.sh — verify image-gen skill dependencies before run.
# Exits 0 if all OK, non-zero with error message if any dep missing.
set -euo pipefail

errors=0

# 1. codex CLI logged in via ChatGPT OAuth
if ! codex login status 2>&1 | grep -q "Logged in"; then
    echo "[FAIL] codex CLI not logged in. Run: codex login" >&2
    errors=$((errors+1))
else
    echo "[OK] codex login"
fi

# 2. pycairo
if ! python3 -c "import cairo" 2>/dev/null; then
    echo "[FAIL] pycairo missing. Install: pip install --user pycairo" >&2
    errors=$((errors+1))
else
    echo "[OK] pycairo $(python3 -c 'import cairo; print(cairo.version)')"
fi

# 3. Pillow
if ! python3 -c "from PIL import Image" 2>/dev/null; then
    echo "[FAIL] Pillow missing. Install: pip install --user Pillow" >&2
    errors=$((errors+1))
else
    echo "[OK] Pillow $(python3 -c 'import PIL; print(PIL.__version__)')"
fi

# 4. pyyaml
if ! python3 -c "import yaml" 2>/dev/null; then
    echo "[FAIL] pyyaml missing. Install: pip install --user pyyaml" >&2
    errors=$((errors+1))
else
    echo "[OK] pyyaml $(python3 -c 'import yaml; print(yaml.__version__)')"
fi

if [ "$errors" -gt 0 ]; then
    echo "" >&2
    echo "$errors dependency check(s) failed. See messages above." >&2
    exit 1
fi
