#!/usr/bin/env bash
# preflight.sh — verify html-preview skill dependencies.
set -euo pipefail
errors=0

# 1. Python 3.9+
PY_VER=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
if python3 -c "import sys; sys.exit(0 if sys.version_info >= (3, 9) else 1)"; then
    echo "[OK] python3 ${PY_VER}"
else
    echo "[FAIL] python3 ${PY_VER} (need >= 3.9)" >&2
    errors=$((errors+1))
fi

# 2. markdown library
if python3 -c "import markdown" 2>/dev/null; then
    echo "[OK] markdown $(python3 -c 'import markdown; print(markdown.__version__)')"
else
    echo "[FAIL] markdown missing. Install: pip install --user 'markdown>=3.5'" >&2
    errors=$((errors+1))
fi

# 3. jinja2 library
if python3 -c "import jinja2" 2>/dev/null; then
    echo "[OK] jinja2 $(python3 -c 'import jinja2; print(jinja2.__version__)')"
else
    echo "[FAIL] jinja2 missing. Install: pip install --user 'jinja2>=3.0'" >&2
    errors=$((errors+1))
fi

# 4. beautifulsoup4 library
if python3 -c "import bs4" 2>/dev/null; then
    echo "[OK] bs4 $(python3 -c 'import bs4; print(bs4.__version__)')"
else
    echo "[FAIL] bs4 missing. Install: pip install --user 'beautifulsoup4>=4.12'" >&2
    errors=$((errors+1))
fi

# 5. pyyaml (Phase 1 dep)
if python3 -c "import yaml" 2>/dev/null; then
    echo "[OK] pyyaml"
else
    echo "[FAIL] pyyaml missing. Install: pip install --user pyyaml" >&2
    errors=$((errors+1))
fi

# 6. Pillow (Phase 1 dep)
if python3 -c "from PIL import Image" 2>/dev/null; then
    echo "[OK] Pillow"
else
    echo "[FAIL] Pillow missing. Install: pip install --user Pillow" >&2
    errors=$((errors+1))
fi

if [ "$errors" -gt 0 ]; then
    echo "" >&2
    echo "$errors dependency check(s) failed." >&2
    exit 1
fi
