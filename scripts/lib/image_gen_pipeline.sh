#!/usr/bin/env bash
# image_gen_pipeline.sh — Batch image generation pipeline using OpenAI gpt-image-2
#
# Reads prompts from <prompts_dir>/*.txt, generates PNG via OpenAI Images API,
# writes to <output_dir>/<prompt_basename>.png. Idempotent (skips if PNG exists).
#
# Usage:
#   image_gen_pipeline.sh <prompts_dir> <output_dir> [--dry-run] [--quality high|medium|low] [--size 1024x1024|1024x1536]
#
# Examples:
#   # Real run (requires OPENAI_API_KEY + openai python pkg)
#   ./scripts/lib/image_gen_pipeline.sh \
#     docs/presentations/validated/2026/05/self_phasing_synthesis_PI/prompts \
#     docs/presentations/validated/2026/05/self_phasing_synthesis_PI/figures
#
#   # Dry run (show commands without executing)
#   ./scripts/lib/image_gen_pipeline.sh prompts/ figures/ --dry-run
#
# Prerequisites (real run):
#   1. pip install --user openai      # or: pip3 install openai
#   2. export OPENAI_API_KEY="sk-..."
#   3. (optional) Python 3.8+ with openai>=1.0
#
# Cost estimate (gpt-image-2 high quality, 1024x1024): ~$0.15-0.20 per image
#
# After generation: open each PNG in viewer + use Claude Read tool for vision verify
#
set -euo pipefail

# ─── Parse args ────────────────────────────────────────────────────────────
DRY_RUN=0
QUALITY="high"
SIZE="1024x1024"
PROMPTS_DIR=""
OUTPUT_DIR=""

while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run) DRY_RUN=1; shift ;;
        --quality) QUALITY="$2"; shift 2 ;;
        --size)    SIZE="$2"; shift 2 ;;
        --help|-h)
            head -25 "$0" | grep -E "^#" | sed 's/^# \?//'
            exit 0 ;;
        *)
            if [ -z "$PROMPTS_DIR" ]; then PROMPTS_DIR="$1"
            elif [ -z "$OUTPUT_DIR" ]; then OUTPUT_DIR="$1"
            else echo "ERROR: unexpected arg: $1" >&2; exit 1
            fi
            shift ;;
    esac
done

if [ -z "$PROMPTS_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 <prompts_dir> <output_dir> [--dry-run] [--quality high|medium|low] [--size 1024x1024|1024x1536]" >&2
    exit 1
fi

if [ ! -d "$PROMPTS_DIR" ]; then
    echo "ERROR: prompts_dir does not exist: $PROMPTS_DIR" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# ─── Pre-flight checks (real run only) ─────────────────────────────────────
if [ "$DRY_RUN" -eq 0 ]; then
    if [ -z "${OPENAI_API_KEY:-}" ]; then
        echo "ERROR: OPENAI_API_KEY not set in environment." >&2
        echo "  Run: export OPENAI_API_KEY=\"sk-...\"" >&2
        echo "  Or use --dry-run to preview commands without executing." >&2
        exit 1
    fi
    if ! python3 -c "import openai" 2>/dev/null; then
        echo "ERROR: openai Python package not installed." >&2
        echo "  Run: pip3 install --user 'openai>=1.0'" >&2
        echo "  Or use --dry-run to preview commands without executing." >&2
        exit 1
    fi
    OPENAI_VER=$(python3 -c "import openai; print(openai.__version__)")
    echo "[ok] openai pkg $OPENAI_VER, API key set (${OPENAI_API_KEY:0:7}...)"
fi

# ─── Log file ──────────────────────────────────────────────────────────────
LOG="$OUTPUT_DIR/.gen_log.jsonl"
echo "[log] $LOG"

# ─── Process each prompt ───────────────────────────────────────────────────
COUNT=0
SKIPPED=0
FAILED=0

for prompt_file in "$PROMPTS_DIR"/*.txt; do
    [ -e "$prompt_file" ] || { echo "[warn] no .txt files found in $PROMPTS_DIR"; break; }

    name=$(basename "$prompt_file" .txt)
    out="$OUTPUT_DIR/${name}.png"

    if [ -f "$out" ]; then
        echo "[skip] $name already exists ($out)"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    if [ "$DRY_RUN" -eq 1 ]; then
        echo "[dry-run] would generate: $name"
        echo "  prompt: $prompt_file ($(wc -c < "$prompt_file") chars)"
        echo "  output: $out"
        echo "  size:   $SIZE"
        echo "  qual:   $QUALITY"
        echo "  cmd: python3 -c \"...openai.images.generate(prompt=..., model='gpt-image-2')\""
        echo ""
        COUNT=$((COUNT + 1))
        continue
    fi

    echo "[gen] $name → $out ($SIZE, $QUALITY)"
    PROMPT=$(cat "$prompt_file")

    # Inline Python to avoid creating extra .py file
    if python3 - "$out" "$SIZE" "$QUALITY" <<EOF
import os, sys, base64, json
from datetime import datetime
from openai import OpenAI

out_path, size, quality = sys.argv[1], sys.argv[2], sys.argv[3]
prompt = """$PROMPT"""

client = OpenAI()
try:
    resp = client.images.generate(
        model="gpt-image-2",
        prompt=prompt,
        size=size,
        quality=quality,
    )
    b64 = resp.data[0].b64_json
    with open(out_path, "wb") as f:
        f.write(base64.b64decode(b64))
    revised = getattr(resp.data[0], "revised_prompt", None) or ""
    log = {
        "ts": datetime.utcnow().isoformat(),
        "name": "$name",
        "out": out_path,
        "size": size,
        "quality": quality,
        "model": "gpt-image-2",
        "prompt_chars": len(prompt),
        "revised_prompt_chars": len(revised),
        "status": "ok",
    }
    print(json.dumps(log))
except Exception as e:
    log = {
        "ts": datetime.utcnow().isoformat(),
        "name": "$name",
        "out": out_path,
        "status": "error",
        "error": str(e)[:200],
    }
    print(json.dumps(log))
    sys.exit(1)
EOF
    then
        echo "[ok] $name done"
        COUNT=$((COUNT + 1))
    else
        echo "[err] $name failed (see log)"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "═════════════════════════════════════════════════"
echo "  Generated: $COUNT"
echo "  Skipped:   $SKIPPED (already exist)"
echo "  Failed:    $FAILED"
echo "═════════════════════════════════════════════════"

if [ "$DRY_RUN" -eq 0 ] && [ "$COUNT" -gt 0 ]; then
    echo ""
    echo "Next step: vision-verify each PNG:"
    echo "  ls $OUTPUT_DIR/*.png | xargs -I{} echo 'Read tool: {}'"
    echo "  → Compare against $PROMPTS_DIR/*.txt"
fi
