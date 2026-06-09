#!/bin/bash
set -x
BT="$(cd "$(dirname "$0")" && pwd)"
export TMPDIR=/big7_disk/liaoyoyo2001/tmp
cd "$BT"
echo "=== Approach 1: gh release binary (nanoporetech/modkit latest) ==="
gh release download --repo nanoporetech/modkit --pattern '*x86_64.tar.gz' --dir "$BT" --clobber 2>&1 || echo "gh-download-FAILED"
ls -la "$BT"
for f in "$BT"/*.tar.gz; do [ -e "$f" ] && tar xzf "$f" -C "$BT" 2>&1; done
MK=$(find "$BT" -name modkit -type f 2>/dev/null | head -1)
if [ -n "$MK" ] && "$MK" --version 2>/dev/null; then echo "MODKIT_OK=$MK"; exit 0; fi
echo "=== Approach 2: conda bioconda ==="
conda create -y -p "$BT/modkit_env" -c bioconda -c conda-forge modkit 2>&1 | tail -25
MK2="$BT/modkit_env/bin/modkit"
if [ -x "$MK2" ] && "$MK2" --version 2>/dev/null; then echo "MODKIT_OK=$MK2"; exit 0; fi
echo "MODKIT_INSTALL_FAILED_BOTH"
