#!/usr/bin/env bash
# Canonical raw-all producer-to-layered-v3 contract gate.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${SM_PYTHON:-/bip7_disk/liaoyoyo2001/miniconda3/bin/python3}"
METHOD="$ROOT/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts"

"$PY" "$METHOD/test_raw_all_production_contract.py"
"$PY" "$METHOD/test_build_longphase_raw_all_capture_receipt_v2.py"
"$PY" "$METHOD/test_finalize_longphase_raw_all_capture_receipts.py"
"$PY" "$METHOD/test_read_tag_sidecar_contract.py"
"$PY" "$METHOD/test_ssnv_site_ledger.py"
"$PY" "$METHOD/test_build_region_view_contract.py"
"$PY" "$ROOT/scripts/test_audit_canonical_longphase_bam_immutability.py"
"$PY" "$METHOD/test_backbone_sensitivity_contract.py"
"$PY" "$ROOT/scripts/test_layered_v3_lifecycle.py"
"$PY" "$ROOT/scripts/test_run_layered_v3.py"
"$PY" "$METHOD/test_layered_v3_contract.py"
"$PY" "$ROOT/scripts/test_verify_layered_v3.py"

printf '%s\n' "LAYERED V3 RAW-ALL CONTRACT TESTS: 103/103 PASS"
