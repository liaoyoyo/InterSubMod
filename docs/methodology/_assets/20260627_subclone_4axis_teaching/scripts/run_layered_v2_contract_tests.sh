#!/usr/bin/env bash
# Run the deterministic layered-v2 data-flow, read-tag, and lifecycle contracts.
set -euo pipefail

SCD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${SM_PYTHON:-/bip7_disk/liaoyoyo2001/miniconda3/bin/python3}"

"$PY" "$SCD/test_layered_reconstruction_v2.py"
"$PY" "$SCD/test_read_tag_sidecar_contract.py"
"$PY" "$SCD/test_ssnv_site_ledger.py"
"$PY" "$SCD/test_build_region_view_contract.py"
"$PY" "$SCD/test_raw_all_production_contract.py"
"$PY" "$SCD/test_build_longphase_raw_all_capture_receipt_v2.py"
"$PY" "$SCD/test_layered_runner_fail_closed.py"
"$PY" "$SCD/test_backbone_sensitivity_contract.py"
"$PY" "$SCD/test_longphase_production_closeout.py"
"$PY" "$SCD/test_parameter_sensitivity_contract.py"

printf '%s\n' "LAYERED V2 CONTRACT TESTS: 46/46 PASS"
