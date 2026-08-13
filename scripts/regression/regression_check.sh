#!/bin/bash
# ISM frozen regression check: schema-v2 provenance + numeric output must match reviewed goldens.
# Golden artifacts are intentionally immutable here; updates require a separate reviewed procedure.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
BIN="$ROOT/build/bin/inter_sub_mod"
TUMOR=/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam
NORMAL=/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam
REF=/big8_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa
VCF=/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp_chr1.vcf.gz
GOLDEN_SKIP="$HERE/golden_chr1_w5000_bernoulli_skip.tsv"
GOLDEN_MAXDIST="$HERE/golden_chr1_w5000_bernoulli.tsv"

COLUMNS=(
  RegionID NumReads NumCpGs GlobalP CramersV PassedGating ClusterPermanovaValid Significant
  VerificationSchemaVersion VerificationClass VerificationClass_V1_Deprecated VerificationClass_Legacy
  LabelFirstSupport ClusterFirstSupport WithinHPSupport DispersionWarning EvidencePath EvidenceDerivation
)
EXPECTED_HEADER="$(IFS=$'\t'; echo "${COLUMNS[*]}")"

if [ "${1:-}" = "--update-golden" ]; then
  echo "[regression] REFUSED: golden artifacts are frozen and cannot be rebuilt by regression_check.sh." >&2
  echo "[regression] Use a separate, code-reviewed golden migration with recorded provenance." >&2
  exit 2
fi
if [ "$#" -ne 0 ]; then
  echo "Usage: $0" >&2
  exit 2
fi

validate_contract_file() {  # $1=file $2=label
  local file="$1"
  local label="$2"
  if [ ! -f "$file" ]; then
    echo "[regression] $label missing: $file" >&2
    exit 2
  fi
  local observed_header
  observed_header="$(head -n 1 "$file")"
  if [ "$observed_header" != "$EXPECTED_HEADER" ]; then
    echo "[regression] $label has a stale/non-v2 header; frozen golden migration is required." >&2
    echo "[regression] expected: $EXPECTED_HEADER" >&2
    echo "[regression] observed: $observed_header" >&2
    exit 2
  fi

  python3 - "$file" "$ROOT" <<'PY'
import sys
from pathlib import Path

import pandas as pd

path = Path(sys.argv[1])
sys.path.insert(0, sys.argv[2])
from scripts.lib.verification_schema_contract import (  # noqa: E402
    CURRENT_TO_EVIDENCE_PATH,
    SchemaContractError,
    read_evidence,
    select_current_view,
    select_legacy_view,
)

try:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    current = select_current_view(frame)
    legacy = select_legacy_view(frame)
    evidence = read_evidence(frame)
    if current.unknown_counts:
        raise SchemaContractError(f"unknown canonical classes: {current.unknown_counts}")

    expected_v1 = current.values.map(
        lambda value: "Strong" if value in {"Strong_Bidirectional", "ClusterFirstOnly"} else value
    )
    mismatch = frame["VerificationClass_V1_Deprecated"].astype(str) != expected_v1.astype(str)
    if mismatch.any():
        raise SchemaContractError("VerificationClass_V1_Deprecated is not the exact frozen alias")

    expected_path = current.values.map(CURRENT_TO_EVIDENCE_PATH)
    if (evidence["EvidencePath"].astype(str) != expected_path.astype(str)).any():
        raise SchemaContractError("EvidencePath conflicts with canonical VerificationClass")
    expected_label = legacy.values.isin({"Strong", "Weak"})
    expected_cluster = legacy.values.isin({"Strong", "Subclone"})
    if (evidence["LabelFirstSupport"].astype(bool) != expected_label).any():
        raise SchemaContractError("LabelFirstSupport conflicts with VerificationClass_Legacy")
    if (evidence["ClusterFirstSupport"].astype(bool) != expected_cluster).any():
        raise SchemaContractError("ClusterFirstSupport conflicts with VerificationClass_Legacy")
    legacy_derived = evidence["EvidenceDerivation"] == "LEGACY_CLASS"
    if evidence.loc[legacy_derived, "WithinHPSupport"].notna().any():
        raise SchemaContractError("LEGACY_CLASS requires WithinHPSupport=NA")
    if evidence.loc[legacy_derived, "DispersionWarning"].notna().any():
        raise SchemaContractError("LEGACY_CLASS requires DispersionWarning=NA")
except (SchemaContractError, ValueError) as exc:
    raise SystemExit(f"[regression][schema-contract] {path}: {exc}") from exc
PY
}

# Validate both frozen goldens before launching the expensive C++ process.
validate_contract_file "$GOLDEN_SKIP" "SKIP golden"
validate_contract_file "$GOLDEN_MAXDIST" "MAX_DIST golden"

# The frozen schema guards above must fail before any host-specific output setup.
# This keeps the contract check side-effect free on clean runners.
export TMPDIR=/big7_disk/liaoyoyo2001/tmp
mkdir -p "$TMPDIR"

run_extract() {  # $1=nan_strategy $2=outdir
  local strategy="$1"
  local outdir="$2"
  if [ -e "$outdir" ]; then
    echo "[regression] refusing stale output directory: $outdir" >&2
    exit 2
  fi
  mkdir -p "$outdir"
  "$BIN" -t "$TUMOR" -n "$NORMAL" -r "$REF" -v "$VCF" -w 5000 -j 16 \
    --distance-metric BERNOULLI --nan-distance-strategy "$strategy" \
    --no-output-distance-matrix -o "$outdir" > "$outdir/run.log" 2>&1
  python3 - "$outdir/significance_summary.csv" "$outdir/current.tsv" <<'PY'
import sys
import pandas as pd

columns = [
    "RegionID", "NumReads", "NumCpGs", "GlobalP", "CramersV", "PassedGating",
    "ClusterPermanovaValid", "Significant", "VerificationSchemaVersion", "VerificationClass",
    "VerificationClass_V1_Deprecated", "VerificationClass_Legacy", "LabelFirstSupport",
    "ClusterFirstSupport", "WithinHPSupport", "DispersionWarning", "EvidencePath",
    "EvidenceDerivation",
]
frame = pd.read_csv(sys.argv[1], keep_default_na=False)
serialized = pd.read_csv(sys.argv[1], dtype=str, keep_default_na=False)
missing = [column for column in columns if column not in frame.columns]
if missing:
    raise SystemExit("[regression][schema-contract] generated summary missing: " + ", ".join(missing))
# Retain the C++ writer's canonical lowercase true/false/NA evidence values.
# Other historical regression columns keep their established pandas formatting.
for column in (
    "LabelFirstSupport", "ClusterFirstSupport", "WithinHPSupport",
    "DispersionWarning", "EvidencePath", "EvidenceDerivation",
):
    frame[column] = serialized[column]
frame[columns].sort_values("RegionID").to_csv(
    sys.argv[2], sep="\t", index=False, float_format="%.6g"
)
PY
  validate_contract_file "$outdir/current.tsv" "$strategy current output"
}

OUT_SKIP="$ROOT/output/_tmp_regression_skip.$$"
OUT_MAXD="$ROOT/output/_tmp_regression_maxdist.$$"
FAIL=0

check_one() {  # $1=label $2=strategy $3=golden $4=outdir
  local label="$1"
  local strategy="$2"
  local golden="$3"
  local outdir="$4"
  run_extract "$strategy" "$outdir"
  if diff -q "$golden" "$outdir/current.tsv" >/dev/null; then
    echo "[regression] ✅ $label PASS — schema/provenance and values match frozen golden ($(($(wc -l < "$golden")-1)) rows)"
  else
    echo "[regression] 🔴 $label FAIL — frozen result changed; first differences:"
    diff -u "$golden" "$outdir/current.tsv" | head -15 || true
    echo "[regression] full current output retained at $outdir"
    FAIL=1
  fi
}

check_one "SKIP(default/main)" SKIP "$GOLDEN_SKIP" "$OUT_SKIP"
check_one "MAX_DIST(control)" MAX_DIST "$GOLDEN_MAXDIST" "$OUT_MAXD"

if [ "$FAIL" -eq 0 ]; then
  echo "[regression] ✅ dual frozen guard PASS (outputs retained: $OUT_SKIP, $OUT_MAXD)"
  exit 0
fi
echo "[regression] 🔴 dual frozen guard FAIL"
exit 1
