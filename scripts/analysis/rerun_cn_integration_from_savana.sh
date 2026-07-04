#!/bin/bash
# Per-sample downstream wiring: SAVANA CN -> SM_CNBED -> refill cn in sm_region_integration.
#
# Run AFTER `savana STATUS=DONE` for the sample. It:
#   1) converts the SAVANA segmented-CN TSV to a gain/loss/loh SM_CNBED (savana_to_smcnbed.py)
#   2) backs up the existing (cn=unknown) integration, then re-runs sm_region_integration.py
#      with SM_CNBED set + SM_CN_SPARSE=0 (SAVANA is genome-wide, non-sparse)
#   3) prints the new by_cn
#
# STOPS at sm_region_integration (fills cn=unknown). The rest of the chain
# (topology_analysis -> candidate_scoring -> build_topology_workstation) is a
# separate downstream step (parallel-session territory) -- printed as NEXT.
#
# Usage: rerun_cn_integration_from_savana.sh <SAMPLE>
set -uo pipefail
S=${1:?usage: $0 <SAMPLE>}
REPO=/big7_disk/liaoyoyo2001/InterSubMod
MSROOT=/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone
SAVDIR=/big7_disk/liaoyoyo2001/cnv_sv_work/$S/savana_wgs
TSV=$SAVDIR/cna_normalhet/${S}_segmented_absolute_copy_number.tsv
BED=$SAVDIR/${S}_smcnbed.bed
SMINT=$REPO/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_region_integration.py
WORK=$MSROOT/$S

# --- preflight (§13.7: verify inputs before claiming) ---
[ -f "$TSV" ] || { echo "SKIP $S: SAVANA CN tsv not found ($TSV) — finish SAVANA first"; exit 1; }
st=$(cat "$SAVDIR/STATUS" 2>/dev/null || echo NONE)
[ "$st" = "DONE" ] || echo "WARN $S: savana STATUS=$st (expected DONE) — proceeding, verify results"
[ -f "$WORK/sm_linkage_genomewide.json" ] || { echo "SKIP $S: no sm_linkage in $WORK"; exit 1; }

# --- 1) SAVANA CN -> SM_CNBED ---
python3 "$REPO/scripts/analysis/savana_to_smcnbed.py" "$TSV" -o "$BED" || exit 2

# --- 2) backup cn=unknown version, refill cn ---
[ -f "$WORK/sm_region_integration.json" ] && \
  cp -f "$WORK/sm_region_integration.json" "$WORK/sm_region_integration.cn_unknown.bak"
SM_WORKDIR="$WORK" SM_CNBED="$BED" SM_CN_SPARSE=0 python3 "$SMINT" || exit 3

# --- 3) report ---
python3 -c "import json; d=json.load(open('$WORK/sm_region_integration.json')); \
print('$S new by_cn:', d['aggregate']['by_cn'])"
echo "OK $S: cn filled. NEXT (downstream, separate step): topology_analysis -> candidate_scoring -> build_topology_workstation"
