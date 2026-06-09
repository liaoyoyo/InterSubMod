#!/bin/bash
set -e
B="$(cd "$(dirname "$0")/.." && pwd)"
export TMPDIR=/big7_disk/liaoyoyo2001/tmp
MK=$(find "$B/tools" -name modkit -type f | head -1)
REF=/big8_disk/liaoyoyo2001/InterSubMod/data/ref/GRCh38_no_alt_analysis_set.fasta
LOCI=$B/loci.bed
for mode in somatic family; do
  echo "########## MODE=$mode ##########"
  BAMM=$B/runs/HCC1395_${mode}_hp.bam
  OUT=$B/runs/pileup_$mode; rm -rf "$OUT"; mkdir -p "$OUT"
  "$MK" pileup --phased --cpg --combine-strands --modified-bases C:m C:h -r "$REF" --threads 8 "$BAMM" "$OUT/" 2>&1 | tail -3
  echo "--- $mode pileup files ---"; for f in "$OUT"/*.bedmethyl; do echo "$(basename $f): $(wc -l <$f) rows"; done
  for hp in hp1 hp2; do
    if [ -s "$OUT/$hp.bedmethyl" ]; then
      sort -k1,1 -k2,2n "$OUT/$hp.bedmethyl" | bgzip -f > "$OUT/$hp.bedmethyl.gz"; tabix -f -p bed "$OUT/$hp.bedmethyl.gz"
    fi
  done
  if [ -s "$OUT/hp1.bedmethyl.gz" ] && [ -s "$OUT/hp2.bedmethyl.gz" ]; then
    echo "--- $mode dmr pair (a=HP1, b=HP2) ---"
    "$MK" dmr pair -a "$OUT/hp1.bedmethyl.gz" -b "$OUT/hp2.bedmethyl.gz" --ref "$REF" -r "$LOCI" -o "$B/runs/dmr_$mode.tsv" --header 2>&1 | tail -4 || echo "dmr_$mode FAILED"
  else
    echo "$mode: hp1/hp2 empty, skip dmr"
  fi
done
echo "MODKIT_CHAIN_DONE"
