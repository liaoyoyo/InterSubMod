#!/bin/bash
# V6 per-chr ALT HP1:HP2 ratio via parallel awk + intermediate files
set -e

REPO=/big7_disk/liaoyoyo2001/InterSubMod
WORK=/tmp/v6_per_chr_work
mkdir -p $WORK
rm -f $WORK/*

process_source() {
    local NAME=$1
    local BASE=$2
    [ -d "$BASE" ] || { echo "  [SKIP] $NAME"; return; }
    local OUTFILE=$WORK/${NAME}.txt
    echo "  $NAME → $OUTFILE"
    # Use awk to count per (chr, hp_class) tuples and emit lines, then aggregate
    find "$BASE" -name reads.tsv -print0 | \
        xargs -0 -P 16 awk -F$'\t' \
            'NR>1 && $8=="ALT" && $9=="1" {
                if ($7=="1" || $7=="1-1" || $7=="11") print $3 "\tHP1"
                else if ($7=="2" || $7=="2-1" || $7=="21") print $3 "\tHP2"
            }' 2>/dev/null | \
        sort | uniq -c > $OUTFILE
    wc -l $OUTFILE
}

# Process all 7 sources
process_source V3F_HCC1395 $REPO/research/paired_priority_bug_audit/phaseC_genome_three_way/V3F_off_tp
process_source V5_HCC1395 $REPO/research/paired_priority_bug_audit/phaseC_genome_three_way/V5_off_tp
process_source V6_HCC1395 $REPO/research/paired_priority_bug_audit/phaseC_genome_three_way/V6_off_tp
process_source V6_H1437 $REPO/research/paired_priority_bug_audit/phaseD_v6_5sample/H1437/off_tp
process_source V6_H2009 $REPO/research/paired_priority_bug_audit/phaseD_v6_5sample/H2009/off_tp
process_source V6_HCC1954 $REPO/research/paired_priority_bug_audit/phaseD_v6_5sample/HCC1954/off_tp
process_source V6_HCC1937 $REPO/research/paired_priority_bug_audit/phaseD_v6_5sample/HCC1937/off_tp

echo
echo "=== Aggregating final results ==="
python3 <<'PYEOF'
from collections import defaultdict
import os, glob
data = defaultdict(lambda: defaultdict(lambda: [0,0]))
for f in glob.glob("/tmp/v6_per_chr_work/*.txt"):
    name = os.path.basename(f).replace(".txt", "")
    with open(f) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) != 3: continue
            n, chrom, hp = int(parts[0]), parts[1], parts[2]
            if hp == "HP1": data[name][chrom][0] += n
            elif hp == "HP2": data[name][chrom][1] += n

chrs = ["chr"+str(i) for i in range(1,23)] + ["chrX","chrY"]

print("\n=== Per-chr ALT HP1:HP2 ratio (somatic ALT reads, tumor only) ===")
print(f"{'chr':<8}", end="")
for name in sorted(data):
    print(f" {name:>22}", end="")
print()
for c in chrs:
    print(f"{c:<8}", end="")
    for name in sorted(data):
        h1, h2 = data[name].get(c, (0,0))
        if h1 + h2 == 0:
            print(f" {'n/a':>22}", end="")
        else:
            r = h1/max(h2,1)
            print(f" {h1:>8}:{h2:<8}={r:>4.2f}", end="")
    print()

print("\n=== Summary ===")
print(f"{'Source':<22} {'genome HP1':>11} {'HP2':>10} {'ratio':>7} {'偏HP1>5x':>10} {'平衡0.5-2':>11} {'偏HP2<0.2':>11}")
for name in sorted(data):
    gh1 = sum(v[0] for v in data[name].values())
    gh2 = sum(v[1] for v in data[name].values())
    gr = gh1/max(gh2,1)
    valid = [(c, data[name].get(c,(0,0))) for c in chrs if sum(data[name].get(c,(0,0)))>100]
    rs = [v[0]/max(v[1],1) for _,v in valid]
    n_hp1 = sum(1 for r in rs if r > 5)
    n_bal = sum(1 for r in rs if 0.5 <= r <= 2)
    n_hp2 = sum(1 for r in rs if r < 0.2)
    print(f"{name:<22} {gh1:>11} {gh2:>10} {gr:>7.3f} {n_hp1:>9}/{len(valid)} {n_bal:>10}/{len(valid)} {n_hp2:>10}/{len(valid)}")
PYEOF
