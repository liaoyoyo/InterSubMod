#!/usr/bin/env python3
"""
Step 3: HP:Z tag coverage check on 6 candidates.
For each candidate position, count paired BAM reads by HP:Z tag.
Required: HP1 >= 10, HP2 >= 10, HPn-1 >= 5
"""
import subprocess, sys, json, os
from collections import Counter

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
TSV = f"{BASE}/output/tp_in_tsg_promoter_nonLOH.tsv"
OUT = f"{BASE}/output"

# Read candidates
candidates = []
with open(TSV) as f:
    next(f)  # header
    for line in f:
        fields = line.rstrip("\n").split("\t")
        candidates.append({
            "chrom": fields[0],
            "pos": int(fields[1]),
            "ref": fields[2],
            "alt": fields[3],
            "tvaf": float(fields[4]),
            "nvaf": float(fields[5]),
            "gene": fields[6],
        })

print(f"[INFO] Processing {len(candidates)} candidates")

results = []
for c in candidates:
    region = f"{c['chrom']}:{c['pos']}-{c['pos']}"
    # samtools view to extract reads overlapping the position
    p = subprocess.run(
        ["samtools", "view", BAM, region],
        capture_output=True, text=True
    )
    if p.returncode != 0:
        print(f"[ERROR] samtools view failed for {region}: {p.stderr}")
        continue

    counts = Counter()
    counts["total"] = 0
    counts["primary"] = 0
    for line in p.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        counts["total"] += 1
        flag = int(fields[1])
        # Skip secondary + supplementary + duplicates
        if flag & 0x100 or flag & 0x800 or flag & 0x400:
            continue
        counts["primary"] += 1
        # Find HP tag
        hp = None
        for tag in fields[11:]:
            if tag.startswith("HP:Z:"):
                hp = tag[5:]
                break
        if hp is None:
            counts["untag"] += 1
        else:
            counts[f"HP:{hp}"] += 1

    c["coverage"] = dict(counts)

    # Compute summary
    hp1 = counts.get("HP:1", 0)
    hp2 = counts.get("HP:2", 0)
    hp1_1 = counts.get("HP:1-1", 0)
    hp2_1 = counts.get("HP:2-1", 0)
    hp3 = counts.get("HP:3", 0)
    hpn_1 = hp1_1 + hp2_1
    untag = counts.get("untag", 0)

    c["hp_summary"] = {
        "HP1": hp1, "HP2": hp2,
        "HP1-1": hp1_1, "HP2-1": hp2_1,
        "HPn-1_total": hpn_1,
        "HP3": hp3, "untag": untag,
        "primary": counts["primary"], "total": counts["total"],
    }

    # Pass criteria
    c["pass_strict"] = (hp1 >= 10 and hp2 >= 10 and hpn_1 >= 5)
    c["pass_relaxed"] = (hp1 >= 5 and hp2 >= 5 and hpn_1 >= 3)

    print(f"\n{c['gene']} {c['chrom']}:{c['pos']} {c['ref']}>{c['alt']} TVAF={c['tvaf']}")
    print(f"  primary={counts['primary']}  HP1={hp1} HP2={hp2}  HP1-1={hp1_1} HP2-1={hp2_1}  HP3={hp3}  untag={untag}")
    print(f"  pass_strict (HP1>=10 & HP2>=10 & HPn-1>=5): {c['pass_strict']}")
    print(f"  pass_relaxed (HP1>=5 & HP2>=5 & HPn-1>=3):  {c['pass_relaxed']}")

# Save results
with open(f"{OUT}/step3_hp_coverage.json", "w") as f:
    json.dump(candidates, f, indent=2)

# Summary
print(f"\n{'='*60}")
print(f"[SUMMARY]")
n_strict = sum(1 for c in candidates if c["pass_strict"])
n_relaxed = sum(1 for c in candidates if c["pass_relaxed"])
print(f"  Strict pass (HP1>=10, HP2>=10, HPn-1>=5): {n_strict}/{len(candidates)}")
print(f"  Relaxed pass (HP1>=5, HP2>=5, HPn-1>=3):  {n_relaxed}/{len(candidates)}")
print(f"\n  Strict passers: {[c['gene'] for c in candidates if c['pass_strict']]}")
print(f"  Relaxed passers: {[c['gene'] for c in candidates if c['pass_relaxed']]}")
