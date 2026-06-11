#!/usr/bin/env python3
"""
Negative control: pick 5 random non-TSG TP sSNV positions (matched for TVAF + chr distribution),
run same ISM differential pipeline. If random regions also show |delta| > 0.05, then BRCA2 result is artifact.
"""
import json, os, random, subprocess
from collections import defaultdict
import pysam, numpy as np
from scipy import stats

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
TP_VCF = "/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
TSG_PROMOTER = f"{BASE}/output/tsg_promoter_merged.bed"
LOH_BED = f"{BASE}/output/seqc2_loh_only.bed"
OUT = f"{BASE}/output"

random.seed(42)
WINDOW_BP = 1000
ML_METH_THR = 200
ML_UNMETH_THR = 50

# Step 1: get non-TSG non-LOH TP sSNV with TVAF in [0.05, 0.5] range matching our positives
print("[INFO] Building random control candidate set...")

# Filter VCF: non-TSG, non-LOH, TVAF in range
control_bed = f"{OUT}/control_candidates_raw.bed"
subprocess.run(
    f"bcftools view -f PASS {TP_VCF} 2>/dev/null | "
    f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/TVAF\\n' | "
    f"awk 'BEGIN{{OFS=\"\\t\"}} $5+0 >= 0.05 && $5+0 <= 0.5 {{print $1, $2-1, $2, $3, $4, $5}}' "
    f"> {OUT}/all_tp_filtered_tvaf.bed",
    shell=True, check=True
)
# Remove TSG promoter overlap
subprocess.run(
    f"bedtools intersect -a {OUT}/all_tp_filtered_tvaf.bed -b {TSG_PROMOTER} -v > {OUT}/tp_nontsg.bed",
    shell=True, check=True
)
# Remove LOH overlap
subprocess.run(
    f"bedtools intersect -a {OUT}/tp_nontsg.bed -b {LOH_BED} -v > {OUT}/tp_nontsg_nonloh.bed",
    shell=True, check=True
)

with open(f"{OUT}/tp_nontsg_nonloh.bed") as f:
    candidates_all = [l.rstrip().split("\t") for l in f]
print(f"[INFO] Non-TSG non-LOH TP sSNV pool: {len(candidates_all)}")

# Randomly sample 5
sampled = random.sample(candidates_all, 5)
print(f"[INFO] Sampled 5 random controls:")
for s in sampled:
    print(f"  {s[0]}:{s[2]} TVAF={s[5]}")

bam = pysam.AlignmentFile(BAM, "rb")

def extract_hp(read):
    for tag, val in read.tags:
        if tag == "HP":
            return str(val)
    return None

def get_methylation_calls(read):
    mods = {}
    try:
        modified = read.modified_bases
    except Exception:
        return mods
    if not modified:
        return mods
    aligned_pairs = read.get_aligned_pairs(matches_only=False)
    rd_to_ref = {rd: rf for rd, rf in aligned_pairs if rd is not None and rf is not None}
    for key, calls in modified.items():
        base, strand, mod_code = key
        if mod_code != 'm':
            continue
        for read_pos, ml in calls:
            ref_pos = rd_to_ref.get(read_pos)
            if ref_pos is None:
                continue
            mods[ref_pos] = ml
    return mods


def analyze_pos(chrom, var_pos, tvaf):
    start = max(0, var_pos - 1 - WINDOW_BP)
    end = var_pos - 1 + WINDOW_BP

    cpg_hp_counts = defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))
    read_count = defaultdict(int)
    for read in bam.fetch(chrom, start, end):
        if read.flag & 0x100 or read.flag & 0x800 or read.flag & 0x400:
            continue
        hp = extract_hp(read) or "untag"
        read_count[hp] += 1
        mods = get_methylation_calls(read)
        for ref_pos, ml in mods.items():
            if ref_pos < start or ref_pos > end:
                continue
            if ml >= ML_METH_THR:
                cpg_hp_counts[ref_pos][hp][0] += 1
            elif ml <= ML_UNMETH_THR:
                cpg_hp_counts[ref_pos][hp][1] += 1
            else:
                cpg_hp_counts[ref_pos][hp][2] += 1

    # Pair HP1 vs HP1-1 (must have both)
    paired_g, paired_s = [], []
    for cpg, hp_dict in cpg_hp_counts.items():
        ng_m, ng_u, _ = hp_dict.get("1", [0,0,0])
        ns_m, ns_u, _ = hp_dict.get("1-1", [0,0,0])
        ng = ng_m + ng_u
        ns = ns_m + ns_u
        if ng >= 3 and ns >= 3:
            paired_g.append(ng_m / ng)
            paired_s.append(ns_m / ns)
    return {
        "n_paired": len(paired_g),
        "mean_g": float(np.mean(paired_g)) if paired_g else None,
        "mean_s": float(np.mean(paired_s)) if paired_s else None,
        "delta": float(np.mean(paired_s) - np.mean(paired_g)) if paired_g else None,
        "w_p": float(stats.wilcoxon(paired_g, paired_s).pvalue) if len(paired_g) >= 5 else None,
        "read_count": dict(read_count),
    }

print(f"\n[INFO] Running ISM analysis on 5 random controls...")
controls = []
for s in sampled:
    chrom = s[0]
    pos = int(s[2])
    tvaf = float(s[5])
    res = analyze_pos(chrom, pos, tvaf)
    res["chrom"] = chrom
    res["pos"] = pos
    res["tvaf"] = tvaf
    controls.append(res)
    print(f"\n  {chrom}:{pos} TVAF={tvaf}")
    print(f"    n_paired_cpg={res['n_paired']}  reads={res['read_count']}")
    if res["n_paired"] >= 5:
        print(f"    germ_beta={res['mean_g']:.3f}  som_beta={res['mean_s']:.3f}  delta={res['delta']:+.3f}  p={res['w_p']:.2e}")
    else:
        print(f"    insufficient paired CpG")

# Summary stats
deltas = [c["delta"] for c in controls if c["delta"] is not None]
n_signif = sum(1 for c in controls if c["w_p"] is not None and c["w_p"] < 0.05 and abs(c["delta"] or 0) >= 0.05)
print(f"\n{'='*60}")
print(f"[CONTROL SUMMARY]")
print(f"  Random control regions analyzed: {len(controls)}")
print(f"  With sufficient paired CpG: {sum(1 for c in controls if c['n_paired'] >= 5)}")
if deltas:
    print(f"  Mean abs |delta|: {np.mean(np.abs(deltas)):.3f}")
    print(f"  Max abs |delta|: {np.max(np.abs(deltas)):.3f}")
    print(f"  Regions with p<0.05 AND |delta|>=0.05: {n_signif}/{len(controls)}")
    print(f"\n  BRCA2 |delta| = 0.122 — is this within null distribution?")
    if abs(0.122) > np.max(np.abs(deltas)):
        print(f"  -> BRCA2 result EXCEEDS max control delta. Real signal.")
    else:
        print(f"  -> BRCA2 result IS WITHIN null. May be artifact.")

with open(f"{OUT}/step4_negative_control.json", "w") as f:
    json.dump(controls, f, indent=2)
print(f"\nSaved -> {OUT}/step4_negative_control.json")
