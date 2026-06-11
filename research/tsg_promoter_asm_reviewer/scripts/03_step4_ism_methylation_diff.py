#!/usr/bin/env python3
"""
Step 4: ISM per-haplotype methylation differential.

For each strict-pass candidate (BRCA2 + KMT2C):
  - Define +/- 1kb window around variant
  - Fetch reads with MM/ML mod tags + HP:Z
  - Parse 5mC modification probabilities
  - Per-CpG per-HP-group methylation rate (beta)
  - Compare HP1 vs HP1-1 (germline-HP1 vs somatic-reconstructed-on-HP1)
                HP2 vs HP2-1 (symmetrical)
  - Stat: per-CpG Fisher exact + window-level Mann-Whitney on betas
  - Effect size: |delta_beta| mean, max
  - Caveat: also stratify by AF, by CpG position to AF site

ML threshold: 256-bin probability, >=200 (~0.78) = methylated, <=50 (~0.20) = unmethylated, else ambiguous
"""
import json, os, sys
from collections import defaultdict, Counter
import pysam
import numpy as np
from scipy import stats

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
HP_JSON = f"{BASE}/output/step3_hp_coverage.json"
OUT = f"{BASE}/output"

WINDOW_BP = 1000
ML_METH_THR = 200   # >=200 (78%) = methylated
ML_UNMETH_THR = 50  # <=50 (20%) = unmethylated

with open(HP_JSON) as f:
    candidates = json.load(f)

passers = [c for c in candidates if c.get("pass_strict")]
print(f"[INFO] Processing {len(passers)} strict-pass candidates: {[c['gene'] for c in passers]}")

bam = pysam.AlignmentFile(BAM, "rb")

def extract_hp(read):
    for tag, val in read.tags:
        if tag == "HP":
            return str(val)
    return None

def get_methylation_calls(read):
    """Return dict {ref_pos: ml_score} for 5mCG modifications on this read."""
    mods = {}
    try:
        modified = read.modified_bases  # {(base, strand, mod): [(read_pos, ml_score), ...]}
    except Exception:
        return mods
    if not modified:
        return mods

    # Build read_pos -> ref_pos map
    aligned_pairs = read.get_aligned_pairs(matches_only=False)
    rd_to_ref = {rd: rf for rd, rf in aligned_pairs if rd is not None and rf is not None}

    for key, calls in modified.items():
        base, strand, mod_code = key
        # We want 5mC: mod_code 'm' or 'C+m' or 'm'
        if mod_code != 'm':
            continue
        for read_pos, ml in calls:
            ref_pos = rd_to_ref.get(read_pos)
            if ref_pos is None:
                continue
            mods[ref_pos] = ml
    return mods


def analyze_candidate(c):
    chrom = c["chrom"]
    var_pos = c["pos"]   # 1-based
    gene = c["gene"]
    start = max(0, var_pos - 1 - WINDOW_BP)
    end = var_pos - 1 + WINDOW_BP

    print(f"\n{'='*60}")
    print(f"{gene} {chrom}:{var_pos} {c['ref']}>{c['alt']} TVAF={c['tvaf']}  Window ±{WINDOW_BP}bp")

    # per-CpG-position per-HP-group counts: cpg_pos -> hp_group -> [meth, unmeth, amb]
    cpg_hp_counts = defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))
    read_count = defaultdict(int)

    for read in bam.fetch(chrom, start, end):
        if read.flag & 0x100 or read.flag & 0x800 or read.flag & 0x400:
            continue
        hp = extract_hp(read)
        if hp is None:
            hp = "untag"
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

    print(f"  Reads in window by HP: {dict(read_count)}")
    print(f"  Unique CpG positions with calls: {len(cpg_hp_counts)}")

    # Per-CpG beta per group
    groups = ["1", "2", "1-1", "2-1", "3", "untag"]
    cpg_records = []
    for cpg, hp_dict in sorted(cpg_hp_counts.items()):
        rec = {"cpg_pos": cpg, "dist_to_var": cpg - (var_pos - 1)}
        for g in groups:
            m, u, a = hp_dict.get(g, [0, 0, 0])
            tot = m + u
            rec[f"HP{g}_meth"] = m
            rec[f"HP{g}_unmeth"] = u
            rec[f"HP{g}_amb"] = a
            rec[f"HP{g}_beta"] = (m / tot) if tot > 0 else None
            rec[f"HP{g}_n"] = tot
        cpg_records.append(rec)

    # Group-level comparison: collect betas per group, then Mann-Whitney
    def collect_betas(g, min_n=3):
        return [r[f"HP{g}_beta"] for r in cpg_records if r[f"HP{g}_n"] >= min_n and r[f"HP{g}_beta"] is not None]

    beta_HP1 = collect_betas("1")
    beta_HP2 = collect_betas("2")
    beta_HP1_1 = collect_betas("1-1")
    beta_HP2_1 = collect_betas("2-1")

    def fmt(x):
        return f"{np.mean(x):.3f}" if x else "NA"
    print(f"  CpGs with >=3 reads:  HP1={len(beta_HP1)}  HP2={len(beta_HP2)}  HP1-1={len(beta_HP1_1)}  HP2-1={len(beta_HP2_1)}")
    print(f"  Mean beta:  HP1={fmt(beta_HP1)}  HP2={fmt(beta_HP2)}  HP1-1={fmt(beta_HP1_1)}  HP2-1={fmt(beta_HP2_1)}")

    # Critical comparison: HP1 vs HP1-1 (same germline allele, germline-tag vs somatic-reconstructed-tag)
    stats_results = {}
    for (g_germ, g_som) in [("1", "1-1"), ("2", "2-1")]:
        betas_g = collect_betas(g_germ)
        betas_s = collect_betas(g_som)
        if len(betas_g) >= 5 and len(betas_s) >= 5:
            # Paired comparison: same CpG with both groups
            paired_g, paired_s = [], []
            for r in cpg_records:
                bg, bs = r.get(f"HP{g_germ}_beta"), r.get(f"HP{g_som}_beta")
                ng, ns = r.get(f"HP{g_germ}_n", 0), r.get(f"HP{g_som}_n", 0)
                if bg is not None and bs is not None and ng >= 3 and ns >= 3:
                    paired_g.append(bg)
                    paired_s.append(bs)

            if len(paired_g) >= 5:
                # Wilcoxon signed-rank on paired CpGs
                try:
                    w_stat, w_p = stats.wilcoxon(paired_g, paired_s, zero_method="wilcox")
                except Exception as e:
                    w_stat, w_p = None, None
                # Mean Delta beta
                deltas = [s - g for g, s in zip(paired_g, paired_s)]
                mean_delta = float(np.mean(deltas))
                max_abs_delta = float(np.max(np.abs(deltas)))
                # Cohen's d
                pooled_sd = float(np.sqrt((np.var(paired_g, ddof=1) + np.var(paired_s, ddof=1)) / 2))
                cohen_d = mean_delta / pooled_sd if pooled_sd > 0 else None
                # Unpaired Mann-Whitney as sanity
                mw_stat, mw_p = stats.mannwhitneyu(betas_g, betas_s, alternative="two-sided")
                stats_results[f"HP{g_germ}_vs_HP{g_som}"] = {
                    "n_paired_cpg": len(paired_g),
                    "mean_germ_beta": float(np.mean(paired_g)),
                    "mean_som_beta": float(np.mean(paired_s)),
                    "mean_delta_som_minus_germ": mean_delta,
                    "max_abs_delta": max_abs_delta,
                    "cohen_d": cohen_d,
                    "wilcoxon_p": float(w_p) if w_p is not None else None,
                    "mannwhitney_p": float(mw_p),
                    "interpretation": (
                        "somatic HYPOMETHYLATED vs germline" if mean_delta < -0.05 else
                        "somatic HYPERMETHYLATED vs germline" if mean_delta > 0.05 else
                        "no significant difference"
                    ),
                }
                print(f"\n  -- HP{g_germ} vs HP{g_som} (germline vs somatic-reconstructed, same allele) --")
                print(f"     n paired CpG: {len(paired_g)}")
                print(f"     mean beta germ={np.mean(paired_g):.3f}  som={np.mean(paired_s):.3f}  delta={mean_delta:+.3f}")
                print(f"     max |delta|: {max_abs_delta:.3f}")
                print(f"     Cohen's d: {cohen_d:.3f}" if cohen_d is not None else "     Cohen's d: NA")
                print(f"     Wilcoxon p: {w_p:.3e}" if w_p is not None else "     Wilcoxon p: NA")
                print(f"     Mann-Whitney p: {mw_p:.3e}")
                print(f"     Interpretation: {stats_results[f'HP{g_germ}_vs_HP{g_som}']['interpretation']}")
            else:
                print(f"\n  -- HP{g_germ} vs HP{g_som}: insufficient paired CpG (n={len(paired_g)}, need >=5) --")
        else:
            print(f"\n  -- HP{g_germ} vs HP{g_som}: insufficient CpG-level coverage (germ={len(betas_g)}, som={len(betas_s)}) --")

    return {
        "candidate": c,
        "n_reads_by_hp": dict(read_count),
        "n_cpg_positions": len(cpg_hp_counts),
        "cpg_records": cpg_records,
        "stats": stats_results,
    }


all_results = []
for c in passers:
    res = analyze_candidate(c)
    all_results.append(res)

# Save
with open(f"{OUT}/step4_ism_results.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)

# Summary
print(f"\n{'='*60}")
print(f"[FINAL SUMMARY]")
for r in all_results:
    gene = r["candidate"]["gene"]
    print(f"\n  {gene}:")
    for key, s in r["stats"].items():
        if "cohen_d" in s and s["cohen_d"] is not None:
            actionable = "✓ ACTIONABLE" if abs(s["cohen_d"]) >= 0.5 or abs(s["mean_delta_som_minus_germ"]) >= 0.2 else "○ marginal"
            print(f"    {key}: delta={s['mean_delta_som_minus_germ']:+.3f}, d={s['cohen_d']:+.2f}, p_wilcox={s['wilcoxon_p']:.2e} {actionable}")
print(f"\nSaved -> {OUT}/step4_ism_results.json")
