#!/usr/bin/env python3
"""
觀察 O14: LOH 區域內 haplotag 偏向性分析
---
問題：在 SEQC2 確認的 TP LOH 區域內，TO haplotagged reads 的 HP tag 是否
      正確偏向同一個 haplotype？

分析角度：
1. HP_Ratio 分佈（TP LOH vs TN 非 LOH vs FP/FN）
2. HP1FamilyN vs HP2FamilyN 散佈圖（偏向性）
3. dominant haplotype 一致性（同一 LOH region 內的 sites 是否偏向同一 HP？）
4. HP0 reads 比例（未分配 haplotype 的 reads 比例）
5. TO vs Paired 對比
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from pathlib import Path
import gzip
import warnings
warnings.filterwarnings('ignore')

# --- paths ---
MASTER = Path("/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
              "20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz")
LOH_ZONE = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/loh_investigation/"
                "data/seqc2_vs_to_validation/hcc1395_variant_loh_zone.tsv")
LOH_TP_BED = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/loh_investigation/"
                  "data/seqc2_vs_to_validation/loh_classified_TP.bed")
LOH_FP_BED = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/loh_investigation/"
                  "data/seqc2_vs_to_validation/loh_classified_FP.bed")
LOH_FN_BED = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/loh_investigation/"
                  "data/seqc2_vs_to_validation/loh_classified_FN.bed")
OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/loh_investigation/figures")
DATA_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/loh_investigation/data")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# --- helpers ---
def load_bed_intervals(bed_path):
    """Load BED file into list of (chr, start, end)."""
    intervals = {}
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            intervals.setdefault(chrom, []).append((start, end))
    # sort by start
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals

def point_in_intervals(chrom, pos, intervals):
    """Check if a point falls in any interval (binary search)."""
    if chrom not in intervals:
        return False
    ivs = intervals[chrom]
    lo, hi = 0, len(ivs) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        s, e = ivs[mid]
        if pos < s:
            hi = mid - 1
        elif pos >= e:
            lo = mid + 1
        else:
            return True
    return False

def assign_loh_zone(chrom, pos, tp_ivs, fp_ivs, fn_ivs):
    """Assign LOH zone label."""
    if point_in_intervals(chrom, pos, tp_ivs):
        return 'TP'
    elif point_in_intervals(chrom, pos, fp_ivs):
        return 'FP'
    elif point_in_intervals(chrom, pos, fn_ivs):
        return 'FN'
    else:
        return 'TN'

# --- load master dataset (HCC1395 only) ---
print("Loading master dataset...")
cols_needed = ['Chr', 'Pos', 'HP1FamilyN', 'HP2FamilyN', 'HP_Ratio', 'Potential_LOH',
               'NHP0', 'effective_hp_reads', 'truth_label', 'sample', 'mode',
               'VerificationClass', 'DominantLabel']
df_all = pd.read_csv(MASTER, sep='\t', usecols=cols_needed)
print(f"  Total rows: {len(df_all):,}")

# Filter HCC1395 only (both modes)
df_hcc = df_all[df_all['sample'].isin(['HCC1395', 'HCC1395_5kHz'])].copy()
print(f"  HCC1395 rows: {len(df_hcc):,}")

# --- assign LOH zones ---
print("Loading LOH zone BED files...")
tp_ivs = load_bed_intervals(LOH_TP_BED)
fp_ivs = load_bed_intervals(LOH_FP_BED)
fn_ivs = load_bed_intervals(LOH_FN_BED)

print("Assigning LOH zones...")
df_hcc['loh_zone'] = df_hcc.apply(
    lambda r: assign_loh_zone(r['Chr'], r['Pos'], tp_ivs, fp_ivs, fn_ivs), axis=1)

# Derived columns
df_hcc['hp_max'] = df_hcc[['HP1FamilyN', 'HP2FamilyN']].max(axis=1)
df_hcc['hp_min'] = df_hcc[['HP1FamilyN', 'HP2FamilyN']].min(axis=1)
df_hcc['hp_total'] = df_hcc['HP1FamilyN'] + df_hcc['HP2FamilyN']
df_hcc['hp_bias'] = np.where(df_hcc['hp_total'] > 0,
                             df_hcc['hp_max'] / df_hcc['hp_total'], 0.5)
df_hcc['dominant_hp'] = np.where(df_hcc['HP1FamilyN'] >= df_hcc['HP2FamilyN'], 'HP1', 'HP2')
df_hcc['hp0_ratio'] = np.where(
    df_hcc['hp_total'] + df_hcc['NHP0'] > 0,
    df_hcc['NHP0'] / (df_hcc['hp_total'] + df_hcc['NHP0']),
    0)

# Split by mode
df_to = df_hcc[df_hcc['mode'] == 'to'].copy()
df_paired = df_hcc[df_hcc['mode'] == 'paired'].copy()
print(f"  TO: {len(df_to):,}, Paired: {len(df_paired):,}")

for mode_name, dft in [('TO', df_to), ('Paired', df_paired)]:
    print(f"\n  {mode_name} LOH zone distribution:")
    for zone in ['TP', 'FP', 'FN', 'TN']:
        n = (dft['loh_zone'] == zone).sum()
        print(f"    {zone}: {n:,}")

# ============================================================================
# Fig01: HP_Ratio distribution by LOH zone (TO vs Paired)
# ============================================================================
print("\n--- Fig01: HP_Ratio distribution by LOH zone ---")
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

for ax, (mode_name, dft) in zip(axes, [('TO', df_to), ('Paired', df_paired)]):
    for zone, color in [('TP', '#2ca02c'), ('TN', '#7f7f7f'), ('FP', '#d62728'), ('FN', '#ff7f0e')]:
        subset = dft[dft['loh_zone'] == zone]['HP_Ratio'].dropna()
        if len(subset) > 10:
            ax.hist(subset, bins=50, range=(0, 1), alpha=0.5, label=f'{zone} (n={len(subset):,})',
                    color=color, density=True)
    ax.set_xlabel('HP_Ratio (HP1 / (HP1+HP2))')
    ax.set_ylabel('Density')
    ax.set_title(f'{mode_name} — HP_Ratio by LOH Zone')
    ax.legend(fontsize=8)
    ax.axvline(0.1, color='red', linestyle='--', alpha=0.3, linewidth=0.8)
    ax.axvline(0.9, color='red', linestyle='--', alpha=0.3, linewidth=0.8)

plt.suptitle('Fig01: HP_Ratio Distribution — TP LOH vs TN (SEQC2 zones)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(OUT_DIR / 'o14_fig01_hp_ratio_by_loh_zone.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved fig01")

# --- Print summary stats ---
print("\n  HP_Ratio summary (TO):")
for zone in ['TP', 'FP', 'FN', 'TN']:
    sub = df_to[df_to['loh_zone'] == zone]['HP_Ratio'].dropna()
    if len(sub) > 0:
        extreme = ((sub < 0.1) | (sub > 0.9)).mean() * 100
        print(f"    {zone}: n={len(sub):,}, median={sub.median():.3f}, "
              f"mean={sub.mean():.3f}, extreme(<0.1|>0.9)={extreme:.1f}%")

print("\n  HP_Ratio summary (Paired):")
for zone in ['TP', 'FP', 'FN', 'TN']:
    sub = df_paired[df_paired['loh_zone'] == zone]['HP_Ratio'].dropna()
    if len(sub) > 0:
        extreme = ((sub < 0.1) | (sub > 0.9)).mean() * 100
        print(f"    {zone}: n={len(sub):,}, median={sub.median():.3f}, "
              f"mean={sub.mean():.3f}, extreme(<0.1|>0.9)={extreme:.1f}%")

# ============================================================================
# Fig02: HP1 vs HP2 scatter — TP LOH vs TN (TO mode)
# ============================================================================
print("\n--- Fig02: HP1 vs HP2 scatter ---")
fig, axes = plt.subplots(1, 4, figsize=(22, 5))

for ax, zone in zip(axes, ['TP', 'FP', 'FN', 'TN']):
    sub = df_to[df_to['loh_zone'] == zone].copy()
    if len(sub) > 5000:
        sub = sub.sample(5000, random_state=42)
    if len(sub) > 0:
        colors_map = sub['truth_label'].map({'TP': '#2ca02c', 'FP': '#d62728'}).fillna('#999999')
        ax.scatter(sub['HP1FamilyN'], sub['HP2FamilyN'], alpha=0.15, s=5, c=colors_map)
        ax.plot([0, sub[['HP1FamilyN', 'HP2FamilyN']].max().max()],
                [0, sub[['HP1FamilyN', 'HP2FamilyN']].max().max()],
                'k--', alpha=0.3, linewidth=0.8)
    ax.set_xlabel('HP1 Family N')
    ax.set_ylabel('HP2 Family N')
    ax.set_title(f'{zone} zone (n={len(df_to[df_to["loh_zone"]==zone]):,})')
    ax.set_aspect('equal')
    max_val = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.set_xlim(0, min(max_val, 300))
    ax.set_ylim(0, min(max_val, 300))

plt.suptitle('Fig02: HP1 vs HP2 Family N — TO mode (green=TP, red=FP)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(OUT_DIR / 'o14_fig02_hp1_vs_hp2_scatter.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved fig02")

# ============================================================================
# Fig03: HP bias (max/total) distribution by zone — TO vs Paired
# ============================================================================
print("\n--- Fig03: HP bias distribution ---")
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

for ax, (mode_name, dft) in zip(axes, [('TO', df_to), ('Paired', df_paired)]):
    # filter hp_total > 0
    dft_valid = dft[dft['hp_total'] > 0]
    for zone, color in [('TP', '#2ca02c'), ('TN', '#7f7f7f')]:
        subset = dft_valid[dft_valid['loh_zone'] == zone]['hp_bias'].dropna()
        if len(subset) > 10:
            ax.hist(subset, bins=50, range=(0.5, 1.0), alpha=0.6,
                    label=f'{zone} (n={len(subset):,}, median={subset.median():.3f})',
                    color=color, density=True)
    ax.set_xlabel('HP Bias = max(HP1,HP2) / (HP1+HP2)')
    ax.set_ylabel('Density')
    ax.set_title(f'{mode_name} — Haplotype Bias')
    ax.legend(fontsize=9)
    ax.axvline(0.9, color='red', linestyle='--', alpha=0.4, linewidth=1)

plt.suptitle('Fig03: Haplotype Bias — TP LOH vs TN non-LOH\n'
             'Bias=1.0 means all reads on one HP; Bias=0.5 means even split',
             fontsize=13, y=1.05)
plt.tight_layout()
plt.savefig(OUT_DIR / 'o14_fig03_hp_bias_distribution.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved fig03")

# Print bias stats
print("\n  HP Bias summary (hp_total > 0):")
for mode_name, dft in [('TO', df_to), ('Paired', df_paired)]:
    dft_valid = dft[dft['hp_total'] > 0]
    print(f"  {mode_name}:")
    for zone in ['TP', 'TN']:
        sub = dft_valid[dft_valid['loh_zone'] == zone]['hp_bias']
        if len(sub) > 0:
            pct_90 = (sub >= 0.9).mean() * 100
            pct_95 = (sub >= 0.95).mean() * 100
            print(f"    {zone}: median={sub.median():.3f}, ≥0.9={pct_90:.1f}%, ≥0.95={pct_95:.1f}%")

# ============================================================================
# Fig04: HP0 ratio by LOH zone (TO vs Paired)
# ============================================================================
print("\n--- Fig04: HP0 (untagged) ratio by zone ---")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, (mode_name, dft) in zip(axes, [('TO', df_to), ('Paired', df_paired)]):
    zone_data = []
    for zone in ['TP', 'FP', 'FN', 'TN']:
        sub = dft[(dft['loh_zone'] == zone) & (dft['hp_total'] + dft['NHP0'] > 0)]
        if len(sub) > 0:
            zone_data.append({
                'zone': zone,
                'hp0_ratio_median': sub['hp0_ratio'].median(),
                'hp0_ratio_mean': sub['hp0_ratio'].mean(),
                'n': len(sub)
            })
    if zone_data:
        zd = pd.DataFrame(zone_data)
        bars = ax.bar(zd['zone'], zd['hp0_ratio_median'],
                      color=['#2ca02c', '#d62728', '#ff7f0e', '#7f7f7f'])
        for i, row in zd.iterrows():
            ax.text(i, row['hp0_ratio_median'] + 0.005,
                    f"{row['hp0_ratio_median']:.3f}\nn={row['n']:,}", ha='center', fontsize=8)
    ax.set_ylabel('HP0 Ratio (median)')
    ax.set_xlabel('LOH Zone')
    ax.set_title(f'{mode_name} — HP0 (untagged) Read Ratio')
    ax.set_ylim(0, max(0.5, ax.get_ylim()[1] * 1.2))

plt.suptitle('Fig04: Untagged (HP0) Read Ratio by LOH Zone\n'
             'Low HP0 = most reads assigned to HP1 or HP2', fontsize=13, y=1.05)
plt.tight_layout()
plt.savefig(OUT_DIR / 'o14_fig04_hp0_ratio_by_zone.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved fig04")

# ============================================================================
# Fig05: Within-region dominant HP consistency
# Group sites by LOH TP region, check if all sites agree on dominant HP
# ============================================================================
print("\n--- Fig05: Dominant HP consistency within LOH TP regions ---")

# Assign each TO site in TP zone to a specific TP region
tp_regions = []
with open(LOH_TP_BED) as f:
    for line in f:
        parts = line.strip().split('\t')
        tp_regions.append((parts[0], int(parts[1]), int(parts[2])))

def find_region_idx(chrom, pos, regions):
    """Find which region a site belongs to."""
    for i, (c, s, e) in enumerate(regions):
        if c == chrom and s <= pos < e:
            return i
    return -1

# Only analyze TO TP zone sites with sufficient reads
to_tp = df_to[(df_to['loh_zone'] == 'TP') & (df_to['effective_hp_reads'] >= 10)].copy()
print(f"  TO TP sites with ehp≥10: {len(to_tp):,}")

# Sample if too many
if len(to_tp) > 50000:
    to_tp_sample = to_tp.sample(50000, random_state=42)
else:
    to_tp_sample = to_tp

# Assign region index
to_tp_sample['region_idx'] = to_tp_sample.apply(
    lambda r: find_region_idx(r['Chr'], r['Pos'], tp_regions), axis=1)
to_tp_sample = to_tp_sample[to_tp_sample['region_idx'] >= 0]
print(f"  Assigned to regions: {len(to_tp_sample):,}")

# Per-region dominant HP consistency
region_stats = []
for ridx, grp in to_tp_sample.groupby('region_idx'):
    if len(grp) < 3:
        continue
    hp1_dominant = (grp['dominant_hp'] == 'HP1').sum()
    hp2_dominant = (grp['dominant_hp'] == 'HP2').sum()
    total = len(grp)
    dominant_frac = max(hp1_dominant, hp2_dominant) / total
    region_chrom = tp_regions[ridx][0]
    region_size = tp_regions[ridx][2] - tp_regions[ridx][1]
    region_stats.append({
        'region_idx': ridx,
        'chrom': region_chrom,
        'size_mb': region_size / 1e6,
        'n_sites': total,
        'hp1_dominant': hp1_dominant,
        'hp2_dominant': hp2_dominant,
        'consistency': dominant_frac,
        'dominant_hp': 'HP1' if hp1_dominant >= hp2_dominant else 'HP2'
    })

rs = pd.DataFrame(region_stats)
print(f"  Regions with ≥3 sites: {len(rs)}")

if len(rs) > 0:
    print(f"  Consistency median: {rs['consistency'].median():.3f}")
    print(f"  Consistency ≥0.9: {(rs['consistency'] >= 0.9).mean()*100:.1f}%")
    print(f"  Consistency ≥0.8: {(rs['consistency'] >= 0.8).mean()*100:.1f}%")
    print(f"  Consistency ≥0.7: {(rs['consistency'] >= 0.7).mean()*100:.1f}%")

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Left: consistency histogram
    ax = axes[0]
    ax.hist(rs['consistency'], bins=30, range=(0.5, 1.0), color='#2ca02c', edgecolor='white', alpha=0.8)
    ax.axvline(rs['consistency'].median(), color='red', linestyle='--',
               label=f'median={rs["consistency"].median():.3f}')
    ax.set_xlabel('Dominant HP Consistency\n(fraction of sites agreeing on dominant HP)')
    ax.set_ylabel('Number of LOH Regions')
    ax.set_title(f'Within-Region HP Consistency (n={len(rs)} TP LOH regions)')
    ax.legend()

    # Right: consistency vs region size
    ax = axes[1]
    ax.scatter(rs['size_mb'], rs['consistency'], alpha=0.3, s=15, c='#2ca02c')
    ax.set_xlabel('LOH Region Size (Mb)')
    ax.set_ylabel('Dominant HP Consistency')
    ax.set_title('HP Consistency vs Region Size')
    ax.axhline(0.9, color='red', linestyle='--', alpha=0.4)
    ax.set_ylim(0.45, 1.05)

    plt.suptitle('Fig05: Within-Region Dominant HP Consistency — TO TP LOH Regions\n'
                 'Q: Do all sites within one LOH region agree on which HP is dominant?',
                 fontsize=13, y=1.05)
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'o14_fig05_within_region_hp_consistency.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved fig05")

# ============================================================================
# Fig06: TO vs Paired HP_Ratio in same LOH TP zone — site-matched
# ============================================================================
print("\n--- Fig06: TO vs Paired HP_Ratio (site-matched, TP zone) ---")

# Merge TO and Paired on (Chr, Pos, sample) — but sample names differ
# Use the same positions — HCC1395 only
to_for_match = df_to[df_to['loh_zone'] == 'TP'][['Chr', 'Pos', 'HP_Ratio', 'hp_bias',
                                                    'HP1FamilyN', 'HP2FamilyN', 'truth_label']].copy()
to_for_match = to_for_match.rename(columns={
    'HP_Ratio': 'to_hp_ratio', 'hp_bias': 'to_hp_bias',
    'HP1FamilyN': 'to_hp1', 'HP2FamilyN': 'to_hp2', 'truth_label': 'to_label'
})

paired_for_match = df_paired[df_paired['loh_zone'] == 'TP'][['Chr', 'Pos', 'HP_Ratio', 'hp_bias',
                                                               'HP1FamilyN', 'HP2FamilyN', 'truth_label']].copy()
paired_for_match = paired_for_match.rename(columns={
    'HP_Ratio': 'paired_hp_ratio', 'hp_bias': 'paired_hp_bias',
    'HP1FamilyN': 'paired_hp1', 'HP2FamilyN': 'paired_hp2', 'truth_label': 'paired_label'
})

matched = pd.merge(to_for_match, paired_for_match, on=['Chr', 'Pos'], how='inner')
print(f"  Matched sites in TP LOH zone: {len(matched):,}")

if len(matched) > 0:
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # Left: TO HP_Ratio vs Paired HP_Ratio
    ax = axes[0]
    sample_n = min(len(matched), 10000)
    ms = matched.sample(sample_n, random_state=42) if len(matched) > sample_n else matched
    ax.scatter(ms['paired_hp_ratio'], ms['to_hp_ratio'], alpha=0.1, s=3, c='#2ca02c')
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)
    corr = matched[['to_hp_ratio', 'paired_hp_ratio']].corr().iloc[0, 1]
    ax.set_xlabel('Paired HP_Ratio')
    ax.set_ylabel('TO HP_Ratio')
    ax.set_title(f'HP_Ratio: TO vs Paired (r={corr:.3f})')

    # Middle: HP bias comparison
    ax = axes[1]
    ax.hist(matched['to_hp_bias'], bins=50, range=(0.5, 1), alpha=0.5,
            label=f'TO (median={matched["to_hp_bias"].median():.3f})', color='#d62728', density=True)
    ax.hist(matched['paired_hp_bias'], bins=50, range=(0.5, 1), alpha=0.5,
            label=f'Paired (median={matched["paired_hp_bias"].median():.3f})', color='#1f77b4', density=True)
    ax.set_xlabel('HP Bias = max(HP1,HP2)/(HP1+HP2)')
    ax.set_ylabel('Density')
    ax.set_title('HP Bias in TP LOH Zone')
    ax.legend()

    # Right: agreement on dominant HP direction
    ax = axes[2]
    # Check if TO and Paired agree on which HP is dominant
    matched['to_dom'] = np.where(matched['to_hp1'] >= matched['to_hp2'], 'HP1', 'HP2')
    matched['paired_dom'] = np.where(matched['paired_hp1'] >= matched['paired_hp2'], 'HP1', 'HP2')
    agree = (matched['to_dom'] == matched['paired_dom']).mean() * 100

    # By ehp tier
    matched['to_total'] = matched['to_hp1'] + matched['to_hp2']
    matched['paired_total'] = matched['paired_hp1'] + matched['paired_hp2']
    matched['min_total'] = matched[['to_total', 'paired_total']].min(axis=1)

    tiers = [(0, 10), (10, 30), (30, 50), (50, 100), (100, 200), (200, 9999)]
    tier_labels = []
    tier_agree = []
    tier_n = []
    for lo, hi in tiers:
        sub = matched[(matched['min_total'] >= lo) & (matched['min_total'] < hi)]
        if len(sub) > 0:
            a = (sub['to_dom'] == sub['paired_dom']).mean() * 100
            tier_labels.append(f'{lo}-{hi}' if hi < 9999 else f'{lo}+')
            tier_agree.append(a)
            tier_n.append(len(sub))

    bars = ax.bar(range(len(tier_labels)), tier_agree, color='#2ca02c', edgecolor='white')
    ax.set_xticks(range(len(tier_labels)))
    ax.set_xticklabels(tier_labels, fontsize=9)
    ax.set_xlabel('min(TO, Paired) total HP reads')
    ax.set_ylabel('TO-Paired Agreement (%)')
    ax.set_title(f'Dominant HP Direction Agreement\n(overall={agree:.1f}%)')
    ax.set_ylim(0, 105)
    for i, (a, n) in enumerate(zip(tier_agree, tier_n)):
        ax.text(i, a + 1.5, f'{a:.0f}%\nn={n:,}', ha='center', fontsize=7)

    plt.suptitle('Fig06: TO vs Paired Haplotag in TP LOH Zone — Same Sites\n'
                 'Q: Do TO and Paired agree on which haplotype is dominant?',
                 fontsize=13, y=1.05)
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'o14_fig06_to_vs_paired_hp_matched.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved fig06 (overall TO-Paired agreement: {agree:.1f}%)")

# ============================================================================
# Fig07: Effective HP reads distribution in LOH vs non-LOH
# ============================================================================
print("\n--- Fig07: Effective HP reads in LOH vs non-LOH ---")
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

for ax, (mode_name, dft) in zip(axes, [('TO', df_to), ('Paired', df_paired)]):
    for zone, color in [('TP', '#2ca02c'), ('TN', '#7f7f7f')]:
        sub = dft[(dft['loh_zone'] == zone) & (dft['effective_hp_reads'] > 0)]['effective_hp_reads']
        if len(sub) > 10:
            ax.hist(sub.clip(upper=300), bins=60, range=(0, 300), alpha=0.5,
                    label=f'{zone} (n={len(sub):,}, med={sub.median():.0f})',
                    color=color, density=True)
    ax.set_xlabel('Effective HP Reads')
    ax.set_ylabel('Density')
    ax.set_title(f'{mode_name} — Effective HP Reads')
    ax.legend()

plt.suptitle('Fig07: Effective HP Reads Distribution — TP LOH vs TN non-LOH', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(OUT_DIR / 'o14_fig07_effective_hp_by_zone.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved fig07")

# ============================================================================
# Summary table
# ============================================================================
print("\n" + "="*80)
print("O14 觀察摘要 — LOH 區域 Haplotag 偏向性")
print("="*80)

summary_rows = []
for mode_name, dft in [('TO', df_to), ('Paired', df_paired)]:
    dft_valid = dft[dft['hp_total'] > 0]
    for zone in ['TP', 'TN']:
        sub = dft_valid[dft_valid['loh_zone'] == zone]
        if len(sub) > 0:
            summary_rows.append({
                'mode': mode_name,
                'zone': zone,
                'n': len(sub),
                'hp_bias_median': sub['hp_bias'].median(),
                'hp_ratio_extreme_pct': ((sub['HP_Ratio'] < 0.1) | (sub['HP_Ratio'] > 0.9)).mean() * 100,
                'effective_hp_median': sub['effective_hp_reads'].median(),
                'hp0_ratio_median': sub['hp0_ratio'].median(),
            })

summary = pd.DataFrame(summary_rows)
print(summary.to_string(index=False))

# Save summary
summary.to_csv(DATA_DIR / 'o14_haplotag_bias_summary.tsv', sep='\t', index=False)
print(f"\nSaved summary to {DATA_DIR / 'o14_haplotag_bias_summary.tsv'}")

if len(rs) > 0:
    rs.to_csv(DATA_DIR / 'o14_within_region_consistency.tsv', sep='\t', index=False)
    print(f"Saved region consistency to {DATA_DIR / 'o14_within_region_consistency.tsv'}")

print("\nDone! 7 figures saved to", OUT_DIR)
