#!/usr/bin/env python3
"""
Updated pipeline waterfall figure with proper LongPhase-TO stage
and ISM middle block FP analysis.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
FIGDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures")

# ═══════════════════════════════════════
# Figure 8: Corrected Pipeline Waterfall
# ═══════════════════════════════════════
print("[Fig 8] Pipeline waterfall (corrected)...", flush=True)

fig, axes = plt.subplots(1, 2, figsize=(18, 7))

# Panel A: Full pipeline waterfall
stages = [
    ('All ClairS-TO\ncalls', 30247, 3157028),
    ('After PoN\n(NonSomatic)', 29342, 24459),
    ('After Quality\nfilters', 28633, 19419),
    ('After Phasing\nfilters', 28509, 19289),
    ('After ISM\nimplicit filter', 28013, 6100),
    ('After ISM\nSuggestFilter', 27978, 6092),
]

x = np.arange(len(stages))
tp_vals = [s[1] for s in stages]
fp_vals = [s[2] for s in stages]

# Use log scale for FP (huge range)
ax = axes[0]
ax2 = ax.twinx()

bars_tp = ax.bar(x - 0.2, tp_vals, 0.35, color='#2ecc71', alpha=0.8, label='TP')
bars_fp = ax2.bar(x + 0.2, fp_vals, 0.35, color='#e74c3c', alpha=0.8, label='FP')
ax2.set_yscale('log')

ax.set_xticks(x)
ax.set_xticklabels([s[0] for s in stages], fontsize=8)
ax.set_ylabel('TP count', color='#2ecc71')
ax2.set_ylabel('FP count (log)', color='#e74c3c')
ax.set_title('A. Pipeline Waterfall: TP/FP at Each Stage')

# Add annotations for key drops
for i in range(1, len(stages)):
    dtp = tp_vals[i-1] - tp_vals[i]
    dfp = fp_vals[i-1] - fp_vals[i]
    if dfp > 100:
        ax.annotate(f'-{dfp:,} FP\n-{dtp} TP', xy=(i, tp_vals[i]),
                    fontsize=6, ha='center', va='bottom', color='#c0392b',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='#fadbd8', alpha=0.7))

# Legends
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1+lines2, labels1+labels2, loc='upper right', fontsize=8)

# Panel B: Middle block analysis — can ISM features remove the 6,100 FP?
ax = axes[1]

strategies = [
    ('ISM baseline\n(implicit filter)', 28013, 6100, '#3498db'),
    ('+ CNV\n(LOH|Gain)', 26930, 2481, '#2ecc71'),
    ('+ AF < 0.55', 18748, 4203, '#e67e22'),
    ('+ CNV +\nAF<0.55', 18014, 2261, '#9b59b6'),
    ('+ AlleleSig\nonly (ISM)', 15679, 3021, '#e74c3c'),
]

x = np.arange(len(strategies))
for i, (name, tp, fp, color) in enumerate(strategies):
    fn = 39447 - tp
    prec = tp / max(tp + fp, 1)
    rec = tp / 39447
    f1 = 2*prec*rec/(prec+rec) if (prec+rec)>0 else 0

    ax.bar(i, f1, 0.6, color=color, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax.text(i, f1 + 0.01, f'F1={f1:.3f}\nP={prec:.3f}\nR={rec:.3f}', ha='center', fontsize=7, va='bottom')

ax.set_xticks(x)
ax.set_xticklabels([s[0] for s in strategies], fontsize=8)
ax.set_ylabel('F1 Score')
ax.set_ylim(0, 0.95)
ax.set_title('B. Can We Remove the Surviving 6,100 FP?\n(ISM features vs CNV annotation)')
ax.axhline(0.7616, color='gray', linestyle='--', alpha=0.5, label='ISM baseline F1')
ax.legend(fontsize=8)

plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "08_corrected_waterfall_and_middle_block.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: 08_corrected_waterfall_and_middle_block.png")

# ═══════════════════════════════════════
# Figure 9: ISM Feature Discrimination in Middle Block
# ═══════════════════════════════════════
print("\n[Fig 9] ISM features in middle block...", flush=True)

df = pd.read_csv(os.path.join(DATADIR, "hcc1395_pass_multimodal.csv"))
ism = df[df['has_ism']==1].copy()

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
tp = ism[ism['is_tp']==1]
fp = ism[ism['is_tp']==0]

# Panel A: AlleleSig rate by TP/FP × CNV
cnv_types = ['neutral', 'gain', 'loh']
x = np.arange(len(cnv_types))
tp_sig_rates = []
fp_sig_rates = []
for ct in cnv_types:
    tp_sub = tp[tp['cnv_type']==ct]
    fp_sub = fp[fp['cnv_type']==ct]
    tp_sig_rates.append(tp_sub['ism_AlleleSig'].mean() if len(tp_sub) > 10 else 0)
    fp_sig_rates.append(fp_sub['ism_AlleleSig'].mean() if len(fp_sub) > 10 else 0)

axes[0,0].bar(x-0.15, tp_sig_rates, 0.3, color='#2ecc71', label='TP')
axes[0,0].bar(x+0.15, fp_sig_rates, 0.3, color='#e74c3c', label='FP')
axes[0,0].set_xticks(x)
axes[0,0].set_xticklabels(cnv_types)
axes[0,0].set_ylabel('AlleleSig Rate')
axes[0,0].set_title('A. AlleleSig by CNV × TP/FP')
axes[0,0].legend(fontsize=8)

# Panel B: AlleleDelta by CNV
for ct, marker in [('neutral', 'o'), ('gain', 's'), ('loh', 'D')]:
    tp_sub = tp[tp['cnv_type']==ct]
    fp_sub = fp[fp['cnv_type']==ct]
    ad_tp = tp_sub['ism_AlleleDelta'].dropna().values
    ad_fp = fp_sub['ism_AlleleDelta'].dropna().values
    if len(ad_tp) > 10:
        axes[0,1].hist(ad_tp, bins=30, alpha=0.3, color='#2ecc71', density=True, label=f'{ct} TP' if ct=='neutral' else '')
    if len(ad_fp) > 10:
        axes[0,1].hist(ad_fp, bins=30, alpha=0.3, color='#e74c3c', density=True)
axes[0,1].set_xlabel('AlleleDelta')
axes[0,1].set_title('B. AlleleDelta Distribution\n(all CNV types overlaid)')
axes[0,1].set_xlim(0, 0.3)

# Panel C: Quality_Score by TP/FP
qs_tp = tp['ism_Quality_Score'].dropna()
qs_fp = fp['ism_Quality_Score'].dropna()
axes[0,2].hist(qs_tp[qs_tp>0], bins=50, alpha=0.6, color='#2ecc71', label=f'TP (mean={qs_tp.mean():.1f})', density=True)
axes[0,2].hist(qs_fp[qs_fp>0], bins=50, alpha=0.6, color='#e74c3c', label=f'FP (mean={qs_fp.mean():.1f})', density=True)
axes[0,2].set_xlabel('Quality_Score')
axes[0,2].set_title('C. Quality_Score (FP slightly higher)')
axes[0,2].legend(fontsize=8)

# Panel D: NumReads by TP/FP
nr_tp = tp['ism_NumReads'].clip(0, 500)
nr_fp = fp['ism_NumReads'].clip(0, 500)
axes[1,0].hist(nr_tp, bins=50, alpha=0.6, color='#2ecc71', density=True, label=f'TP (mean={tp["ism_NumReads"].mean():.0f})')
axes[1,0].hist(nr_fp, bins=50, alpha=0.6, color='#e74c3c', density=True, label=f'FP (mean={fp["ism_NumReads"].mean():.0f})')
axes[1,0].set_xlabel('NumReads')
axes[1,0].set_title('D. NumReads (TP slightly more)')
axes[1,0].legend(fontsize=8)

# Panel E: CNV × TP rate (key comparison)
cnv_types_full = ['neutral', 'gain', 'loss', 'loh']
tp_rates = []
for ct in cnv_types_full:
    sub = ism[ism['cnv_type']==ct]
    tp_rates.append(sub['is_tp'].mean())

colors = ['#bdc3c7', '#2ecc71', '#e74c3c', '#3498db']
bars = axes[1,1].bar(cnv_types_full, tp_rates, color=colors, edgecolor='black', linewidth=0.5)
axes[1,1].axhline(ism['is_tp'].mean(), color='red', linestyle='--', alpha=0.5,
                   label=f'Overall={ism["is_tp"].mean():.3f}')
axes[1,1].set_ylabel('TP Rate')
axes[1,1].set_title('E. CNV Status (ISM subset)')
axes[1,1].legend(fontsize=8)
for bar, rate in zip(bars, tp_rates):
    axes[1,1].text(bar.get_x()+bar.get_width()/2, rate+0.01, f'{rate:.3f}', ha='center', fontsize=9, fontweight='bold')

# Panel F: Summary AUC bar chart for middle block
features_auc = {
    'in_gain\n(CNV)': 0.651,
    'in_loh\n(CNV)': 0.608,
    'af\n(VCF)': 0.565,
    'NumReads\n(ISM)': 0.559,
    'AlleleDelta\n(ISM)': 0.520,
    'AlleleSig\n(ISM)': 0.475,  # inverted
    'QualityScore\n(ISM)': 0.444,  # inverted
}

x_f = np.arange(len(features_auc))
vals = list(features_auc.values())
colors_f = ['#e67e22' if v >= 0.6 else '#3498db' if v > 0.52 else '#e74c3c' for v in vals]
axes[1,2].barh(x_f, vals, color=colors_f, alpha=0.8)
axes[1,2].set_yticks(x_f)
axes[1,2].set_yticklabels(list(features_auc.keys()), fontsize=8)
axes[1,2].axvline(0.5, color='gray', linestyle='--', alpha=0.5)
axes[1,2].set_xlabel('AUC (within ISM subset)')
axes[1,2].set_title('F. Feature AUC for Middle Block\n(ISM features ≤ 0.56)')
for i, v in enumerate(vals):
    axes[1,2].text(v+0.005, i, f'{v:.3f}', va='center', fontsize=8)

plt.suptitle('Middle Block Analysis: Can ISM Features Remove the Surviving 6,100 FP?', fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "09_middle_block_ism_features.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: 09_middle_block_ism_features.png")

print("\n[DONE]")
