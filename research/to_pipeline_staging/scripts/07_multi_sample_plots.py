#!/usr/bin/env python3
"""
Multi-Sample ClairS-TO Comparison Plots
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
FIGDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures")
os.makedirs(FIGDIR, exist_ok=True)

with open(os.path.join(DATADIR, "multi_sample_to_summary.json")) as f:
    summary = json.load(f)
df = pd.read_csv(os.path.join(DATADIR, "multi_sample_to_features.csv"))

samples_data = summary['samples']
sample_order = ['HCC1395', 'HCC1395_DORADO', 'COLO829', 'H1437', 'H2009', 'HCC1937', 'HCC1954']

# ═══════════════════════════════════════
# Figure M1: Cross-Sample F1 Comparison
# ═══════════════════════════════════════
print("[Fig M1] Cross-sample F1 comparison...", flush=True)

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Panel A: F1 bar chart
ax = axes[0]
f1_vals = [samples_data[s]['f1'] for s in sample_order]
prec_vals = [samples_data[s]['precision'] for s in sample_order]
rec_vals = [samples_data[s]['recall'] for s in sample_order]
x = np.arange(len(sample_order))

ax.bar(x - 0.25, prec_vals, 0.25, color='#3498db', label='Precision', alpha=0.8)
ax.bar(x, rec_vals, 0.25, color='#e67e22', label='Recall', alpha=0.8)
ax.bar(x + 0.25, f1_vals, 0.25, color='#2ecc71', label='F1', alpha=0.8, edgecolor='black', linewidth=0.5)

for i, f1 in enumerate(f1_vals):
    ax.text(i + 0.25, f1 + 0.01, f'{f1:.3f}', ha='center', fontsize=7, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([s.replace('_', '\n') for s in sample_order], fontsize=7)
ax.set_ylabel('Score')
ax.set_ylim(0, 1.05)
ax.legend(fontsize=7, loc='upper left')
ax.set_title('A. ClairS-TO F1 per Sample')

# Panel B: TP/FP counts
ax = axes[1]
tp_vals = [samples_data[s]['tp'] for s in sample_order]
fp_vals = [samples_data[s]['fp'] for s in sample_order]
ax.bar(x - 0.2, tp_vals, 0.4, color='#2ecc71', label='TP', alpha=0.8)
ax.bar(x + 0.2, fp_vals, 0.4, color='#e74c3c', label='FP', alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels([s.replace('_', '\n') for s in sample_order], fontsize=7)
ax.set_ylabel('Count')
ax.legend(fontsize=8)
ax.set_title('B. TP/FP Counts per Sample')

# Panel C: TP rate (Precision)
ax = axes[2]
tp_rates = [samples_data[s]['tp'] / (samples_data[s]['tp'] + samples_data[s]['fp']) for s in sample_order]
colors_pr = ['#2ecc71' if r > 0.6 else '#e67e22' if r > 0.4 else '#e74c3c' for r in tp_rates]
bars = ax.bar(x, tp_rates, color=colors_pr, edgecolor='black', linewidth=0.5)
for i, (bar, rate) in enumerate(zip(bars, tp_rates)):
    ax.text(bar.get_x()+bar.get_width()/2, rate+0.01, f'{rate:.3f}', ha='center', fontsize=8, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([s.replace('_', '\n') for s in sample_order], fontsize=7)
ax.set_ylabel('TP Rate (Precision)')
ax.axhline(0.5, color='red', linestyle='--', alpha=0.3, label='50% (random)')
ax.set_ylim(0, 1.0)
ax.legend(fontsize=8)
ax.set_title('C. Precision (TP Rate) per Sample')

plt.suptitle('ClairS-TO PASS Benchmark: 7-Sample Comparison', fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "M01_multi_sample_f1_comparison.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: M01_multi_sample_f1_comparison.png")

# ═══════════════════════════════════════
# Figure M2: AF Distribution per Sample
# ═══════════════════════════════════════
print("[Fig M2] AF distributions...", flush=True)

fig, axes = plt.subplots(2, 4, figsize=(20, 8))
axes_flat = axes.flatten()

for idx, sample in enumerate(sample_order):
    ax = axes_flat[idx]
    sub = df[df['sample'] == sample]
    tp_sub = sub[sub['is_tp'] == 1]
    fp_sub = sub[sub['is_tp'] == 0]

    ax.hist(tp_sub['af'], bins=40, alpha=0.5, color='#2ecc71', density=True, label=f'TP (μ={tp_sub["af"].mean():.3f})')
    ax.hist(fp_sub['af'], bins=40, alpha=0.5, color='#e74c3c', density=True, label=f'FP (μ={fp_sub["af"].mean():.3f})')
    ax.set_xlabel('AF')
    f1 = samples_data[sample]['f1']
    ax.set_title(f'{sample}\nF1={f1:.3f}', fontsize=9)
    ax.legend(fontsize=6)
    ax.set_xlim(0, 1)

# Hide last subplot
axes_flat[7].set_visible(False)

plt.suptitle('AF Distribution: TP vs FP per Sample (ClairS-TO)', fontsize=13, y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "M02_multi_sample_af_distribution.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: M02_multi_sample_af_distribution.png")

# ═══════════════════════════════════════
# Figure M3: Feature AUC Heatmap
# ═══════════════════════════════════════
print("[Fig M3] Feature AUC heatmap...", flush=True)

per_sample_auc = summary['per_sample_auc']
features = ['af', 'qual', 'gq', 'dp', 'sb', 'haplotype_flag']
feature_labels = ['AF\n(inv)', 'QUAL', 'GQ', 'DP', 'SB', 'H-flag']

auc_matrix = np.zeros((len(sample_order), len(features)))
for i, sample in enumerate(sample_order):
    for j, feat in enumerate(features):
        auc = per_sample_auc[sample].get(feat, 0.5)
        # Use effective AUC (max(auc, 1-auc))
        auc_matrix[i, j] = max(auc, 1 - auc)

fig, ax = plt.subplots(figsize=(10, 6))
im = ax.imshow(auc_matrix, cmap='YlOrRd', aspect='auto', vmin=0.5, vmax=0.85)
ax.set_xticks(range(len(features)))
ax.set_xticklabels(feature_labels, fontsize=9)
ax.set_yticks(range(len(sample_order)))
ax.set_yticklabels(sample_order, fontsize=9)

# Add text values
for i in range(len(sample_order)):
    for j in range(len(features)):
        val = auc_matrix[i, j]
        color = 'white' if val > 0.7 else 'black'
        ax.text(j, i, f'{val:.3f}', ha='center', va='center', fontsize=8, color=color)

plt.colorbar(im, label='Effective AUC', shrink=0.8)
ax.set_title('VCF Feature AUC Heatmap: TP vs FP Discrimination\n(Effective AUC = max(AUC, 1-AUC), 0.5 = random)')
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "M03_multi_sample_feature_auc_heatmap.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: M03_multi_sample_feature_auc_heatmap.png")

# ═══════════════════════════════════════
# Figure M4: TO vs Paired F1 comparison
# ═══════════════════════════════════════
print("[Fig M4] TO vs Paired comparison...", flush=True)

# Paired canonical data (from earlier survey)
paired_f1 = {
    'HCC1395': 0.8551,
    'HCC1395_DORADO': 0.8662,
    'COLO829': 0.8869,
    'H1437': 0.9195,
    'H2009': 0.9710,
    'HCC1937': 0.3692,  # paired_pileup (full=0.8415)
    'HCC1954': 0.9083,
}

fig, ax = plt.subplots(figsize=(12, 6))
x = np.arange(len(sample_order))
to_f1 = [samples_data[s]['f1'] for s in sample_order]
pa_f1 = [paired_f1.get(s, 0) for s in sample_order]
delta_f1 = [pa - to for pa, to in zip(pa_f1, to_f1)]

w = 0.35
ax.bar(x - w/2, pa_f1, w, color='#3498db', alpha=0.8, label='Paired (ClairS)', edgecolor='black', linewidth=0.3)
ax.bar(x + w/2, to_f1, w, color='#e67e22', alpha=0.8, label='Tumor-Only (ClairS-TO)', edgecolor='black', linewidth=0.3)

for i in range(len(sample_order)):
    ax.annotate(f'Δ={delta_f1[i]:+.3f}', xy=(i, max(pa_f1[i], to_f1[i]) + 0.02),
                ha='center', fontsize=7, color='#c0392b', fontweight='bold')

ax.set_xticks(x)
ax.set_xticklabels([s.replace('_', '\n') for s in sample_order], fontsize=8)
ax.set_ylabel('F1 Score')
ax.set_ylim(0, 1.1)
ax.legend(fontsize=10)
ax.set_title('ClairS-TO vs ClairS Paired: F1 Comparison\n(paired_pileup for all samples)')
ax.axhline(0.5, color='gray', linestyle=':', alpha=0.3)
plt.tight_layout()
fig.savefig(os.path.join(FIGDIR, "M04_to_vs_paired_f1.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: M04_to_vs_paired_f1.png")

print("\n[DONE] Multi-sample plots generated.")
