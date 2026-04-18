#!/usr/bin/env python3
"""
Generate LOH Concordance Summary Figure for weekly report Slide 10.

5-panel figure:
  A) Global 2x2 concordance matrix (Paired LOH vs TO LOH)
  B) Per-sample 2x2 concordance breakdown (stacked bar, 7 samples)
  C) Per-sample TO excess LOH ratio
  D) LOH rate by Mode × Truth (Paired/TO × TP/FP)
  E) Enrichment direction summary
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ──────────────────────────────────────────────
# Data from q1c_same_locus_concordance_corrected.tsv (TP only)
# ──────────────────────────────────────────────
samples = {
    'HCC1395':       {'both': 13107, 'paired_only': 242,  'to_only': 3898,  'neither': 10844},
    'HCC1395\nDORADO': {'both': 12514, 'paired_only': 232,  'to_only': 4918,  'neither': 10697},
    'COLO829':       {'both': 6417,  'paired_only': 207,  'to_only': 5205,  'neither': 20980},
    'H1437':         {'both': 12850, 'paired_only': 112,  'to_only': 5802,  'neither': 26079},
    'H2009':         {'both': 34700, 'paired_only': 903,  'to_only': 16334, 'neither': 72883},
    'HCC1937':       {'both': 6505,  'paired_only': 64,   'to_only': 1146,  'neither': 4111},
    'HCC1954':       {'both': 1645,  'paired_only': 96,   'to_only': 2421,  'neither': 12180},
}

sample_names = list(samples.keys())

# Aggregate
both_total = sum(s['both'] for s in samples.values())
paired_only_total = sum(s['paired_only'] for s in samples.values())
to_only_total = sum(s['to_only'] for s in samples.values())
neither_total = sum(s['neither'] for s in samples.values())
grand_total = both_total + paired_only_total + to_only_total + neither_total

# HP balance data
hp_data = {
    'Paired TP': {'extreme_loh': 0.3230, 'moderate': 0.1761, 'balanced': 0.5009, 'n': 295397},
    'Paired FP': {'extreme_loh': 0.3646, 'moderate': 0.2423, 'balanced': 0.3931, 'n': 3294},
    'TO TP':     {'extreme_loh': 0.4460, 'moderate': 0.1574, 'balanced': 0.3967, 'n': 290578},
    'TO FP':     {'extreme_loh': 0.3589, 'moderate': 0.2260, 'balanced': 0.4151, 'n': 128044},
}

# ──────────────────────────────────────────────
# Color scheme
# ──────────────────────────────────────────────
C_BOTH_LOH = '#3498db'       # Both agree LOH — blue
C_NEITHER = '#2ecc71'        # Both agree nonLOH — green
C_TO_ONLY = '#e74c3c'        # TO-only LOH — red (problematic)
C_PAIRED_ONLY = '#f39c12'    # Paired-only LOH — orange

# ──────────────────────────────────────────────
# Create figure: 3 rows, 2 columns
# ──────────────────────────────────────────────
fig = plt.figure(figsize=(18, 20))
gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.3,
                       height_ratios=[1.0, 1.2, 0.8])

fig.suptitle('LOH Concordance Analysis: Paired vs TO\n288,609 Matched Loci (TP Only, Post HP-Fix)',
             fontsize=17, fontweight='bold', y=0.99)

# ═══════════════════════════════════════════════
# Panel A: Global 2x2 Matrix
# ═══════════════════════════════════════════════
ax = fig.add_subplot(gs[0, 0])
ax.set_title('A) 全域 Concordance Matrix', fontsize=14, fontweight='bold', pad=12)

matrix = np.array([
    [neither_total, to_only_total],     # Paired=nonLOH
    [paired_only_total, both_total],     # Paired=LOH
])
pct = matrix / grand_total * 100
colors_mat = [[C_NEITHER, C_TO_ONLY], [C_PAIRED_ONLY, C_BOTH_LOH]]
descs = [
    ['共識：非LOH', 'TO 單獨判 LOH'],
    ['Paired 單獨判 LOH', '共識：LOH'],
]

for i in range(2):
    for j in range(2):
        rect = plt.Rectangle((j, 1-i), 1, 1, facecolor=colors_mat[i][j], alpha=0.75,
                               edgecolor='white', linewidth=4)
        ax.add_patch(rect)
        text_color = 'white' if colors_mat[i][j] in [C_TO_ONLY, C_BOTH_LOH] else 'black'
        ax.text(j+0.5, 1.55-i, f'{matrix[i,j]:,}', ha='center', va='center',
                fontsize=16, fontweight='bold', color=text_color)
        ax.text(j+0.5, 1.38-i, f'({pct[i,j]:.1f}%)', ha='center', va='center',
                fontsize=12, color=text_color, alpha=0.9)
        ax.text(j+0.5, 1.2-i, descs[i][j], ha='center', va='center',
                fontsize=9, color=text_color, alpha=0.8)

ax.set_xlim(0, 2); ax.set_ylim(0, 2)
ax.set_xticks([0.5, 1.5]); ax.set_xticklabels(['nonLOH', 'LOH'], fontsize=12, fontweight='bold')
ax.set_yticks([0.5, 1.5]); ax.set_yticklabels(['LOH', 'nonLOH'], fontsize=12, fontweight='bold')
ax.set_xlabel('← TO LOH Status →', fontsize=12, fontweight='bold')
ax.set_ylabel('← Paired LOH Status →', fontsize=12, fontweight='bold')
ax.set_aspect('equal')

concordance = (both_total + neither_total) / grand_total * 100
ratio_global = to_only_total / paired_only_total
ax.text(1.0, -0.12,
        f'一致率: {concordance:.1f}%  |  不一致中 TO-only 佔 {to_only_total/(to_only_total+paired_only_total)*100:.1f}%\n'
        f'TO excess ratio: {ratio_global:.0f}× (TO 系統性 over-call LOH)',
        ha='center', va='center', fontsize=10, transform=ax.transData,
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#fff3e0', edgecolor='#e65100'))

# ═══════════════════════════════════════════════
# Panel B: Per-Sample Stacked Bar (percentages)
# ═══════════════════════════════════════════════
ax = fig.add_subplot(gs[0, 1])
ax.set_title('B) Per-Sample LOH Concordance Breakdown', fontsize=14, fontweight='bold', pad=12)

x = np.arange(len(sample_names))
width = 0.65

pct_neither = []
pct_both = []
pct_to_only = []
pct_paired_only = []

for name in sample_names:
    d = samples[name]
    total = d['both'] + d['paired_only'] + d['to_only'] + d['neither']
    pct_neither.append(d['neither'] / total * 100)
    pct_both.append(d['both'] / total * 100)
    pct_to_only.append(d['to_only'] / total * 100)
    pct_paired_only.append(d['paired_only'] / total * 100)

b1 = ax.bar(x, pct_neither, width, label='共識：非LOH', color=C_NEITHER, alpha=0.85)
bottom1 = pct_neither
b2 = ax.bar(x, pct_both, width, bottom=bottom1, label='共識：LOH', color=C_BOTH_LOH, alpha=0.85)
bottom2 = [a+b for a,b in zip(bottom1, pct_both)]
b3 = ax.bar(x, pct_to_only, width, bottom=bottom2, label='TO-only LOH', color=C_TO_ONLY, alpha=0.85)
bottom3 = [a+b for a,b in zip(bottom2, pct_to_only)]
b4 = ax.bar(x, pct_paired_only, width, bottom=bottom3, label='Paired-only LOH', color=C_PAIRED_ONLY, alpha=0.85)

# Add percentage labels for key segments
for i in range(len(sample_names)):
    # Neither (center)
    if pct_neither[i] > 10:
        ax.text(i, pct_neither[i]/2, f'{pct_neither[i]:.0f}%', ha='center', va='center',
                fontsize=8, fontweight='bold', color='black')
    # Both LOH (center)
    mid_both = bottom1[i] + pct_both[i]/2
    if pct_both[i] > 8:
        ax.text(i, mid_both, f'{pct_both[i]:.0f}%', ha='center', va='center',
                fontsize=8, fontweight='bold', color='white')
    # TO-only
    mid_to = bottom2[i] + pct_to_only[i]/2
    if pct_to_only[i] > 5:
        ax.text(i, mid_to, f'{pct_to_only[i]:.0f}%', ha='center', va='center',
                fontsize=8, fontweight='bold', color='white')

ax.set_xticks(x)
ax.set_xticklabels(sample_names, fontsize=9, fontweight='bold')
ax.set_ylabel('Percentage (%)', fontsize=11)
ax.set_ylim(0, 105)
ax.legend(loc='upper right', fontsize=9, ncol=2)

# ═══════════════════════════════════════════════
# Panel C: Per-Sample TO Excess Ratio
# ═══════════════════════════════════════════════
ax = fig.add_subplot(gs[1, 0])
ax.set_title('C) Per-Sample TO/Paired LOH Excess Ratio', fontsize=14, fontweight='bold', pad=12)

ratios = []
to_counts = []
paired_counts = []
for name in sample_names:
    d = samples[name]
    r = d['to_only'] / d['paired_only'] if d['paired_only'] > 0 else 0
    ratios.append(r)
    to_counts.append(d['to_only'])
    paired_counts.append(d['paired_only'])

bar_colors = [C_TO_ONLY for _ in ratios]
bars = ax.barh(range(len(sample_names)), ratios, color=bar_colors, edgecolor='white',
               height=0.55, alpha=0.85)

for i, (r, to_c, pa_c) in enumerate(zip(ratios, to_counts, paired_counts)):
    ax.text(r + 0.8, i, f'{r:.0f}×  ({to_c:,} vs {pa_c:,})',
            va='center', fontsize=10, fontweight='bold')

ax.set_yticks(range(len(sample_names)))
ax.set_yticklabels(sample_names, fontsize=10, fontweight='bold')
ax.set_xlabel('TO-only LOH ÷ Paired-only LOH', fontsize=11)
ax.axvline(x=1, color='gray', linestyle='--', alpha=0.5)
ax.set_xlim(0, max(ratios) * 1.5)
ax.invert_yaxis()

ax.text(0.95, 0.95, f'全部樣本 TO 過判 LOH\n是 Paired 的 {min(ratios):.0f}-{max(ratios):.0f} 倍',
        transform=ax.transAxes, ha='right', va='top', fontsize=11,
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#ffe0e0', edgecolor=C_TO_ONLY, linewidth=2))

# ═══════════════════════════════════════════════
# Panel D: LOH Rate by Mode × Truth (stacked)
# ═══════════════════════════════════════════════
ax = fig.add_subplot(gs[1, 1])
ax.set_title('D) Extreme LOH Rate: Mode × Truth Label', fontsize=14, fontweight='bold', pad=12)

categories = ['Paired\nTP', 'Paired\nFP', 'TO\nTP', 'TO\nFP']
keys = ['Paired TP', 'Paired FP', 'TO TP', 'TO FP']
extreme_rates = [hp_data[k]['extreme_loh'] * 100 for k in keys]
moderate_rates = [hp_data[k]['moderate'] * 100 for k in keys]
balanced_rates = [hp_data[k]['balanced'] * 100 for k in keys]
ns = [hp_data[k]['n'] for k in keys]

x = np.arange(len(categories))
width = 0.55

b1 = ax.bar(x, extreme_rates, width, label='Extreme LOH\n(HP_Ratio <0.1 or >0.9)',
            color=C_TO_ONLY, alpha=0.85)
b2 = ax.bar(x, moderate_rates, width, bottom=extreme_rates,
            label='Moderate Imbalance', color=C_PAIRED_ONLY, alpha=0.85)
bottom2 = [e + m for e, m in zip(extreme_rates, moderate_rates)]
b3 = ax.bar(x, balanced_rates, width, bottom=bottom2,
            label='Balanced (0.3-0.7)', color=C_NEITHER, alpha=0.85)

for i, (v, n) in enumerate(zip(extreme_rates, ns)):
    ax.text(i, v/2, f'{v:.1f}%', ha='center', va='center',
            fontsize=13, fontweight='bold', color='white')
    ax.text(i, 102, f'n={n:,}', ha='center', va='bottom', fontsize=8, color='gray')

ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=11, fontweight='bold')
ax.set_ylabel('Percentage (%)', fontsize=11)
ax.set_ylim(0, 110)
ax.legend(loc='upper left', fontsize=9)

# Arrows showing direction flip
ax.annotate('', xy=(0, extreme_rates[0]+3), xytext=(1, extreme_rates[1]+3),
            arrowprops=dict(arrowstyle='<->', color='#3498db', lw=2.5))
ax.text(0.5, max(extreme_rates[0], extreme_rates[1])+7,
        f'Paired: FP > TP\n(+{extreme_rates[1]-extreme_rates[0]:.1f}pp)',
        ha='center', fontsize=10, fontweight='bold', color='#3498db')

ax.annotate('', xy=(2, extreme_rates[2]+3), xytext=(3, extreme_rates[3]+3),
            arrowprops=dict(arrowstyle='<->', color=C_TO_ONLY, lw=2.5))
ax.text(2.5, max(extreme_rates[2], extreme_rates[3])+7,
        f'TO: TP > FP ← 翻轉!\n(+{extreme_rates[2]-extreme_rates[3]:.1f}pp)',
        ha='center', fontsize=10, fontweight='bold', color=C_TO_ONLY)

# ═══════════════════════════════════════════════
# Panel E: Enrichment Direction (bottom span)
# ═══════════════════════════════════════════════
ax = fig.add_subplot(gs[2, :])
ax.set_title('E) LOH Enrichment Direction — Paired vs TO 完全翻轉', fontsize=14, fontweight='bold', pad=12)

# Per-sample enrichments
per_sample_paired = [1.02, 1.155, 1.260, 1.505, 1.795, 2.685, 3.185]
per_sample_to = [0.956, 0.940, 0.926, 0.909, 0.895, 0.866, 0.852]
per_sample_names = ['HCC1954', 'HCC1937', 'H1437', 'H2009', 'HCC1395', 'HCC1395\nDORADO', 'COLO829']

y = np.arange(len(per_sample_names))
height = 0.35

# Paired bars (right side of 1.0)
ax.barh(y + height/2, per_sample_paired, height, label='Paired (FP-enriched)',
        color=C_TO_ONLY, alpha=0.75, edgecolor='white')
# TO bars (left side of 1.0)
ax.barh(y - height/2, per_sample_to, height, label='TO (TP-enriched)',
        color=C_BOTH_LOH, alpha=0.75, edgecolor='white')

# Labels
for i, (p, t) in enumerate(zip(per_sample_paired, per_sample_to)):
    ax.text(p + 0.02, i + height/2, f'{p:.2f}×', va='center', fontsize=9, fontweight='bold', color=C_TO_ONLY)
    ax.text(t - 0.02, i - height/2, f'{t:.3f}×', va='center', ha='right', fontsize=9, fontweight='bold', color=C_BOTH_LOH)

ax.axvline(x=1.0, color='black', linestyle='-', linewidth=2.5, alpha=0.8)
ax.text(1.0, len(per_sample_names) + 0.3, 'Enrichment = 1.0\n(無差異線)', ha='center', fontsize=10, style='italic')

ax.set_yticks(y)
ax.set_yticklabels(per_sample_names, fontsize=10, fontweight='bold')
ax.set_xlabel('LOH Enrichment (FP LOH% ÷ TP LOH%)', fontsize=11)
ax.set_xlim(0.7, 3.5)
ax.legend(loc='upper right', fontsize=11)

# Pooled values
ax.text(0.5, 0.05, f'Pooled: Paired = 1.194× | TO = 0.805×\n7/7 樣本方向一致：Paired FP-enriched, TO TP-enriched',
        transform=ax.transAxes, ha='center', va='bottom', fontsize=11, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow', edgecolor='gray', linewidth=1.5))

plt.savefig('/big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/validated/2026/04/20260401_LOH_weekly_report_draft/images/concordance_summary_5panel.png',
            dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print('Done: concordance_summary_5panel.png')
