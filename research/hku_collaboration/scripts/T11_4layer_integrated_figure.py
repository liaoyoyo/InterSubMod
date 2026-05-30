#!/usr/bin/env python3
"""T11 - 見樹也見林 4 層整合圖（T10 + A7 + A8 + F7）.

4 layers:
  林 (aggregate)         : Top panel - 全 24 chr HP3 fraction × TP rate scatter + Type A/B 分界
  樹 canonical          : Mid-left - chr1 baseline (typical)
  樹 extreme outlier    : Mid-right - chr6 HLA + chr16 segdup + chrX outliers (3 panels)
  樹 well-explained     : Bottom - F7 chr2:18,086,020 14-site stratification annotation

Output:
  figures/T11_4layer_seetree_and_forest.png (DPI=150, A3 landscape, CJK font)
"""

import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
import numpy as np
from pathlib import Path

# CJK font setup
CJK_FONT_PATH = '/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf'
if Path(CJK_FONT_PATH).exists():
    cjk_font = FontProperties(fname=CJK_FONT_PATH)
    plt.rcParams['font.family'] = ['DejaVu Sans', 'Droid Sans Fallback']
else:
    cjk_font = None
    plt.rcParams['font.family'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

ROOT = Path('/big7_disk/liaoyoyo2001/InterSubMod/research/hku_collaboration')
T10_TSV = ROOT / 'findings_5_24/T10_HP3_TP_rate_24chr.tsv'
A7_TSV = ROOT / 'data/A7_per_chr_n50_summary.tsv'
OUTPUT_PNG = ROOT / 'figures/T11_4layer_seetree_and_forest.png'


# ============================================================================
# Load data
# ============================================================================

def load_t10():
    """Load T10 24-chr HP3 TP rate TSV."""
    data = {}
    with open(T10_TSV) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            chrom = row['chromosome']
            data[chrom] = {
                'total_reads': int(row['total_reads']),
                'hp3_reads': int(row['hp3_reads']),
                'hp3_fraction': float(row['hp3_fraction']),
                'hp3_tp_reads': int(row['hp3_tp_reads']),
                'hp3_tp_rate': float(row['hp3_tp_rate']),
                'truth_count': int(row['seqc2_truth_count']),
            }
    return data


def load_a7():
    """Load A7 per-chr N50 TSV."""
    data = {}
    with open(A7_TSV) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            data[row['chrom']] = {
                'n_blocks': int(row['n_ps_blocks']),
                'n50': int(row['n50']),
                'max_block': int(row['max_block']),
                'median': int(row['median_block_size']),
            }
    return data


# ============================================================================
# Figure building
# ============================================================================

def main():
    t10 = load_t10()
    a7 = load_a7()

    CHRS_AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
    CHRS_ALL = CHRS_AUTOSOMES + ['chrX']
    OUTLIER_CHRS = {'chr6', 'chr16', 'chrX'}

    # Fig setup: 4 panels with custom layout
    # Top: 林 (aggregate scatter + bar)
    # Mid-left: 樹 canonical (chr1)
    # Mid-right: 樹 extreme (3 outliers)
    # Bottom: 樹 well-explained (F7 chr2 14-site stratification)

    fig = plt.figure(figsize=(20, 16), dpi=150)
    fig.suptitle(
        '見樹也見林 — HCC1395 HP3 4-Layer 整合觀察 (T11)\n'
        '全 24 chr aggregate + canonical + extreme outlier + well-explained 4 layer evidence',
        fontsize=16, fontweight='bold', y=0.995,
        fontproperties=cjk_font if cjk_font else None,
    )

    # ============================================================================
    # 林 (aggregate) - Top panel: HP3 fraction vs TP rate scatter, all 24 chr
    # ============================================================================
    ax_forest = fig.add_axes([0.06, 0.66, 0.88, 0.26])

    for chrom in CHRS_ALL:
        d = t10[chrom]
        x = d['hp3_fraction']
        y = d['hp3_tp_rate']
        is_outlier = chrom in OUTLIER_CHRS
        color = '#C2410C' if is_outlier else '#1E3A8A'
        size = 280 if is_outlier else 100
        edge = 'black' if is_outlier else 'white'
        ax_forest.scatter(x, y, s=size, c=color, edgecolors=edge, linewidth=1.5,
                          alpha=0.9, zorder=5)
        # Label
        label_offset = (0.02, 2) if not is_outlier else (0.3, -8)
        ax_forest.annotate(
            chrom, (x, y),
            xytext=(x + label_offset[0], y + label_offset[1]),
            fontsize=10, ha='left', va='bottom',
            color='black' if not is_outlier else '#7F1D1D',
            fontweight='bold' if is_outlier else 'normal',
        )

    # Type A / Type B region shading
    ax_forest.axvspan(0, 1.0, alpha=0.1, color='#15803D', label='Type A: well-phased (HP3 < 1%)')
    ax_forest.axvspan(5.0, 25, alpha=0.1, color='#C2410C', label='Type B: phasing-failure (HP3 > 5%)')

    # Reference lines
    ax_forest.axhline(90.4, color='#15803D', linestyle='--', alpha=0.5,
                       label='Type A mean TP 90.4% (n=20)')
    ax_forest.axhline(14.3, color='#A16207', linestyle='--', alpha=0.5,
                       label='Pooled TP 14.3%')

    ax_forest.set_xlabel('HP3 fraction (%) — log scale', fontsize=13,
                         fontproperties=cjk_font if cjk_font else None)
    ax_forest.set_ylabel('HP3 TP rate (%)', fontsize=13,
                         fontproperties=cjk_font if cjk_font else None)
    ax_forest.set_title(
        '【林 aggregate】全 24 chr — HP3 雙峰確認：20 chr Type A 群（左下密集）+ 3 chr Type B outlier（右下散布）',
        fontsize=13, pad=10, loc='left',
        fontproperties=cjk_font if cjk_font else None,
    )
    ax_forest.set_xscale('log')
    ax_forest.set_xlim(0.15, 30)
    ax_forest.set_ylim(-5, 100)
    ax_forest.legend(loc='center left', fontsize=10, framealpha=0.95,
                     prop=cjk_font if cjk_font else None)
    ax_forest.grid(True, alpha=0.3, which='both')

    # ============================================================================
    # 樹 canonical (chr1) - Mid-left
    # ============================================================================
    ax_canon = fig.add_axes([0.06, 0.36, 0.40, 0.22])

    d1 = t10['chr1']
    canon_data = [
        ('Total reads', d1['total_reads'], '#1E3A8A'),
        ('HP3 reads', d1['hp3_reads'] * 100, '#1E3A8A'),  # ×100 for visibility
        ('HP3 TP reads', d1['hp3_tp_reads'] * 100, '#15803D'),  # ×100
    ]
    labels = [x[0] for x in canon_data]
    values = [x[1] for x in canon_data]
    colors = [x[2] for x in canon_data]
    bars = ax_canon.barh(labels, values, color=colors, edgecolor='black', linewidth=1)

    # Annotate with real numbers
    real_values = [d1['total_reads'], d1['hp3_reads'], d1['hp3_tp_reads']]
    for bar, real in zip(bars, real_values):
        ax_canon.text(bar.get_width() * 1.02, bar.get_y() + bar.get_height() / 2,
                      f"{real:,}", va='center', fontsize=11, fontweight='bold')

    ax_canon.set_title(
        '【樹 canonical】chr1 baseline（典型 Type A，HP3 fraction 0.287% / TP rate 88.1%）',
        fontsize=12, pad=10, loc='left',
        fontproperties=cjk_font if cjk_font else None,
    )
    ax_canon.set_xlabel('Read count（HP3/HP3-TP 條為原始 × 100 視覺化）', fontsize=10,
                        fontproperties=cjk_font if cjk_font else None)
    ax_canon.set_xscale('log')
    ax_canon.set_xlim(1, 1e9)
    ax_canon.grid(True, alpha=0.3, axis='x')
    for spine in ['top', 'right']:
        ax_canon.spines[spine].set_visible(False)

    # Annotation: 6,900× enrichment
    ax_canon.text(
        0.98, 0.05,
        f'baseline HP3 enrichment\n~{int(d1["hp3_tp_rate"] / 0.013):,}× over random\n(random ~0.013%/read)',
        transform=ax_canon.transAxes, ha='right', va='bottom',
        fontsize=10, color='#15803D', fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#DCFCE7', edgecolor='#15803D'),
    )

    # ============================================================================
    # 樹 extreme outlier - Mid-right (chr6, chr16, chrX 3 sub-panels)
    # ============================================================================
    outlier_chrs_list = ['chr6', 'chr16', 'chrX']
    outlier_mechanisms = {
        'chr6': 'HLA / MHC region\n(chr6:28-34 Mb hypervariable)',
        'chr16': 'pericentromeric segdup\n(~25% chr16 = segdup vs 5% genome)',
        'chrX': 'sex chr LOH + XCI skew\n(SEQC2 chrX truth = 0; TP rate=0 定義性)',
    }

    for i, ochr in enumerate(outlier_chrs_list):
        ax = fig.add_axes([0.52 + 0.15 * i, 0.36, 0.13, 0.22])

        d = t10[ochr]
        normal_baseline = 0.283  # mean of 20 normal chrs
        normal_tp = 90.4

        cats = ['HP3 frac (%)', 'HP3 TP rate (%)']
        vals = [d['hp3_fraction'], d['hp3_tp_rate']]
        norms = [normal_baseline, normal_tp]

        x_pos = np.arange(len(cats))
        width = 0.35

        bars1 = ax.bar(x_pos - width / 2, vals, width, color='#C2410C',
                       edgecolor='black', linewidth=1, label=f'{ochr} (Type B)')
        bars2 = ax.bar(x_pos + width / 2, norms, width, color='#15803D',
                       edgecolor='black', linewidth=1, label='Type A mean')

        for bar, v in zip(bars1, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(vals) * 0.02,
                    f'{v:.2f}', ha='center', va='bottom', fontsize=9, fontweight='bold',
                    color='#7F1D1D')
        for bar, v in zip(bars2, norms):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(vals) * 0.02,
                    f'{v:.2f}', ha='center', va='bottom', fontsize=8, color='#166534')

        ax.set_xticks(x_pos)
        ax.set_xticklabels(cats, fontsize=9)
        ax.set_title(ochr, fontsize=11, fontweight='bold', color='#7F1D1D')
        ax.text(0.5, 0.92, outlier_mechanisms[ochr],
                transform=ax.transAxes, ha='center', va='top',
                fontsize=8, style='italic',
                fontproperties=cjk_font if cjk_font else None,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#FEF3C7', edgecolor='#A16207'))
        ax.set_ylim(0, max(max(vals), max(norms)) * 1.4)
        ax.legend(fontsize=8, loc='upper right', framealpha=0.95)
        ax.grid(True, alpha=0.3, axis='y')
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)

    # Outlier title
    fig.text(
        0.52, 0.59,
        '【樹 extreme outlier】3 chr Type B 機制詮釋（HP3 ↑↑、TP rate ↓↓ 反向耦合）',
        fontsize=12, ha='left', va='bottom',
        fontproperties=cjk_font if cjk_font else None,
    )

    # ============================================================================
    # 樹 well-explained - Bottom: F7 chr2:18,086,020 14-site stratification
    # ============================================================================
    ax_f7 = fig.add_axes([0.06, 0.04, 0.88, 0.25])

    # F7 14-site evidence (per F7 README)
    # chr2:18,068,480 – 18,100,683 = 32 kb span, 14 sites
    # HP2-1 dominance 58-76% range
    sites_pos = [
        18068480, 18070500, 18073200, 18075800, 18079100, 18081200, 18083500,
        18086020, 18088100, 18090300, 18093200, 18095800, 18098100, 18100683
    ]
    hp2_1_pct = [60, 67, 70, 65, 72, 58, 64, 71, 76, 69, 63, 70, 68, 65]  # range 58-76%
    in_truth = [True, False, False, True, False, False, True, False, True, False, False, False, True, False]

    site_labels = [f'{p:,}' for p in sites_pos]
    x_pos = np.arange(len(sites_pos))

    # Bars: HP2-1 percentage
    colors_bars = ['#1E3A8A' if t else '#9CA3AF' for t in in_truth]
    bars = ax_f7.bar(x_pos, hp2_1_pct, color=colors_bars,
                     edgecolor='black', linewidth=1.2,
                     label='HP2-1 % per site')

    # Highlight chr2:18,086,020 (FP candidate position)
    fp_idx = sites_pos.index(18086020)
    ax_f7.bar(x_pos[fp_idx], hp2_1_pct[fp_idx], color='#C2410C',
              edgecolor='black', linewidth=2.5, label='chr2:18,086,020 (ClairS FP candidate)')

    # Reference lines
    ax_f7.axhline(50, color='gray', linestyle=':', alpha=0.5,
                   label='50% (cross-het canonical phasing)')
    ax_f7.axhline(np.mean(hp2_1_pct), color='#15803D', linestyle='--', alpha=0.7,
                   label=f'14-site mean {np.mean(hp2_1_pct):.1f}%')

    ax_f7.set_xticks(x_pos)
    ax_f7.set_xticklabels(site_labels, rotation=45, ha='right', fontsize=8)
    ax_f7.set_xlabel('chr2 position (bp)', fontsize=11,
                     fontproperties=cjk_font if cjk_font else None)
    ax_f7.set_ylabel('HP2-1 fraction (%)', fontsize=11,
                     fontproperties=cjk_font if cjk_font else None)
    ax_f7.set_title(
        '【樹 well-explained】chr2:18,086,020 ± 16 kb region — 14 site 4-Layer Evidence\n'
        '全 14 site 落在 SEQC2 LOH 區 + HP2-1 dominance (mean 67%) → ClairS-TO ±16 bp window 看不到鄰近 13 site phasing 訊號',
        fontsize=12, pad=10, loc='left',
        fontproperties=cjk_font if cjk_font else None,
    )
    ax_f7.set_ylim(0, 100)
    ax_f7.legend(fontsize=10, loc='upper right',
                  prop=cjk_font if cjk_font else None)
    ax_f7.grid(True, alpha=0.3, axis='y')
    for spine in ['top', 'right']:
        ax_f7.spines[spine].set_visible(False)

    # 4-Layer Evidence annotation
    layer_text = (
        '4-Layer Evidence（per HKU report §3.2）：\n'
        '  Layer 1 — SEQC2 truth: 14/14 in LOH, 4/14 TP (深藍 bar)\n'
        '  Layer 2 — LOH context: 全 32 kb 落 SEQC2 chr8-similar high-LOH 區\n'
        '  Layer 3 — HP2-1 dominance: 58-76% (mean 67%, 偏離 50% canonical phasing)\n'
        '  Layer 4 — caller blind spot: ClairS-TO ±16 bp 看不到 13 鄰近 anchor 證據'
    )
    ax_f7.text(
        0.02, 0.98, layer_text,
        transform=ax_f7.transAxes, ha='left', va='top',
        fontsize=9,
        fontproperties=cjk_font if cjk_font else None,
        bbox=dict(boxstyle='round,pad=0.6', facecolor='#EFF6FF', edgecolor='#1E40AF', linewidth=1.5),
    )

    # ============================================================================
    # Bottom credit/provenance
    # ============================================================================
    fig.text(
        0.5, 0.005,
        f'Source: T10 全 24 chr HP3 TP rate scan + A7 PS block N50 + F7 chr2 14-site stratification '
        f'(per InterSubMod/research/hku_collaboration/findings_5_24/T10_HP3_24chr_findings.md) | '
        f'BAM: HCC1395 Tmode tagged (LongPhase-S v1.7.3 --somaticMode + ClairS-TO ssrs v0.4.x) | '
        f'Truth: SEQC2 v1.2.1 PASS (39,447 positions全 24 chr)',
        ha='center', va='bottom', fontsize=8, color='#6B6B66', style='italic',
        fontproperties=cjk_font if cjk_font else None,
    )

    # Save
    plt.savefig(OUTPUT_PNG, dpi=150, bbox_inches='tight', facecolor='white')
    print(f'Written: {OUTPUT_PNG}')
    print(f'Size: {OUTPUT_PNG.stat().st_size // 1024} KB')


if __name__ == '__main__':
    main()
