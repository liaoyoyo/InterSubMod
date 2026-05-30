#!/usr/bin/env python3
"""Build dashboard v2 with base64-inlined figures + 5 visual cards + 23 Q list matrix.

Replaces relative <img src="..."> paths with data:image/png;base64,... data URIs
for true standalone portability (HKU can open .html anywhere without figures/ dir).

Output: docs/reports/pi_reports/2026/05/20260524_ClairS_TO_HKU_completion_dashboard_v2.standalone.html
"""

import base64
from pathlib import Path

ROOT = Path('/big7_disk/liaoyoyo2001/InterSubMod')
FIGURES_DIR = ROOT / 'research/hku_collaboration/figures'
DASHBOARD_V1 = ROOT / 'docs/reports/pi_reports/2026/05/20260524_ClairS_TO_HKU_completion_dashboard_01.standalone.html'
OUTPUT = ROOT / 'docs/reports/pi_reports/2026/05/20260524_ClairS_TO_HKU_completion_dashboard_02.standalone.html'

# Figures to inline
FIGURES = {
    'A8_per_chr_hp_summary': FIGURES_DIR / 'A8_per_chr_hp_summary_table.png',
    'A8_chr_loh_cnv': FIGURES_DIR / 'A8_chr_loh_cnv_coverage.png',
    'A7_n50_per_chr': FIGURES_DIR / 'A7_n50_per_chr_bar.png',
    'T1_HP3_3chr': FIGURES_DIR / 'T1_HP3_TP_rate_per_chr.png',
    'T2_block_len': FIGURES_DIR / 'T2_somatic_block_length_dist.png',
    # T10 24-chr will be added if exists
}

T10_FIGURE = FIGURES_DIR / 'T10_HP3_TP_rate_24chr.png'
if T10_FIGURE.exists():
    FIGURES['T10_HP3_24chr'] = T10_FIGURE


def to_data_uri(png_path: Path) -> str:
    """Convert PNG to data:image/png;base64,... URI."""
    b64 = base64.b64encode(png_path.read_bytes()).decode('ascii')
    return f"data:image/png;base64,{b64}"


def build_5_visual_cards() -> str:
    """5 strong-correspondence visual cards (T12)."""
    cards = [
        {
            'id': 'R1-Curiosity', 'icon': '🧬',
            'title': 'TSG × Methylation Phasing',
            'number': '+0.787',
            'unit': 'ΔNG (7/7 samples, p<10⁻⁶⁵)',
            'desc': 'Reviewer 1 §1.7 自行 cue InterSubMod；TP53/CDKN2A/RB1 IGV demo ready',
            'tier': '[O+F]', 'color': '#15803D'
        },
        {
            'id': 'R3-M6', 'icon': '🌳',
            'title': 'Multi-Subclone HP Linear Depth Schema',
            'number': '6-state',
            'unit': 'R/A/X matrix + linear depth tree',
            'desc': 'HP1-1-1 / HP1-1-2 / HP1-2-1 schema 補 LongPhase-S 自承 future work',
            'tier': '[I] design', 'color': '#1E3A8A'
        },
        {
            'id': 'R1-M1', 'icon': '📊',
            'title': 'Region-level Subclonal GHIR',
            'number': '93-99%',
            'unit': 'Inner same-hap NG=2 (6/6 TO samples)',
            'desc': 'TP gap +0.37 跨樣本一致；chr8 96% chr-level / 99.1% bin-level LOH × HP2-1',
            'tier': '[O]', 'color': '#A16207'
        },
        {
            'id': 'R3-M5', 'icon': '🎯',
            'title': 'HP3 是 Somatic-Evidence Bucket',
            'number': '80-94%',
            'unit': 'HP3 TP rate (chr1 88.1% / chr8 80.4% / chr19 93.7%)',
            'desc': 'HP3 fraction 0.299% 但 4-8× enrichment vs random — 非 unphased fallback',
            'tier': '[F]', 'color': '#C2410C'
        },
        {
            'id': 'R3-M1', 'icon': '🔍',
            'title': 'Rescued Variant Class Framework',
            'number': '4-Layer',
            'unit': 'F7 chr2:18,086,020 14-site Evidence',
            'desc': 'truth × LOH × HP dominance × caller blind spot 範例 stratification',
            'tier': '[O]', 'color': '#B85A3F'
        },
    ]

    html = '\n    <h3>§3.0 5 強對應 Visual Cards (T12 dashboard 補強)</h3>\n'
    html += '    <div class="visual-cards">\n'
    for c in cards:
        html += f'''      <div class="vcard" style="border-left-color: {c["color"]};">
        <div class="vcard-header">
          <span class="vcard-icon">{c["icon"]}</span>
          <span class="vcard-id">{c["id"]}</span>
          <span class="vcard-tier">{c["tier"]}</span>
        </div>
        <div class="vcard-title">{c["title"]}</div>
        <div class="vcard-number" style="color: {c["color"]};">{c["number"]}</div>
        <div class="vcard-unit">{c["unit"]}</div>
        <div class="vcard-desc">{c["desc"]}</div>
      </div>
'''
    html += '    </div>\n'
    return html


def build_q_matrix() -> str:
    """23 Q list visual matrix (T13)."""
    external = [
        ('Layer 1: linear depth schema 命名 + retraining 整合 (5)',
         '#1E3A8A', [
            '1.1 HP1-1-1 / HP1-1-2 / HP1-2-1 schema 採納？',
            '1.2 N 上限建議（IQR 下界=2, N=3 稀疏）？',
            '1.3 增獨立 channel vs 共用 channel 6？',
            '1.4 BAM aux field encoding format round-trip？',
            '1.5 HKU sub-clone ground truth dataset？',
         ]),
        ('Layer 2: ASM evidence 進 ClairS retraining (4)',
         '#15803D', [
            '2.1 加 ASM evidence 進 model？',
            '2.2 獨立 channel vs 融合 phasing_info？',
            '2.3 cell-line ASM pretrain 進去？',
            '2.4 HP1 vs HP2 imbalance ratio 內部 metric？',
         ]),
        ('Layer 3: LongPhase-S Rebuttal 對齊 (3)',
         '#C2410C', [
            '3.1 R1-Curiosity TSG IGV demo 採納？',
            '3.2 R3-M6 F7 chr2 14-site reference？',
            '3.3 R3-M5 HP3 panel — 我們幫做 8 dataset？',
         ]),
        ('Layer 4: 跨 lab data sharing + future work (3)',
         '#A16207', [
            '4.1 LongPhase-S 7+1 dataset BAM 分享？',
            '4.2 PacBio HiFi 跨平台 (COLO829 HiFi)？',
            '4.3 共同作者 / 致謝 安排？',
         ]),
    ]

    internal = [
        ('內 Layer 1: ISM 既有 evidence 待補強 (5)',
         '#7C3AED', [
            'A.1 HP3 cross-sample 6/7 樣本一致？',
            'A.2 T2 block length cross-validate A7？',
            'A.3 linear depth N=3 後 IQR 下界？',
            'A.4 Inner NG=2 跨 7 樣本 (COLO829 pending)',
            'A.5 HP3-stratified TP rate vs baseline?',
         ]),
        ('內 Layer 2: production readiness (4)',
         '#BE185D', [
            'A.6 cycle 5+ Path A/B/C 優先？',
            'A.7 caller_af normalization 解 LOSO?',
            'A.8 methyl_filter framing 改為 complementary？',
            'A.9 ISM C++ perf optimization 需要？',
         ]),
    ]

    html = '\n    <h3>§6.0 23 Q list Visual Matrix (T13 dashboard 補強)</h3>\n'
    html += '    <div class="q-matrix-section">\n'
    html += '      <h4>對 HKU 外部 collaboration ask (14 questions)</h4>\n'
    html += '      <div class="q-matrix">\n'
    for title, color, qs in external:
        html += f'        <div class="q-layer" style="border-top: 4px solid {color};">\n'
        html += f'          <div class="q-layer-title">{title}</div>\n'
        html += '          <ul>\n'
        for q in qs:
            html += f'            <li>{q}</li>\n'
        html += '          </ul>\n'
        html += '        </div>\n'
    html += '      </div>\n'

    html += '      <h4 style="margin-top: 24px;">對內 audit list (9 questions)</h4>\n'
    html += '      <div class="q-matrix">\n'
    for title, color, qs in internal:
        html += f'        <div class="q-layer" style="border-top: 4px solid {color};">\n'
        html += f'          <div class="q-layer-title">{title}</div>\n'
        html += '          <ul>\n'
        for q in qs:
            html += f'            <li>{q}</li>\n'
        html += '          </ul>\n'
        html += '        </div>\n'
    html += '      </div>\n'
    html += '    </div>\n'
    return html


def build_figures_section() -> str:
    """Inline base64 figures section."""
    has_t10 = 'T10_HP3_24chr' in FIGURES

    html = '\n    <h2><span class="h-num">§3.5</span>視覺化 evidence（T14 base64 inline portable）</h2>\n'

    if has_t10:
        html += f'''
    <h3>§3.5.1 T10 24-chr HP3 fraction + TP rate (R3 M5 完整證據)</h3>
    <figure class="inline-fig">
      <img src="{to_data_uri(FIGURES["T10_HP3_24chr"])}" alt="T10 24-chr HP3 TP rate">
      <figcaption>圖 T10: HCC1395 全 24 chr (22 autosomes + chrX) HP3 fraction（左軸對數）+ HP3 TP rate（右軸線性）；outlier chr16/chrX 紅色背景；chr8 LOH hotspot 淡藍框</figcaption>
    </figure>
'''
    else:
        html += '''
    <div class="callout-warning"><strong>T10 24-chr figure pending</strong>: BAM scan 仍在跑（~30-50 min）。完成後 dashboard rebuild 即可。當前 T1 3-chr figure 仍可顯示主要訊號。</div>
'''

    html += f'''
    <h3>§3.5.2 T1 HP3 TP rate (3-chr precursor of T10)</h3>
    <figure class="inline-fig">
      <img src="{to_data_uri(FIGURES["T1_HP3_3chr"])}" alt="T1 HP3 TP rate 3-chr">
      <figcaption>圖 T1: chr1/chr8/chr19 HP3 fraction + HP3 TP rate 雙條形圖（5/24 first run）</figcaption>
    </figure>

    <h3>§3.5.3 T2 Somatic block length distribution</h3>
    <figure class="inline-fig">
      <img src="{to_data_uri(FIGURES["T2_block_len"])}" alt="T2 block length">
      <figcaption>圖 T2: chr1/chr8/chr19 PS block 長度分布（左 log-x histogram + median/N50/17.4 kb 參考線；右 per-chr boxplot log-y）。N50 1.04 Mb >> ONT read 10-20 kb</figcaption>
    </figure>

    <h3>§3.5.4 A7 全 24 chr PS block N50</h3>
    <figure class="inline-fig">
      <img src="{to_data_uri(FIGURES["A7_n50_per_chr"])}" alt="A7 N50 per chr">
      <figcaption>圖 A7: 全 24 chr PS block N50 bar chart。chr22 max 1.625 Mb / chr16 min 227 kb / chrX 101 kb — phasing-weak chr 對應 short N50</figcaption>
    </figure>

    <h3>§3.5.5 A8 全 24 chr HP composition stacked bar</h3>
    <figure class="inline-fig">
      <img src="{to_data_uri(FIGURES["A8_per_chr_hp_summary"])}" alt="A8 per-chr HP summary">
      <figcaption>圖 A8: 23 chr × HP1/HP2/HP1-1/HP2-1/HP3/no_HP normalized stack。chr6 HP1-1+HP2-1+HP3 ~9% outlier / chr16 no_HP 52.9% / chrX no_HP 61.1% HP3 12% — 與 A7 N50 短的 chr 對應</figcaption>
    </figure>

    <h3>§3.5.6 A8 LOH/CNV coverage per-chr</h3>
    <figure class="inline-fig">
      <img src="{to_data_uri(FIGURES["A8_chr_loh_cnv"])}" alt="A8 LOH/CNV coverage">
      <figcaption>圖 A8 (LOH/CNV): chr8 96% / chr17 93% / chr11 85% Tier 1 LOH；chr16 chr20 chrX = 0% LOH baseline reference</figcaption>
    </figure>
'''
    return html


def main():
    # Read dashboard v1
    v1_html = DASHBOARD_V1.read_text()

    # Inject extra CSS for visual cards + q matrix + inline-fig
    extra_css = '''
    /* Visual cards (T12) */
    .visual-cards {
      display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 16px; margin: 16px 0;
    }
    .vcard {
      background: var(--c-surface); border: 1px solid var(--c-border);
      border-left-width: 5px; border-radius: 6px; padding: 16px;
    }
    .vcard-header {
      display: flex; align-items: center; gap: 8px; margin-bottom: 8px;
      font-size: 12px; color: var(--c-text-soft);
    }
    .vcard-icon { font-size: 24px; }
    .vcard-id { font-family: var(--ff-mono); font-weight: 700; background: var(--c-code-bg); padding: 2px 8px; border-radius: 10px; }
    .vcard-tier { font-family: var(--ff-mono); font-weight: 700; background: #EFF6FF; color: #1E40AF; padding: 2px 8px; border-radius: 10px; }
    .vcard-title { font-weight: 600; font-size: 15px; margin: 4px 0 8px; line-height: 1.35; }
    .vcard-number { font-family: var(--ff-mono); font-weight: 700; font-size: 28px; line-height: 1.1; }
    .vcard-unit { font-size: 13px; color: var(--c-text-soft); margin: 4px 0 8px; }
    .vcard-desc { font-size: 13px; line-height: 1.55; color: var(--c-text); }

    /* Q matrix (T13) */
    .q-matrix-section { margin: 16px 0; }
    .q-matrix-section h4 {
      font-size: 14px; text-transform: uppercase; letter-spacing: 0.5px;
      color: var(--c-text-soft); margin: 16px 0 8px; font-weight: 700;
    }
    .q-matrix {
      display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 12px;
    }
    .q-layer {
      background: var(--c-surface); border: 1px solid var(--c-border);
      border-radius: 6px; padding: 12px;
    }
    .q-layer-title {
      font-weight: 600; font-size: 13px; color: var(--c-text);
      margin-bottom: 8px; line-height: 1.4;
    }
    .q-layer ul { margin: 0; padding-left: 18px; }
    .q-layer li { font-size: 13px; line-height: 1.6; margin: 4px 0; color: var(--c-text); }

    /* Inline figures (T14 base64) */
    figure.inline-fig {
      margin: 16px 0; text-align: center;
      background: var(--c-bg); border: 1px solid var(--c-border);
      border-radius: 6px; padding: 12px;
    }
    figure.inline-fig img {
      max-width: 100%; height: auto; display: block; margin: 0 auto;
      border-radius: 3px;
    }
    figure.inline-fig figcaption {
      font-size: 13px; color: var(--c-text-soft);
      margin-top: 8px; line-height: 1.5; text-align: left;
    }
'''
    # Insert extra CSS before </style>
    v2_html = v1_html.replace('</style>', extra_css + '\n  </style>')

    # Build new sections
    cards_html = build_5_visual_cards()
    q_matrix_html = build_q_matrix()
    figures_html = build_figures_section()

    # Insert 5 visual cards before §3.1 強對應 5 條
    v2_html = v2_html.replace(
        '<h3>§3.1 強對應 5 條（✅ 可主動引用作 LongPhase-S Supp）</h3>',
        cards_html + '\n    <h3>§3.1 強對應 5 條（✅ 可主動引用作 LongPhase-S Supp）— 詳細表</h3>'
    )

    # Insert figures section before §4 T8 Evaluator
    v2_html = v2_html.replace(
        '<!-- T8 Evaluator Fixes -->',
        figures_html + '\n    <!-- T8 Evaluator Fixes -->'
    )

    # Insert Q matrix before §6 Outstanding
    v2_html = v2_html.replace(
        '<!-- Outstanding -->',
        q_matrix_html + '\n    <!-- Outstanding -->'
    )

    # Update version in hero
    v2_html = v2_html.replace(
        '<dt>狀態</dt><dd>D3 完成 + T8 evaluator 修正 + dashboard</dd>',
        '<dt>狀態</dt><dd>D3 v2 完成 + T10-T14 視覺化 + base64 portable</dd>'
    )
    v2_html = v2_html.replace(
        '<title>ClairS-TO HKU Handoff — 5/24 Completion Dashboard</title>',
        '<title>ClairS-TO HKU Handoff — 5/24 Completion Dashboard v2 (portable)</title>'
    )

    # Write output
    OUTPUT.write_text(v2_html)
    print(f"Written: {OUTPUT}")
    print(f"Size: {OUTPUT.stat().st_size // 1024} KB")
    print(f"Figures inlined: {list(FIGURES.keys())}")
    print(f"T10 24-chr included: {'T10_HP3_24chr' in FIGURES}")


if __name__ == '__main__':
    main()
