<!--
建立時間: 2026-04-02 22:30
目標: 擴展 methyl_mean LOH 分析到 TO mode read-level 數據
處理範圍: sample400 shard (HCC1395 + HCC1395_DORADO) × TO + Paired
關聯檔案:
  - research/loh_investigation/scripts/analyze_to_read_methyl_in_loh.py
  - research/loh_investigation/figures/to_read_methyl_fig*.png
  - research/loh_investigation/data/to_read_methyl_stats.tsv
-->

# TO Mode Read-Level Methylation Inside LOH

**Date**: 2026-04-02
**Data**: sample400 shard export (HCC1395 × 2 platforms × 2 modes)
**LOH zones**: per-sample LongPhase-TO LOH.bed (binary inside/outside)

## Data Summary

- **TO inside**: 596 TP reads (4 regions) + 3133 FP reads (21 regions)
- **TO outside**: 15970 TP reads (86 regions) + 11618 FP reads (77 regions)
- **Paired inside**: 0 TP reads (0 regions) + 726 FP reads (7 regions)
- **Paired outside**: 16638 TP reads (90 regions) + 10114 FP reads (86 regions)

## Key Results: methyl_mean Inside LOH

| Mode | Zone | Feature | n_TP | n_FP | TP median | FP median | AUC | Cohen's d |
|------|------|---------|------|------|-----------|-----------|-----|-----------|
| TO | inside | methyl_mean | 596 | 3133 | 0.7706 | 0.8447 | 0.4404 | -0.1201 |
| TO | inside | methyl_high_fraction | 596 | 3133 | 0.7064 | 0.8378 | 0.4185 | -0.2087 |
| TO | inside | methyl_low_fraction | 596 | 3133 | 0.1883 | 0.1455 | 0.5223 | 0.0463 |
| TO | outside | methyl_mean | 15970 | 11618 | 0.5718 | 0.6975 | 0.3735 | -0.4521 |
| TO | outside | methyl_high_fraction | 15970 | 11618 | 0.5473 | 0.6578 | 0.3826 | -0.4205 |
| TO | outside | methyl_low_fraction | 15970 | 11618 | 0.4146 | 0.2634 | 0.6360 | 0.4837 |
| TO | inside | region_methyl_mean | 4 | 21 | 0.6922 | 0.8607 | 0.4167 | -0.0983 |
| TO | outside | region_methyl_mean | 86 | 77 | 0.6035 | 0.7291 | 0.3819 | -0.4284 |
| Paired | outside | methyl_mean | 16638 | 10114 | 0.5438 | 0.8303 | 0.2193 | -1.0842 |
| Paired | outside | methyl_high_fraction | 16638 | 10114 | 0.5185 | 0.8039 | 0.2261 | -1.0475 |
| Paired | outside | methyl_low_fraction | 16638 | 10114 | 0.4444 | 0.1224 | 0.7788 | 1.0895 |
| Paired | outside | region_methyl_mean | 90 | 86 | 0.5800 | 0.8210 | 0.2072 | -1.1439 |

## Full Feature Table

| Mode | Zone | Feature | n_TP | n_FP | TP median | FP median | AUC | Cohen's d |
|------|------|---------|------|------|-----------|-----------|-----|-----------|
| TO | inside | methyl_mean | 596 | 3133 | 0.7706 | 0.8447 | 0.4404 | -0.1201 |
| TO | inside | methyl_std | 596 | 3133 | 0.3695 | 0.3361 | 0.5292 | 0.1340 |
| TO | inside | methyl_median | 596 | 3133 | 0.9569 | 1.0000 | 0.4031 | -0.1882 |
| TO | inside | methyl_high_fraction | 596 | 3133 | 0.7064 | 0.8378 | 0.4185 | -0.2087 |
| TO | inside | methyl_low_fraction | 596 | 3133 | 0.1883 | 0.1455 | 0.5223 | 0.0463 |
| TO | inside | methyl_na_fraction | 596 | 3133 | 0.2500 | 0.3226 | 0.4501 | -0.1131 |
| TO | outside | methyl_mean | 15970 | 11618 | 0.5718 | 0.6975 | 0.3735 | -0.4521 |
| TO | outside | methyl_std | 15970 | 11618 | 0.4410 | 0.3925 | 0.6531 | 0.6400 |
| TO | outside | methyl_median | 15970 | 11618 | 0.9373 | 0.9725 | 0.4632 | -0.2017 |
| TO | outside | methyl_high_fraction | 15970 | 11618 | 0.5473 | 0.6578 | 0.3826 | -0.4205 |
| TO | outside | methyl_low_fraction | 15970 | 11618 | 0.4146 | 0.2634 | 0.6360 | 0.4837 |
| TO | outside | methyl_na_fraction | 15970 | 11618 | 0.2955 | 0.2930 | 0.5009 | -0.0235 |
| TO | inside | region_methyl_mean | 4 | 21 | 0.6922 | 0.8607 | 0.4167 | -0.0983 |
| TO | outside | region_methyl_mean | 86 | 77 | 0.6035 | 0.7291 | 0.3819 | -0.4284 |
| Paired | outside | methyl_mean | 16638 | 10114 | 0.5438 | 0.8303 | 0.2193 | -1.0842 |
| Paired | outside | methyl_std | 16638 | 10114 | 0.4483 | 0.3178 | 0.7804 | 1.1343 |
| Paired | outside | methyl_median | 16638 | 10114 | 0.8862 | 0.9882 | 0.3317 | -0.6062 |
| Paired | outside | methyl_high_fraction | 16638 | 10114 | 0.5185 | 0.8039 | 0.2261 | -1.0475 |
| Paired | outside | methyl_low_fraction | 16638 | 10114 | 0.4444 | 0.1224 | 0.7788 | 1.0895 |
| Paired | outside | methyl_na_fraction | 16638 | 10114 | 0.2807 | 0.2713 | 0.5055 | -0.0112 |
| Paired | outside | region_methyl_mean | 90 | 86 | 0.5800 | 0.8210 | 0.2072 | -1.1439 |

## Conclusions

1. **TO methyl_mean inside LOH**: AUC=0.4404, d=-0.1201 (n_TP=596, n_FP=3133) — **完全無區分力**
2. **TO region-level inside LOH**: AUC=0.4167, d=-0.0983 (n_TP=4, n_FP=21) — 樣本過少，不可靠
3. **Paired inside LOH**: 零個 TP read（HCC1395-only 在 paired 模式下 LOH.bed 內無 TP）
4. **先前 paired AUC=0.785 不可複現**: 來自 full-637 跨樣本 export（9 TP regions），HCC1395-only 無法重現
5. **TO TP 甲基化分佈**: median=0.771（高甲基化），與 paired TP median=0.463 完全不同
6. **Per-sample AUC (Simpson's paradox)**: HCC1395=0.614, HCC1395_DORADO=0.653（池化後降至 0.44）
7. **Potential_LOH 交叉檢查**: LOH.bed inside 3,729 reads 全部 Potential_LOH=False — ISM site-level 與 LOH.bed region-level 不一致

## GO/NO-GO Verdict

**NEGATIVE — methyl_mean 不可用於 TO LOH 區域 TP/FP 區分**

| 指標 | 數值 | 判定 |
|------|------|------|
| TO inside LOH AUC | 0.440 | FAIL (< 0.60) |
| Per-sample consistency | 2/2 < 0.66 | FAIL |
| Region-level AUC | 0.417 (n=25) | FAIL |
| 所有甲基化特徵 AUC | < 0.56 | FAIL |

## Implications

1. **methyl_mean 的「基因組背景」區分（TP 在活性染色質、FP 在靜默染色質）只在 paired 模式有效**，TO 模式因 self-phasing 導致 TP/FP 共享相同甲基化背景
2. **TO LOH 區域內所有甲基化方向正式關閉** — 與 O12、O15 結論一致
3. **QS 設計方向**: TO LOH 區域只能用 caller 特徵（AF, AD），甲基化權重必須歸零
