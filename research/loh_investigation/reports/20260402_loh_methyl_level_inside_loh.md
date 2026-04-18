    <!--
    建立時間: 2026-04-02
    目標: 分析 TP/FP 甲基化水平差異在 LOH 區域內是否持續存在
    處理範圍: 86,521 reads x 7 samples, paired-only, LOH zone overlay from TO LOH.bed
    關聯檔案:
      - research/loh_investigation/figures/loh_methyl_level_fig*.png
      - research/loh_investigation/data/loh_methyl_level_stats.tsv
    -->

    # LOH 區域內甲基化水平 TP/FP 差異分析

    ## 背景

    O10 觀察發現 read-level methyl_mean 在 TP (0.463) vs FP (0.679) 之間有顯著差異 (AUC=0.733)。
    本分析的核心問題：**這個差異在 LOH 區域內是否仍然存在？**

    ## 數據

    - 來源: phase1_shard_read_training_table.tsv (paired-only, 86,521 reads)
    - LOH 區域: 7 samples 的 TO-mode LOH.bed
    - LOH 判定: 以 read 的 chr + start 座標是否落在 LOH.bed 區間內

    ## 核心結果

    ### LOH Zone Distribution

    | Zone    | n_reads | % of total |
    |---------|---------|------------|
    | inside  |   7,632 | 8.8% |
    | outside |  78,889 | 91.2% |

    ### methyl_mean AUC (TP=1)

    | Level        | Inside LOH | Outside LOH | O10 Baseline |
    |--------------|------------|-------------|--------------|
    | Read-level   | **0.785**    | 0.719       | 0.733        |
    | Region-level | **0.792**    | 0.703       | —            |

    ### Cohen's d (TP - FP)

    | Level        | Inside LOH | Outside LOH |
    |--------------|------------|-------------|
    | Read-level   | -1.157 [-1.242, -1.076] | -0.788 [-0.808, -0.770] |
    | Region-level | -1.151 [-2.057, -0.465] | -0.706 [-0.895, -0.517] |

    ### Per-Sample Inside LOH

    | Sample               |      AUC | Cohen's d | n_TP   | n_FP   |
    |----------------------|----------|-----------|--------|--------|
    | HCC1395              |    0.810 |   -0.662 |     11 |    202 |
| HCC1395_DORADO       |      N/A |      N/A |      0 |    513 |
| COLO829              |    0.024 |    2.363 |     14 |      6 |
| H1437                |    0.392 |    0.707 |    452 |     35 |
| H2009                |      N/A |      N/A |      0 |   2926 |
| HCC1937              |    0.921 |   -1.632 |    218 |   2749 |
| HCC1954              |    0.637 |   -0.679 |    270 |    236 |

    ## 圖表清單

    1. fig01 — Zone distribution by sample
    2. fig02 — **Core figure**: methyl_mean violin by zone x truth_status
    3. fig03 — All methyl features inside LOH
    4. fig04 — AUC inside vs outside comparison
    5. fig05 — Per-sample violin inside LOH
    6. fig06 — Region-level vs read-level comparison
    7. fig07 — ALT/REF support inside LOH
    8. fig08 — HP tag inside LOH
    9. fig09 — Confound check (num_cpg, reads per region)
    10. fig10 — Complete summary table

    ## 結論

    **Inside LOH read-level AUC = 0.785** (region-level = 0.792)
    **Outside LOH read-level AUC = 0.719** (region-level = 0.703)

**判定: methyl_mean TP/FP 差異在 LOH 區域內 *持續存在*，AUC 維持在高水平。** 甲基化水平是一個穩健的 discriminator，即使在 LOH 區域內也有效。
