<!--
建立時間: 2026-04-12 19:00
目標: 文獻驗證研究工作區索引
-->

# Literature Validation — 文獻假說 vs ISM 實證數據

對 60+ 篇癌症甲基化文獻的核心假說進行系統性驗證。

## 研究結果

| 假說 | 來源 | 結果 | 報告 |
|------|------|------|------|
| L1 方向性 ASM | epiTRACERx, Gaiti | **NEGATIVE** | `reports/20260412_L1_directional_ASM_report.md` |
| L2 PMD 分層 | Berman, Landau | 生物 POSITIVE, 應用 **NEGATIVE** | `reports/20260412_文獻驗證綜合報告_01.md` §2 |
| L3 Normal baseline | Phase 2A | **FEASIBLE** | 同上 §3 |
| L4 fCpG 選擇 | EVOFLUx | **NEGATIVE** | 同上 §4 |

## 目錄

```
scripts/
├── L1_directional_ASM.py           # L1 主分析（Part A/B/C）
├── L1_supplementary_AF_stratified.py # L1 AF 分層補充
└── L2_L3_feasibility_check.py      # L2 PMD + L3 Normal BAM 可行性

figures/
├── 01_allele_delta_aggregate.png   # AlleleDelta TP vs FP
├── 02_directional_asm_detail.png   # 方向性 ASM + CpG variance
├── 03_alt_ref_scatter.png          # ALT vs REF 甲基化
├── 04_supplementary_af_stratified.png # AF 分層分析
└── 05_L2_pmd_proxy.png             # PMD proxy 分析

data/
├── part_a_aggregate.csv            # Per-sample aggregate 數據
└── part_b_directional_asm.csv      # 5,360 regions raw ASM 數據

reports/
├── 20260412_L1_directional_ASM_report.md
└── 20260412_文獻驗證綜合報告_01.md
```

## 數據規模

- **Part A**: 340,173 區域 × 7 樣本（significance_summary.csv aggregate）
- **Part B**: 5,360 區域抽樣（reads.tsv + methylation.csv raw reads）
- 全部在 paired_pileup complete_matrix canonical 輸出上執行
