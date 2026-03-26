# 亞克隆地圖 Pilot 報告 (20260315)

## 1. Subclone 位點概況

| 項目 | 數值 |
|------|------|
| 總 Subclone 位點 | 1,412 |
| TP 中的 Subclone | 1,402 |
| FP 中的 Subclone | 10 |

## 2. 各染色體分布

| Chr | Subclone Count | Density (per Mbp) | Strong Count |
|---|---|---|---|
| chr1 | 68 | 0.2731 | 1 |
| chr2 | 119 | 0.4913 | 1 |
| chr3 | 163 | 0.822 | 0 |
| chr4 | 109 | 0.5731 | 3 |
| chr5 | 88 | 0.4847 | 1 |
| chr6 | 72 | 0.4215 | 2 |
| chr7 | 58 | 0.364 | 2 |
| chr8 | 156 | 1.0748 | 1 |
| chr9 | 39 | 0.2818 | 0 |
| chr10 | 90 | 0.6726 | 0 |
| chr11 | 104 | 0.7699 | 0 |
| chr12 | 104 | 0.7803 | 0 |
| chr13 | 52 | 0.4547 | 0 |
| chr14 | 27 | 0.2522 | 1 |
| chr15 | 39 | 0.3824 | 0 |
| chr16 | 2 | 0.0221 | 1 |
| chr17 | 58 | 0.6966 | 0 |
| chr18 | 25 | 0.3111 | 0 |
| chr19 | 20 | 0.3412 | 1 |
| chr20 | 8 | 0.1241 | 0 |
| chr21 | 5 | 0.107 | 1 |
| chr22 | 6 | 0.1181 | 0 |

## 3. 後續行動

- Phase 4 Case Study: 選取 CramersV>0.3 AND NumReads>=30 的 top 20-50 個位點
  - methylartist SVG 視覺化
  - modbamtools HTML 互動驗證
- 與 Phase 2 辨識特徵結果交叉分析

## 4. 資料來源

- Feature matrix: `analysis/feature_matrix_20260315.tsv`
- 來源 run: `s-pure/HCC1395/20260307` (新欄位 schema)
