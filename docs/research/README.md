<!--
建立時間: 2026-03-12 10:30
目標: 提供 docs/research/ 的主題入口、邊界說明與查詢順序
處理範圍: docs/research/ 下的主題研究文檔，不含 scripts/ 與 repo 外結果輸出
關聯檔案:
  - docs/README.md
  - docs/research/methylation_methodology/README.md
  - docs/research/methylation_f1_optimization/README.md
-->

# 研究主題入口

## 目的

`docs/research/` 用來放 **長期研究主題** 的文檔入口、研究藍圖、方法學整理與主題式歷史脈絡。

這裡的重點不是每一次實驗的即時結果，而是：

1. 某個研究主題為什麼重要
2. 目前主線如何定義
3. 該主題有哪些正式研究文檔
4. 圖表資產、歷史腳本與正式報告各自放在哪裡

## 與其他目錄的邊界

1. `docs/research/`
   - 長期主題研究脈絡、方法學藍圖、主題式歷史整理
2. `docs/plans/`
   - 當前階段的執行計畫、任務拆解與驗收條件
3. `docs/experiments/`
   - 單次實驗與 round 的進行中/已驗證結果
4. `docs/reports/`
   - 已驗證或 finalized 的正式報告
5. `scripts/analysis/legacy/`
   - 歷史分析腳本，若仍需保留可追溯性但不屬於當前主流程

## 建議查詢順序

1. 先讀 [docs/CURRENT_FOCUS.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md)
2. 再讀 [docs/experiments/INDEX.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/INDEX.md)
3. 若要理解長期研究主題，再從本頁選對應 topic

## 目前主題

### 1. 甲基方法學主線

- 入口：
  [docs/research/methylation_methodology/README.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/research/methylation_methodology/README.md)
- 定位：
  純樣本優先、5kHz 主實驗、DORADO 交叉驗證、label-first / cluster-first 方法學驗證

### 2. 甲基顯著性與 F1 最佳化歷史研究

- 入口：
  [docs/research/methylation_f1_optimization/README.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/research/methylation_f1_optimization/README.md)
- 定位：
  早期 F1 優化研究與圖表、歷史腳本、舊有特徵分析整理

## 目錄原則

1. 每個 topic 目錄都應有 `README.md`
2. topic 目錄只保留：
   - 主題文檔
   - 最小必要資產
3. 若混入腳本或工作區容器：
   - 腳本移到 `scripts/`
   - 舊結構移到 `docs/archive/...pending_review.../`

## 目前 topic 結構

```text
docs/research/
├── README.md
├── methylation_methodology/
│   ├── README.md
│   └── 2026/03/*.md
└── methylation_f1_optimization/
    ├── README.md
    ├── 2026/01/*.md
    └── assets/2026/01/
```
