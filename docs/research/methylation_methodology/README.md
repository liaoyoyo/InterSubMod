<!--
建立時間: 2026-03-12 10:30
目標: 提供甲基方法學主線的研究入口、主文檔與邊界說明
處理範圍: docs/research/methylation_methodology/
關聯檔案:
  - docs/CURRENT_FOCUS.md
  - docs/plans/2026/03/20260307_純樣本甲基研究執行計畫_01.md
  - docs/reports/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.md
-->

# 甲基方法學主線入口

## 主題定位

`docs/research/methylation_methodology/` 是目前 InterSubMod 研究主線的長期方法學入口，重點在：

1. 為何先做純樣本、paired 與 tumor-only
2. 為何 `HCC1395 5kHz` 是主實驗樣本
3. `HCC1395_DORADO` 如何作交叉驗證
4. `caller-first / methylation-support / artifact-triage` 的方法學分工
5. `label-first / cluster-first` 的設計與驗證邏輯

## 主文檔

1. 方法學藍圖：
   [docs/research/methylation_methodology/2026/03/20260307_5kHz主實驗與方法學驗證藍圖_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/research/methylation_methodology/2026/03/20260307_5kHz主實驗與方法學驗證藍圖_01.md)

## 與其他文件的關係

1. 本目錄：
   - 長期方法學與主題主線
2. `docs/plans/`
   - 這一階段實際要做哪些任務
3. `docs/experiments/`
   - 每一輪實驗如何執行、結果如何
4. `docs/reports/`
   - 已整合的正式週報與 validated 結論

## 建議閱讀順序

1. [docs/CURRENT_FOCUS.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md)
2. [docs/experiments/INDEX.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/INDEX.md)
3. [docs/research/methylation_methodology/2026/03/20260307_5kHz主實驗與方法學驗證藍圖_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/research/methylation_methodology/2026/03/20260307_5kHz主實驗與方法學驗證藍圖_01.md)
4. [docs/plans/2026/03/20260307_純樣本甲基研究執行計畫_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260307_純樣本甲基研究執行計畫_01.md)
5. [docs/reports/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.md)

## 目前邊界

```text
docs/research/methylation_methodology/
├── README.md
└── 2026/03/
    └── 20260307_5kHz主實驗與方法學驗證藍圖_01.md
```

若未來此主題新增：

1. 方法學藍圖或長期研究策略，放這裡
2. 實驗 round 結果，放 `docs/experiments/`
3. 整合週報與結論，放 `docs/reports/`
