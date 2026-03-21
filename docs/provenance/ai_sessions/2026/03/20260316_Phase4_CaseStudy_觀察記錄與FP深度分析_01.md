# Phase 4 Case Study 觀察記錄與 FP 深度分析 執行報告

<!--
建立時間: 2026-03-16 13:14
目標: 完成 Phase 4 case study 文件的觀察記錄填寫、HP tag 修正、FP-B1/B2 深度分析、知識點整理
處理範圍: 5 TP + 4 FP 案例的觀察/驗證/推論記錄；FP-B1 SEQC2 INDEL 關聯；FP-B2 MNP 機制
關聯檔案:
  - docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_01.md
  - scripts/utils/igv_batch_screenshot.sh
  - scripts/analysis/screen_mnp_adjacent_fp.py
  - docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_assets/igv_positions.tsv
-->

## 對話資訊

| 項目 | 內容 |
|------|------|
| 日期 | 2026-03-16 |
| 主要目的 | 填入 Phase 4 case study 文件的觀察記錄；修正 HP tag 定義；深度分析 FP-B1 SEQC2 INDEL 關聯與 FP-B2 MNP 機制 |
| AI 模型 | Claude Sonnet 4.6 |
| 資料來源 | HCC1395 5kHz Paired (s-pure/HCC1395/20260307)，9 個 case study 位點 |

## 對話背景

接續前次對話（20260315 Phase 4 建構）的後半段工作：

1. 使用者提供了所有 9 個案例（5 TP + 4 FP）的詳細觀察，要求以 觀察/驗證/推論 三種分類記錄
2. 使用者發現 HP tag 定義錯誤（HP=3 非「第三 phase group」，而是「LongPhase-S 無法區分 HP1/HP2 但跨越 Somatic ALT」）
3. 使用者在 IGV 觀察到 FP-B1 (chr7:52087777) 鄰近有 SEQC2 INDEL TP
4. 使用者在 IGV 觀察到 FP-B2 (chr9:75383880) 的每個 ALT reads 在前一位置有一致的 indel-like mismatch
5. 要求整理可重複應用的知識點，並建立 IGV 自動化截圖工具

## 關鍵決策

| 決策 | 原因 | 影響 |
|------|------|------|
| HP=3 定義修正 | 使用者確認：HP=3 = LongPhase-S 無法區分 HP1/HP2 但 read 跨越 Somatic ALT；非「第三 phase group」 | 修正了文件中所有相關 HP 說明 |
| HP=1-1/2-1 定義修正 | 非「chimeric tag（跨越兩個 phased block）」，而是「HP1/HP2 背景的 read 帶有 Somatic ALT」 | 更正了對 ALT reads 的 HP 解讀邏輯 |
| FP-B2 確認為 MNP 非 indel | BAM 分析顯示 ALT reads 在 75383879 有一致 A（應為 C），CIGAR 無 indel；normal 無此鹼基 | 正確解釋了 IGV 視覺現象；計劃大規模 MNP 篩選 |
| FP-B1 識別為 SEQC2 INDEL 鄰近 | chr7:52087776 有 PASS+HighConf INDEL TP，chr7:52087777 SNV 緊鄰；SNV-only benchmark 無法正確評估 | 新增 benchmark annotation gap 分析章節 |
| 知識點獨立章節 | 使用者要求「整理可重複觀察與推論的知識點」 | 新增 K1-K6 知識點表格（HP 讀法、FP 鑑別、MNP 辨識等） |
| IGV 腳本建立 | 使用者要求基於 template.xml 的批次截圖方案 | 新增 `igv_batch_screenshot.sh` 腳本 |

## 產出成果

### 修改檔案
- `docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_01.md`
  - 修正「資料狀態說明」HP tag 定義表（HP=0/1/2/1-1/2-1/3 全部重新定義）
  - 新增「觀察記錄（TP 案例）」章節：5 個 TP 案例各含 觀察/已驗證/推論 三類記錄
  - 新增「FP 觀察記錄」章節：4 個 FP 案例各含 觀察/已驗證/推論 + 與 TP 的視覺差別
  - 新增「FP-B1 深度分析」章節：SEQC2 INDEL 關聯 + benchmark gap 說明 + pipeline 啟示
  - 新增「FP-B2 深度分析」章節：MNP 機制說明 + BAM 驗證結果 + 大規模驗證計劃
  - 新增「知識點整理」章節（K1-K6）：可重複應用的觀察規則表格
  - 新增「IGV 自動化截圖方案」章節：bash 腳本說明（方案一/二）
  - 修改「下一步」章節：新增大規模 MNP 篩選、IGV 環境確認、FP-B1 INDEL benchmark 重評估

### 新增檔案
- `scripts/utils/igv_batch_screenshot.sh` — 基於 template.xml 的 per-locus IGV 截圖腳本（支援單一位點與批次模式）
- `scripts/analysis/screen_mnp_adjacent_fp.py` — 大規模 MNP FP 篩選腳本（全 627 FP 位點掃描相鄰 mismatch）
- `docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_assets/igv_positions.tsv` — 9 個案例位點的 IGV 批次截圖輸入清單

## 核心發現（觀念更新）

| 更新前認知 | 更新後（本輪發現） |
|---|---|
| HP=3 是「第三個 phase group」暗示 LOH | HP=3 表示 LongPhase-S **無法區分 HP1/HP2**，但 read 確實跨越 Somatic ALT；與 LOH 無直接對應 |
| HP=1-1 是「chimeric tag，跨越兩個 phased block」 | HP=1-1 是「HP1 背景的 read，但同時帶有 Somatic ALT allele」 |
| FP-B2 chr9:75383880 的 ALT reads 有 indel | 實際是 MNP（C>A at 75383879 + T>A at 75383880），CIGAR 無 indel；IGV 的「indel-like」外觀來自 adjacent mismatch |
| FP-B1 chr7:52087777 是純粹的 FP | 緊鄰 chr7:52087776（SEQC2 INDEL TP，TA>T，PASS+HighConf）；SNV-only benchmark 造成 annotation gap |
| CramersV=1.0 一定是真實 somatic 信號 | HP-driven germline imprinting（FP-B1）也可達 CramersV=1.0；需同時確認 AlleleDelta 和 HP_Ratio |

## 後續行動

- [ ] **大規模 MNP 篩選**：執行 `screen_mnp_adjacent_fp.py` 掃描全 627 FP 位點
  - 使用 FP VCF + tumor haplotagged BAM + normal BAM
  - 預期產出：`analysis/mnp_adjacent_fp_screen_20260316.tsv`
- [ ] **IGV 環境確認**：
  - `which igv.sh || ls /opt/igv/igv.sh` — 確認 IGV 安裝
  - `which xvfb-run` — 確認 xvfb 可用
  - 測試：`./scripts/utils/igv_batch_screenshot.sh chr9 75383880 2000 /tmp/test/`
- [ ] **IGV 批次截圖**：執行 9 個位點的截圖
  - 輸入：`docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_assets/igv_positions.tsv`
  - 需先修改 tsv 格式（去掉 comment 行 header）
- [ ] **FP-B1 benchmark 重評估**：用含 sINDEL truth set 的 benchmark 重評估 chr7:52087777
- [ ] **HCC1395 DORADO 交叉驗證**：確認 FP-B2 MNP pattern 在 DORADO BAM 中的重現性

## 對話摘要

本輪接續前次對話，完成 Phase 4 case study 文件的觀察記錄填寫工作。核心貢獻：(1) 修正了 HP tag 定義（HP=3 = ALT+unresolved_HP，HP=1-1/2-1 = HP1/2背景+ALT）；(2) 為 9 個案例填入觀察/驗證/推論分類記錄；(3) 確認 FP-B2 是 MNP（非 indel），C>A+T>A 在 HP1 haplotype 上，CIGAR 無 indel 但 IGV 顯示 adjacent mismatch；(4) 確認 FP-B1 鄰近 SEQC2 INDEL TP，是 SNV-only benchmark 的 annotation gap 案例；(5) 新增 K1-K6 知識點表格整理可重複應用的判斷規則；(6) 建立 igv_batch_screenshot.sh 和 screen_mnp_adjacent_fp.py 兩個自動化工具。
