---
title: Archive TO 資料重跑 ISM 需求記錄
status: pending
owner: liaoyoyo2001
created: 2026-04-22
last_updated: 2026-04-22
priority: P1
related:
  - /big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/08_archive_to_crosssample.md
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/
---

# Archive TO 資料重跑 ISM 需求記錄

## 背景

於 2026-04-22 發現：`research_rounds/archive/202603_early_pilots/` 下已有 5 個樣本的 TO pipeline archive 結果（HCC1395_DORADO、H1437、H2009、HCC1937、HCC1954，COLO829 缺 step05）。這批資料之前**未被納入** canonical master（`output/canonical/README.md` 的 `to_pileup` 列對除 HCC1395 外皆填 `—`）。

本次（2026-04-22）做了**輕量彙整**到 `data/master_extended.tsv.gz`，但檔案**缺少 LOH/CN/KDE 相關欄位**（LOH_Subtype、Coverage_Multiple、CovM_used、baseline_x 等），限制了可做的分析：

- ✅ 可做：L1 AF × NG × NumReads 跨樣本觀察
- ❌ 無法做：完整 L2/L3（涉及 LOH/CN）、Top-17 filter 跨樣本驗證、baseline_x 統一標準化

要完整驗證 ISM-based TP/FP filter 的跨樣本泛化性，**必須重跑 ISM pipeline**（HP-fix 後版本 + KDE baseline）。

---

## 為何 archive 資料缺 LOH/CN 欄位

Archive 執行時間（2026-03-15 至 2026-03-30）**早於**以下 pipeline 升級：

| 升級項目 | 日期 | 影響 |
|---------|------|------|
| KDE Coverage baseline | ~2026-04 | `CovM_used`, `baseline_x`, `kde_status` 欄位誕生 |
| LOH.bed 標註（ISM output）| ~2026-04 | `LOH_Subtype`, `LOH_Bed_Overlap` 欄位誕生 |
| cn_tier_F 定義 | ~2026-04 | CN tier 對應邏輯改版 |
| HP Integer Tag Fix | 2026-03-30 | archive 部分已是 post-fix（`intersubmod_*` 非 `_before_hp_fix`）|

因此 archive 的 `label_first_metrics.tsv` 雖然 HP tag 正確，但**不含**後續 KDE/LOH 衍生欄位。若要做 `S1-S5 scheme` 或 `Top-17 cell` 的跨樣本驗證，必須用**當前 pipeline** 重跑。

---

## 需要重跑 ISM 的樣本與資源需求

### 輸入資料可得性（2026-04-22 確認）

| 樣本 | Tumor tagged BAM (archive step03) | Normal BAM | Truth VCF | 估計 ISM 時間 |
|------|:--:|:--:|:--:|:--:|
| HCC1395_DORADO | ✅ `20260315_hcc1395_dorado_to_pilot/step03_longphase_to/tumor_tagged.bam` | - (TO 不需) | ✅ | ~2-3 hr |
| H1437 | ✅ `20260318_h1437_to_pilot_fastresume/...` | - | ✅ | ~3-4 hr |
| H2009 | ✅ `20260318_h2009_to_pilot_fastresume/...` | - | ✅ | ~8-10 hr（large） |
| HCC1937 | ✅ `20260318_hcc1937_to_pilot_fastresume/...` | - | ✅ | ~2-3 hr |
| HCC1954 | ✅ `20260318_hcc1954_to_pilot_fastresume/...` | - | ✅ | ~4-5 hr |
| COLO829 | ❌ **step05 empty** | - | ✅ | 需重跑 step01-03（~2 hr）+ ISM（~3 hr）|

**總估時（序列）**：~25-30 小時；可 5 樣本 parallel 7-10 小時

### 所需管線步驟

對每樣本：
1. **輸入**：archive 的 `step03_longphase_to/tumor_tagged.bam`（已 HP-tagged）+ archive 的 TP/FP split VCF
2. **ISM**：當前 `run_vcf_all_snv.sh --mode all-with-w5000` + KDE enabled + LOH.bed 標註
3. **輸出**：每 region 有 `LOH_Subtype`, `Coverage_Multiple`, `CovM_used`, `baseline_x` 欄位
4. **合併 master**：重跑 `tpfp_loh_af_kde_discrimination` 的 `step4/4b/6/8` → 擴充 master.tsv.gz

### 前置條件

- ✅ ISM pipeline 已穩定（HP-fix 完成）
- ✅ KDE baseline 已啟用（`expected_coverage` bug 已修?需確認）
- ✅ LOH.bed 產生流程已加入 ISM 輸出

---

## 決策點（留給 PI/主研究者）

**Path 1（全量重跑，~10h parallel）**：
- 跑 6 樣本 → 完整 L3 LOH×AF×CN 跨樣本觀察 → 驗證 Top-17 cell 跨樣本泛化性
- 時間估：7-10 hr（5-6 cores parallel）
- 收益：完整回答「TO filter 跨樣本穩定嗎」

**Path 2（優先 3 樣本，~4h）**：
- COLO829（原缺）+ HCC1395_DORADO（重要對照）+ H2009（最大樣本 + 最高 TP%）
- 時間估：4-5 hr
- 收益：covers 最具差異的 3 種 TO 行為

**Path 3（僅 COLO829）**：
- 先補最缺的一個樣本 → 決定後續是否擴大
- 時間估：~3 hr
- 收益：至少 master 有第 2 個 TO 樣本

### 建議

**Path 2 優先**。理由：
1. COLO829 TO 完全空，補齊後 master 有 2 個 paired+TO 對照組（HCC1395、COLO829）。
2. HCC1395 vs DORADO 的 AF 分佈差異（80.9% vs 0.4% Extreme）顯示 basecaller 版本衝擊巨大 → DORADO 是關鍵對照，重跑 ISM 後可驗證 LOH/CN 分佈是否也跟著 basecaller 變動。
3. H2009 是最大 TO 樣本（137k rows）且 baseline 最高（91%），是 characterization 的重要正面案例。

---

## 補充：為什麼這批 archive 之前被忽略

三層遺忘機制（見 Memory：`project_to_cross_sample_archive_data_exists.md`）：

1. **Canonical README 的 `—`**：`output/canonical/README.md` 的 `to_pileup` 列對 6 樣本寫 `—`，暗示「無資料」，實際應為 `archive-only`。
2. **master_run_manifest 只登記 1 TO run**：其他 6 個 archive pilot 未登記，JSON 視覺與全量 CLI 搜不到。
3. **CURRENT_FOCUS 的定性敘事**：「TO 模式 ISM 無效」的過度斷言掩蓋了「archive 有跨樣本資料但 LOH/CN 欄位過時」這個真實狀態。

**修正動作（本次完成）**：
- ✅ 建立本 .md 記錄
- ✅ 建立 `08_archive_to_crosssample.md` 觀察檔
- ⏳ 更新 `output/canonical/README.md` 將 `—` 改為 `archive-only` + 路徑（下次對話）
- ⏳ 更新 Memory `project_to_cross_sample_archive_data_exists.md`
- ⏳ 更新 `CURRENT_FOCUS.md` 第 60/129/221 行 "TO 模式 ISM 無效" 為 "TO 模式 ISM filter 跨樣本泛化性尚未驗證（需重跑 ISM）"

---

## Exit Criteria

本需求完成當且僅當：
1. 至少完成 Path 2（COLO829 + DORADO + H2009）的 ISM 重跑
2. `master.tsv.gz`（或 `master_extended.tsv.gz`）包含新樣本的 LOH_Subtype / CovM_used
3. 重跑 `obs09-12` 產生跨樣本 L1-L3 + consistency 圖
4. `04_comparison_narrative.md` 更新 Top-17 表增加跨樣本欄位（n_samples_n20, n_samples_high）
