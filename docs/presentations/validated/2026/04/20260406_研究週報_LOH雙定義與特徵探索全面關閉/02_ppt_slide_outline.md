<!--
建立時間: 2026-04-06 23:30
目標: PPTX 投影片段落主軸（13 slides）
處理範圍: 2026-03-31 ~ 2026-04-06
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/validated/2026/04/20260406_研究週報_LOH雙定義與特徵探索全面關閉/01_full_narrative_report.md
-->

# 投影片段落主軸

---

## Slide 1 — 封面

**標題**：InterSubMod 研究週報：LOH 雙定義與特徵探索全面關閉
**核心訊息**：本週完成 LOH 雙定義系統性驗證 + TO 特徵全面關閉，研究方向正式轉向

**要點**：
- 報告期間：2026-03-31 至 2026-04-06
- 主題：LOH 雙定義交叉分析 + 特徵探索全面關閉
- 166 圖表 × 16 判定 × 7 樣本

---

## Slide 2 — 本週重點結論

**標題**：本週重點結論
**核心訊息**：四個決定性結論改變研究走向

**要點**：
- SEQC2 Jaccard=0.928：LOH.bed 與 FDA 金標準幾乎完全吻合
- TO 特徵全面關閉：10/10 filter FAIL + Non-LOH max AUC<0.58 + Voting 0.577
- 甲基化三維度 NEGATIVE：O11-O13 全部因 confound 否決
- 戰略轉向：read-level epigenetic characterization，非 variant filter

---

## Slide 3 — 研究背景與四象限定義

**標題**：LOH 雙定義與四象限分析框架
**核心訊息**：兩種 LOH 定義建立四象限，58.1% 位點非 LOH

**要點**：
- 三層分析：LOH 定義可信度 → ISM 影響 → Filter 可行性
- 四象限：Q1 Both(26.7%) / Q2 HP-only(15.2%) / Q3 LOH.bed-only(0.07%) / Q4 Neither(58.1%)
- 本週範圍：166 圖表 × 16 項判定（J1-J16）
- 注意：所有數據為 HP fix 後新版

---

## Slide 4 — 本週目標與完成度

**標題**：主要目標與完成度
**核心訊息**：7 項目標全部 100% 完成

**要點**：
- LOH 雙定義 Wave 1+2+3：16 判定全完成
- 甲基化三維度 O11-O13：全 NEGATIVE
- TO VCF 特徵 G1-G7：60+ 特徵全 AUC<0.64
- Self-phasing 因果鏈：完整五步驗證 CONFIRMED
- QS mode-aware：已實作並部署
- ASM 定量：5 方法驗證 32-66%

---

## Slide 5 — 研究推進時間軸

**標題**：研究推進時間軸
**核心訊息**：一週內從 LOH 定義驗證到全面關閉

**要點**：
- 03-31：O11-O13 甲基化假說 NEGATIVE
- 04-01~02：G1-G7 NO-GO + Self-phasing 因果鏈 CONFIRMED
- 04-03~05：LOH 雙定義 Wave 1+2 + TO Deep Study
- 04-06：Wave 3 Non-LOH + 多特徵 + 全面關閉

---

## Slide 6 — SEQC2 外部驗證：LOH.bed 可信

**標題**：SEQC2 LOH 外部驗證：Jaccard=0.928
**核心訊息**：完全獨立的 FDA 金標準驗證 LOH.bed 高度可信

**要點**：
- SEQC2：FDA 多平台共識（WGS + array + 多中心），320 entries, 1,490.5 Mb
- LongPhase-TO：ONT long-read phased genotype — 方法學零重疊
- Jaccard=0.928, Sensitivity=0.961, Precision=0.964, F1=0.963
- HP_Ratio 作為 per-variant 預測器 AUC=0.8979

---

## Slide 7 — ISM 在 LOH 區域全面失效

**標題**：LOH 區域：PERMANOVA 失效 + 甲基化 AUC~0.50
**核心訊息**：LOH 區域只有一條 haplotype，ISM HP 比較無統計意義

**要點**：
- PERMANOVA valid rate 僅 5-6%（94% min(HP)<3）
- 7/7 樣本甲基化 AUC~0.50（隨機）
- 根因：LOH 只有一條 haplotype 有 reads
- 已修正：QS mode-aware（TO 停用 LOH penalty + verify bonus）

---

## Slide 8 — LOH 不可作 Filter + Non-LOH 無突破

**標題**：10/10 Filter FAIL + Non-LOH max AUC<0.58
**核心訊息**：LOH 是 TP-enriched；Non-LOH 同樣無有效特徵

**要點**：
- LOH 內 FP rate=0.239 < Non-LOH 0.338（LOH 內 TP 更多）
- 10 種策略全超 TP loss 2% 安全線
- Non-LOH：最高 corrected AUC 無法超過 0.58
- 多特徵 Voting AUC=0.577 < 0.58

---

## Slide 9 — Simpson's Paradox 與 Confound 警告

**標題**：cnLOH 陷阱：表面 0.587 是 Simpson's Paradox
**核心訊息**：整體數值被大樣本主導，per-sample 驗證揭露真相

**要點**：
- cnLOH PairwiseMeanDist overall AUC=0.5865
- Per-sample mean AUC=0.4987（本質隨機）
- H2009 佔 48.6% 主導整體數字
- CramersV：AUC 0.511→0.464 after NumReads residualize
- AlleleDelta 是唯一 confound-free 信號（AUC=0.556，不足）

---

## Slide 10 — 甲基化 NEGATIVE + TO 特徵 NO-GO + ASM 亮點

**標題**：甲基化三維度 NEGATIVE × G1-G7 NO-GO × ASM 唯一亮點
**核心訊息**：所有方向窮盡後，ISM 價值在 characterization 非 filtering

**要點**：
- O11：epipolymorphism 0.845→0.530（n_reads confound）
- O12：AlleleDelta=AF confound + L2 collider bias
- O13：shared read count confound
- G1-G7：48 圖 60+ 特徵全 AUC<0.64
- ASM 32-66%：ISM PERMANOVA 唯一量化 haplotype entropy imbalance

---

## Slide 11 — 14 結論總表

**標題**：結論總表：6 CONFIRMED + 5 NEGATIVE + 1 POSITIVE + 1 CONDITIONAL + 1 CONFIRMED
**核心訊息**：本週 +2 CONFIRMED +3 NEGATIVE，研究圖景大幅收斂

**要點**：
- C5 LOH.bed 可信 NEW
- C6 LOH 不可作 filter NEW
- C7 甲基化三維度 NEGATIVE
- C8 TO VCF 特徵 NO-GO
- C9 Non-LOH 無突破
- C12 多特徵組合不可行
- C14 LOH 價值在 stratification

---

## Slide 12 — 下週行動與風險

**標題**：下一步：Phase 2A 規劃與風險評估
**核心訊息**：PON-only phasing 全量重跑 + normal reference baseline

**要點**：
- 行動 1：Phase 2A PON-only phasing + ISM 重跑
- 行動 2：Normal methylation reference baseline
- 行動 3：研究全景 v2.0 更新
- 風險：PON-only 重跑後仍無改善（中）、低純度樣本不可得（低）

---

## Slide 13 — 複查路徑

**標題**：複查路徑與再生成方式
**核心訊息**：所有來源、數據、設定檔可直接回查

**要點**：
- 來源週報路徑
- LOH 雙定義分析報告路徑
- 研究全景索引
- PPTX config + 重生指令
