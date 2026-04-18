<!--
建立時間: 2026-04-05 23:30
目標: Step 4 — 建立可驗證的癌症生物學解釋框架
處理範圍: 文獻假說建立、觀察驗證、未來研究方向確認
關聯檔案:
  - docs/plans/2026/04/20260405_LOH定義與多維標籤交叉分析研究計劃/00_總覽與執行順序.md
  - docs/reports/research_landscape/00_INDEX.md
  - docs/architecture/20260327_InterSubMod研究願景定錨_01.md
-->

# Step 4: 生物學解釋框架計畫

## 背景

ISM 的定位不只是 variant filter，而是「read-level epigenetic context for variant interpretation」。區分 TP/FP 提升 F1 是其中一種功能；用甲基化與特徵分析數據並以癌症生物學邏輯解釋是另一個核心功能（且是前者的基礎）。

先前 O1-O13 的教訓：看似有效的特徵（epipolymorphism AUC=0.845, L2 AUC=0.80）都被證實為 confound。**必須先有生物學解釋框架，才能判斷真正的信號與統計偽像的差異。**

### 為什麼需要這個步驟

1. **避免重蹈覆轍** — O11 的 epipolymorphism 在校正 n_reads 後 AUC 從 0.845 驟降至 0.530；O12 的 L2 collider bias 產生虛假信號。這些都是缺乏生物學錨點的後果。
2. **建立可引用的推論鏈** — 論文需要從生物學機制出發解釋觀察結果，而非僅報告統計數字。
3. **指導後續研究方向** — Phase 2A（Normal Methylation Reference）和 TO pipeline 的設計決策依賴於此框架的結論。

---

## 子任務 4.1: 文獻假說框架 (P2)

**目標**：建立 4 個可驗證的生物學假說，每個假說具備完整的文獻依據、可操作的預測、以及明確的判定標準。

### 假說總覽

| 假說 | 生物學依據 | 預測 | 所需數據 |
|------|-----------|------|----------|
| H1: LOH → ASM 消失 | LOH 區域只剩一個 allele，cis-regulatory 差異消失 | LOH 內 ASM 比例 < 非 LOH | ASM 定量 × LOH 分層 |
| H2: cnLOH 保留甲基化異質性 | 兩份同源拷貝可能有 epigenetic drift | cnLOH 的 PairwiseMedianDist ≈ non-LOH | Step 3.2 |
| H3: Self-phasing → 假 HP×methylation 關聯 | ALT reads 被推向一個 HP → 假 correlation | TO ALT_HP_skew 與 HPMergedDelta 正相關 | Step 2.2 |
| H4: Second-hit = LOH + methylation silencing | 雙等位基因失活模型 | LOH + low methylation 的 TP enrichment 更高 | 需 gene-level 數據 |

---

### H1: LOH 導致等位基因特異性甲基化 (ASM) 消失

**生物學依據**：
- 等位基因特異性甲基化 (Allele-Specific Methylation, ASM) 源於兩個等位基因的 cis-regulatory element 差異（Kerkel et al., 2008; Schalkwyk et al., 2010）
- LOH 事件消除了其中一個等位基因，使得 cis-regulatory 差異不復存在
- 因此 LOH 區域理論上應該呈現均一的甲基化狀態
- 參考文獻：Jenkinson et al., 2020 (epipolymorphism quantification); Do et al., 2020 (ASM in cancer)

**可驗證的預測**：
- LOH 區域內的 ASM 比例（以 AlleleDelta 或 HPMergedDelta 定量）應顯著低於非 LOH 區域
- 具體數值預測：LOH 內 ASM 比例 < 非 LOH 的 50%

**使用的數據來源**：
- Step 1.2 的 LOH 分佈與密度數據
- Step 2.1 的 HP-甲基化交叉分析結果
- O5/O12 已有的 ASM 定量結果（5 方法驗證 ASM 32-66%）

**判定標準**：

| 結果 | 判定 | 後續動作 |
|------|------|---------|
| LOH 內 ASM 比例 < 非 LOH 的 50%，且效應量 Cohen's d > 0.5 | CONFIRMED | 納入解釋框架，支持 LOH-aware 分析模式 |
| LOH 內 ASM 比例略低但差異不顯著（d < 0.3） | INCONCLUSIVE | 需增加樣本或細化 LOH 定義 |
| LOH 內外 ASM 比例無差異或 LOH 內反而更高 | REJECTED | 檢視 LOH 定義是否正確（是否混入 cnLOH），或 ASM 定量方法是否受 confound |

**潛在陷阱**：
- cnLOH 可能保留甲基化差異，混入分析會稀釋效應（見 H2）
- Self-phasing 可能產生假 ASM（見 H3），需區分 paired 與 TO

---

### H2: cnLOH 保留甲基化異質性

**生物學依據**：
- 拷貝中性 LOH (copy-neutral LOH, cnLOH) 指一個等位基因丟失後另一個等位基因倍增，維持總拷貝數不變
- 倍增的同源拷貝在 DNA 複製過程中可能經歷 epigenetic drift，導致兩份拷貝的甲基化狀態分歧
- 參考文獻：O'Keefe et al., 2010 (cnLOH in cancer); Berman et al., 2012 (epigenetic drift in cancer)

**可驗證的預測**：
- cnLOH 區域的 PairwiseMedianDist 應接近 non-LOH 區域，而非 deletion-LOH
- cnLOH 區域的 reads 數量應正常（不像 deletion-LOH 那樣減少）

**使用的數據來源**：
- Step 3.2 的 cnLOH 識別結果（深度分析區分 deletion-LOH 與 cnLOH）
- Step 1.3 的多維標籤交叉數據

**判定標準**：

| 結果 | 判定 | 後續動作 |
|------|------|---------|
| cnLOH 的 PairwiseMedianDist 與 non-LOH 無顯著差異（p > 0.05），且與 deletion-LOH 有顯著差異 | CONFIRMED | cnLOH 應作為獨立類別處理，不套用 LOH penalty |
| cnLOH 行為介於 deletion-LOH 與 non-LOH 之間 | INCONCLUSIVE | 可能需要更精細的 cnLOH 分型 |
| cnLOH 行為與 deletion-LOH 一致 | REJECTED | epigenetic drift 在短期內不足以產生可觀測差異，或 cnLOH 識別有誤 |

**潛在陷阱**：
- cnLOH 的識別依賴深度資訊，ONT 的覆蓋度變異可能影響判斷
- 需要足夠的 cnLOH 位點數量才能進行統計檢驗

---

### H3: Self-phasing 產生假 HP × methylation 關聯

**生物學依據**：
- 在 Tumor-Only (TO) 模式中，somatic SNV 被用於 phasing，導致 ALT reads 系統性地被分配到同一個 haplotype（self-phasing circular dependency）
- 這種偏差已被量化確認：62% LOH 消失、31.2% self-phasing LOH、somatic bias 17.3:1
- 結果是 HP 與 ALT allele 之間產生人工相關，進而與甲基化狀態產生假關聯
- 參考文獻：專案內部因果鏈驗證（project_self_phasing_causal_chain_confirmed.md）

**可驗證的預測**：
- TO 模式中 ALT_HP_skew（ALT reads 在兩個 HP 之間的不平衡）與 HPMergedDelta（HP 間甲基化差異）應呈正相關
- 此相關在 Paired 模式中應不存在或顯著減弱（因 paired 模式使用 germline SNP phasing）
- 移除 somatic SNV 參與 phasing 後（PON-only phasing），此相關應消失

**使用的數據來源**：
- Step 2.2 的 self-phasing 影響定量
- Step 2.1 的 HP-甲基化交叉分析（Paired vs TO 對比）

**判定標準**：

| 結果 | 判定 | 後續動作 |
|------|------|---------|
| TO: r(ALT_HP_skew, HPMergedDelta) > 0.3 且 Paired: r < 0.1 | CONFIRMED | HPMergedDelta 在 TO 模式中不可信，需使用 PON-only phasing 或移除此特徵 |
| 兩者皆有相關但 TO 更強（差異 > 0.15） | INCONCLUSIVE | 可能存在真實的 allele-specific methylation 效應疊加 self-phasing artifact |
| 兩者相關程度相似 | REJECTED | self-phasing 可能不是主要驅動因素，需重新檢視因果模型 |

**潛在陷阱**：
- Paired 模式也可能有殘餘的 phasing bias（germline heterozygous sites 分佈不均）
- 需控制 AF 作為混淆因子（高 AF 的 SNV 更容易產生 HP skew）

---

### H4: Second-hit = LOH + methylation silencing（雙等位基因失活）

**生物學依據**：
- Knudson 的二次打擊假說 (Knudson, 1971)：腫瘤抑制基因的失活需要兩個等位基因都被破壞
- 第一次打擊可以是突變或 LOH，第二次打擊可以是甲基化沉默 (methylation silencing)
- 這種模式在多種癌症中已被報導：BRCA1 (乳癌)、MLH1 (大腸癌)、VHL (腎癌)
- 參考文獻：Knudson, 1971 (two-hit hypothesis); Jones & Baylin, 2002 (epigenetic silencing); Esteller et al., 2001 (BRCA1 methylation)

**可驗證的預測**：
- 若 LOH 區域內的 TP somatic SNV 伴隨高甲基化（promoter 區域）或低甲基化（gene body），TP enrichment 應更高
- 具體預測：LOH + 異常甲基化（相對於 normal reference）的區域中 TP/FP ratio 應高於 LOH + 正常甲基化

**使用的數據來源**：
- 需要 gene-level annotation（CpG island, promoter, gene body）
- Step 1.2 的 LOH 區域 × TP/FP 分層
- 外部資源：RefSeq gene annotation, CpG island annotation

**判定標準**：

| 結果 | 判定 | 後續動作 |
|------|------|---------|
| LOH + 異常甲基化區域 TP/FP ratio 顯著高於其他組合 (OR > 2.0) | CONFIRMED | 支持將 ISM 定位為「二次打擊偵測工具」，論文核心賣點 |
| 趨勢存在但效應弱 (OR 1.2-2.0) | INCONCLUSIVE | 可能受限於樣本量或 gene annotation 精度 |
| 無差異或反向 | REJECTED | 二次打擊模型在此分析解析度下不適用，或甲基化變化非 driver 而是 passenger |

**潛在陷阱**：
- 需要 gene-level 解析度，但目前 ISM 以 SNV 為中心的分析窗口（通常 1-5 kb）可能不足以涵蓋完整 promoter 區域
- Normal methylation reference 缺失（Phase 2A 尚未完成）會限制「異常甲基化」的定義
- 此假說的完整驗證可能需要延後至 Phase 2A 完成後

---

### 假說間的關係與依賴

```
H1 (LOH → ASM 消失)
 ├── 如果 CONFIRMED → 強化 H2 的必要性（cnLOH 是否例外？）
 └── 如果 REJECTED → 需重新檢視 LOH 定義或 ASM 定量方法

H2 (cnLOH 保留異質性)
 └── 依賴 Step 3.2 的 cnLOH 識別結果

H3 (Self-phasing → 假關聯)
 ├── 與 H1 互補：如果 H1 CONFIRMED，TO 中看到的 ASM 可能是 H3 的 artifact
 └── 驗證結果直接影響 TO pipeline 設計

H4 (二次打擊)
 ├── 依賴 H1 的結果作為前提
 └── 需要額外的 gene annotation 資源，可部分延後
```

---

## 子任務 4.2: 觀察驗證 (P3)

**目標**：用 Steps 1-3 的數據系統性驗證每個假說。

### 驗證流程

對每個假說 (H1-H4) 執行以下結構化驗證：

#### 4.2.1 H1 驗證：LOH 內 ASM 是否消失

**數據來源**：
- Step 1.2 的 LOH 分佈數據（LOH.bed 重疊判定）
- Step 2.1 的 HP-甲基化交叉分析結果
- 已有的 O5 ASM 定量結果（AlleleDelta, HPMergedDelta）

**分析步驟**：
1. 將所有 SNV 分為 LOH / non-LOH 兩組
2. 計算兩組的 ASM 指標分佈（AlleleDelta, HPMergedDelta）
3. 統計檢驗：Wilcoxon rank-sum test + Cohen's d 效應量
4. 分層分析：按 Paired / TO 分別執行（因 H3 的 confound）
5. 繪圖：violin plot 或 box plot，標註 p-value 和效應量

**預期產出**：
- 圖表：`figures/H1_ASM_by_LOH_status.png`
- 數據表：`data/H1_ASM_LOH_statistics.tsv`
- Verdict: CONFIRMED / REJECTED / INCONCLUSIVE

#### 4.2.2 H2 驗證：cnLOH 的甲基化異質性

**數據來源**：
- Step 3.2 的 cnLOH 識別結果
- PairwiseMedianDist 特徵

**分析步驟**：
1. 將 LOH 進一步分為 deletion-LOH / cnLOH（基於 Step 3.2 的深度分析）
2. 比較三組（non-LOH, deletion-LOH, cnLOH）的 PairwiseMedianDist 分佈
3. 成對統計檢驗：Kruskal-Wallis + 事後 Dunn's test
4. 繪圖：三組 violin plot

**注意事項**：
- cnLOH 位點數量可能不足，需報告統計檢力 (statistical power)
- 如果 Step 3.2 無法可靠區分 deletion-LOH 與 cnLOH，此驗證標記為 INCONCLUSIVE 並記錄原因

**預期產出**：
- 圖表：`figures/H2_cnLOH_methylation_heterogeneity.png`
- Verdict: CONFIRMED / REJECTED / INCONCLUSIVE

#### 4.2.3 H3 驗證：Self-phasing 假關聯

**數據來源**：
- Step 2.2 的 self-phasing 影響定量
- ALT_HP_skew 與 HPMergedDelta 的逐位點數據

**分析步驟**：
1. 計算每個 SNV 的 ALT_HP_skew = |HP1_ALT - HP2_ALT| / (HP1_ALT + HP2_ALT)
2. 計算 r(ALT_HP_skew, HPMergedDelta)，分 Paired / TO
3. 散點圖 + 回歸線，分面板顯示 Paired 與 TO
4. 如果有 PON-only phasing 數據：加入第三個面板做對比

**預期產出**：
- 圖表：`figures/H3_self_phasing_false_correlation.png`
- 相關係數表：`data/H3_correlation_paired_vs_TO.tsv`
- Verdict: CONFIRMED / REJECTED / INCONCLUSIVE

#### 4.2.4 H4 驗證：二次打擊模型（初步）

**數據來源**：
- LOH × TP/FP 分層數據
- 甲基化水平（per-SNV 周圍平均甲基化）
- Gene annotation（如可取得）

**分析步驟**：
1. 初步版本（不需 gene annotation）：
   - 將 SNV 分為 4 組：LOH+高甲基化, LOH+低甲基化, non-LOH+高甲基化, non-LOH+低甲基化
   - 計算各組的 TP/FP ratio
   - 2x2 交互作用分析
2. 完整版本（需 gene annotation，可延後至 Phase 2A）：
   - 疊加 CpG island / promoter / gene body 分層
   - 聚焦已知的 tumor suppressor genes

**預期產出**：
- 圖表：`figures/H4_second_hit_preliminary.png`（2x2 bar plot）
- Verdict: 初步版本可能為 INCONCLUSIVE，完整驗證延後

---

### 驗證報告格式

每個假說的驗證報告遵循以下格式：

```markdown
### H[N] 驗證結果

**假說**：[一句話描述]

**數據摘要**：
- 樣本量：N = xxx（LOH: xxx, non-LOH: xxx）
- 關鍵統計量：[效應量、p-value]

**圖表引用**：
- 圖 [N]：[圖表描述]（`figures/H[N]_xxx.png`）

**Verdict**：CONFIRMED / REJECTED / INCONCLUSIVE

**推論**：
- 如果 CONFIRMED：[此結果對 ISM 的具體意義]
- 如果 REJECTED：[是方法問題還是假說錯誤？需要什麼修正？]
- 如果 INCONCLUSIVE：[缺乏什麼數據？需要什麼額外分析？]
```

---

## 子任務 4.3: 方向確認 (P3)

**目標**：基於 H1-H4 的驗證結果，確認後續研究方向，與五個研究目標對齊。

### 需要回答的四個關鍵問題

#### Q1: Phase 2A (Normal Methylation Reference) 是否應該包含 LOH 修正？

**決策邏輯**：
- 如果 H1 CONFIRMED（LOH → ASM 消失）：
  - Normal reference 需區分 LOH / non-LOH 區域
  - LOH 區域的 normal reference 應來自單一 allele 的甲基化水平
- 如果 H1 REJECTED：
  - Normal reference 不需特別處理 LOH
  - 但需理解為什麼 LOH 內的 ASM 沒有消失

**輸出**：Phase 2A 設計規格的 LOH 修正章節

#### Q2: ISM 是否需要 LOH-aware 分析模式？

**決策邏輯**：
- 如果 H1 + H2 的結果顯示 LOH 類型（deletion-LOH vs cnLOH）顯著影響甲基化行為：
  - ISM 需要在分析前先判斷每個 region 的 LOH 狀態
  - 不同 LOH 類型使用不同的 distance metric 或 threshold
- 如果 LOH 類型的影響可忽略：
  - 維持目前的 LOH penalty 但調整權重
  - 或完全移除 LOH 相關的特徵加權

**輸出**：ISM LOH-aware 模式的設計需求文件

#### Q3: cnLOH 是否值得作為獨立生物學類別？

**決策邏輯**：
- 如果 H2 CONFIRMED（cnLOH 保留異質性）：
  - cnLOH 應該在 ISM 中標記為獨立類別
  - cnLOH 區域不應套用 LOH penalty（目前 QS AUC=0.497 的根因之一）
  - 後續研究可探索 cnLOH 內的 epigenetic drift 作為 subclonal marker
- 如果 H2 REJECTED：
  - cnLOH 與 deletion-LOH 可統一處理
  - 簡化分析流程

**輸出**：cnLOH 獨立分類的決策文件

#### Q4: TO 模式是否需要完全不同的 pipeline？

**決策邏輯**：
- 如果 H3 CONFIRMED（self-phasing 產生假關聯）：
  - TO 模式中所有依賴 HP 標籤的特徵（HPMergedDelta, HPSig 等）都不可信
  - 選項 A：TO 使用 PON-only phasing（已驗證 Jaccard=1.0，但需完整重跑）
  - 選項 B：TO 完全不使用 HP-based 特徵，改用 allele-agnostic 甲基化特徵
  - 選項 C：TO 作為獨立 pipeline，使用不同的特徵集和模型
- 如果 H3 REJECTED：
  - TO 可繼續使用現有 pipeline，但需調整 QS 權重

**輸出**：TO pipeline 架構決策文件

---

### 與五個研究目標的對齊

| 研究目標 | 相關假說 | 本步驟的貢獻 |
|---------|---------|-------------|
| 目標 1: per-CpG 評分模型 | H1, H2 | 確認 LOH 類型是否影響 CpG 甲基化的可解釋性；多維標籤（LOH/cnLOH/non-LOH）對評分模型的影響 |
| 目標 2: clone 結構推斷 | H1, H2 | LOH 內外的甲基化異質性差異決定了 clone 結構推斷的適用範圍；cnLOH 內的 epigenetic drift 可能是 subclonal marker |
| 目標 3: 二次打擊偵測 | H4 | 直接驗證 LOH + methylation silencing 的二次打擊模型，決定此功能的可行性 |
| 目標 4: TO normal 補強 | H3 | 確認 self-phasing 的影響範圍，決定 TO 模式需要何種修正（PON-only phasing 或 HP-free features） |
| 目標 5: F1 提升 | H1-H4 | 所有假說的驗證結果都會影響特徵工程和模型設計，間接影響 F1；特別是 LOH penalty 和 TO QS 的修正 |

---

## 輸出

### 主要產出文件

根據驗證結果的成熟度，輸出至不同路徑：

| 條件 | 輸出路徑 | 格式 |
|------|---------|------|
| 4 個假說中 >= 3 個有明確 verdict（CONFIRMED 或 REJECTED） | `docs/reports/research_landscape/07_LOH生物學解釋框架.md` | 正式報告，納入研究推論鏈 |
| 4 個假說中 >= 2 個為 INCONCLUSIVE | `docs/experiments/in_progress/2026/04/` 中的實驗報告 | 階段性報告，標註待解決項目 |

### 附屬產出

- 圖表目錄：對應 `figures/` 子目錄
- 數據表：對應 `data/` 子目錄
- 方向決策文件：納入 `docs/CURRENT_FOCUS.md` 更新

---

## 驗證清單

- [ ] **文獻可查證**：所有引用的文獻（Jenkinson 2020, Knudson 1971, O'Keefe 2010, Jones & Baylin 2002, Esteller 2001 等）均為已發表論文，可在 PubMed 查證
- [ ] **每個假說有明確判定**：H1-H4 各有 CONFIRMED / REJECTED / INCONCLUSIVE 的 verdict，且附有判定依據
- [ ] **方向確認與現有文件一致**：Q1-Q4 的決策結果與 `CURRENT_FOCUS.md` 和 `INDEX.md` 的記錄不矛盾；如有矛盾需更新對應文件
- [ ] **不超出數據支持範圍**：
  - H4 的完整驗證明確標註為需延後至 Phase 2A
  - cnLOH 的統計檢力不足時標記為 INCONCLUSIVE 而非強行判定
  - 所有 confound 分析（特別是 AF、n_reads、self-phasing）已被控制或至少被討論
- [ ] **Paired / TO 分離分析**：所有驗證步驟都分別報告 Paired 和 TO 的結果，不合併分析（遵循系統性觀察的教訓）
- [ ] **與 O1-O13 教訓一致**：
  - 使用 L3 AF-bin 交叉驗證（避免 L2 collider bias）
  - 控制 n_reads 混淆（避免 epipolymorphism 的覆轍）
  - 報告 effect size 而非僅報告 p-value

---

## 時程與優先序

| 子任務 | 優先級 | 預估工作量 | 前置依賴 |
|--------|-------|-----------|---------|
| 4.1 文獻假說框架 | P2 | 0.5 天 | 無（可獨立於 Steps 1-3 先行） |
| 4.2.1 H1 驗證 | P3 | 1 天 | Step 1.2, Step 2.1 |
| 4.2.2 H2 驗證 | P3 | 1 天 | Step 3.2 |
| 4.2.3 H3 驗證 | P3 | 0.5 天 | Step 2.2 |
| 4.2.4 H4 驗證（初步） | P3 | 0.5 天 | Step 1.2 |
| 4.3 方向確認 | P3 | 0.5 天 | 4.2 全部完成 |

**總計**：約 4 天（依 Steps 1-3 完成進度而定）

**關鍵路徑**：4.1 → (Steps 1-3 完成) → 4.2 → 4.3
