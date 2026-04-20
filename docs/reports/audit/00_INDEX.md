<!--
建立時間: 2026-04-20
目標: 22 結論完整 methodology & conclusion audit — 每結論單卡片 + 跨卡分析
處理範圍: 2025-01 ~ 2026-04-20 全部研究方向
用戶決策:
  - CovM hardcoded bug 時機: 立即修 + 重跑
  - 論文定位: 主軸 read-level epigenetic characterization，Subclone Marker Tool 為功能之一
  - 產出物形式: 每結論單卡片（22 × 1）
關聯計劃: /bip7_disk/liaoyoyo2001/.claude/plans/big8-disk-liaoyoyo2001-knowledge-frolicking-hippo.md
-->

# InterSubMod Methodology & Conclusion Audit — 00_INDEX

> **版本**: v1.0 (2026-04-20)
> **審查涵蓋**: 22 結論 × 6 維度（D1–D6）× 跨卡分析
> **來源計劃**: big8-disk-liaoyoyo2001-knowledge-frolicking-hippo.md

---

## 六維度定義

| 維度 | 名稱 | 檢查範圍 |
|------|------|---------|
| **D1** | 內部一致性 | 跨文件數值、術語、定義是否一致 |
| **D2** | 方法論健康度 | pooled OLS trap / L2 collider / n_reads confound / spatial auto-correlation / bootstrap CI / FDR |
| **D3** | 證據鏈完整性 | 依賴結論、被依賴結論、鏈完整度 |
| **D4** | 數據信任度 | dataset 版本、CovM bug 影響、重跑必要性 |
| **D5** | 統計嚴謹度 | effect size、CI、power、multiple testing |
| **D6** | 知識庫交叉驗證 | `/big8_disk/liaoyoyo2001/knowledge/` 支持/挑戰/缺口 |

**評分記號**：✅ 通過 / ⚠️ 有瑕疵但可解釋 / ❌ 明確問題

---

## 22 結論總覽

| # | 簡稱 | 狀態 | 穩定度 | 核心發現 | 所屬樣本 | 卡片 |
|---|------|------|--------|---------|---------|------|
| 1 | Paired/TO 必須分離 | CONFIRMED | ⭐5 | FP rate 1.04% vs 30.6%, HP r=0.001 | 7/7 | [C01](cards/C01_Paired_TO_Separation.md) |
| 2 | PON 移除 99.48% raw FP | CONFIRMED | ⭐3 | 2,717,339 / 2,731,541 | HCC1395-only | [C02](cards/C02_PON_Coverage.md) |
| 3 | TO 無特徵 AUC>0.58 | NEGATIVE | ⭐4 | 60+ 特徵 × 7 樣本天花板 | 7/7 | [C03](cards/C03_TO_AUC_Ceiling.md) |
| 4 | O11 heterogeneity n_reads confound | NEGATIVE | ⭐5 | 原 0.845 → matched 0.560 | 7/7 | [C04](cards/C04_O11_Heterogeneity.md) |
| 5 | O12 LOH 場景不可區分 | NEGATIVE | ⭐4 | AlleleDelta=AF confound, L2 collider | 7/7 | [C05](cards/C05_O12_LOH_Scenarios.md) |
| 6 | O13 跨區域 correlation | NEGATIVE | ⭐5 | shared reads confound | 7/7 | [C06](cards/C06_O13_Cross_Region.md) |
| 7 | G1-G7 germline FP 識別 | NO-GO | ⭐4 | LOSO AUC=0.638, FP removal=0% | 7/7 | [C07](cards/C07_G1_G7_Germline.md) |
| 8 | Read-level FP 識別 | CONDITIONAL NO-GO | ⭐3 | LOSO 0.721, TP loss ≤2% = FP removal 0% | 7/7 | [C08](cards/C08_Read_Level_FP.md) |
| 9 | Self-phasing 是 TO HP_Ratio LOH 主因 | CONFIRMED | ⭐4 | 62% + 17.3:1 bias + r≈0 | 7/7 | [C09](cards/C09_Self_Phasing.md) |
| 10 | PON-only phasing | PARTIAL | ⭐4 | LOH.bed Jaccard=1.0, N50 翻倍 | HCC1395-only | [C10](cards/C10_PON_Only_Phasing.md) |
| 11 | Phase 1A paired F1+0.011 | POSITIVE(弱) | ⭐3 | External +0.0112, Discovery CI 含零 | 7/7 | [C11](cards/C11_Phase1A_F1.md) |
| 12 | ASM 32-66% | POSITIVE | ⭐4 | 5 方法交叉，Jaccard 0.78-0.83 | 7/7 | [C12](cards/C12_ASM.md) |
| 13 | LOH.bed Jaccard=0.847 vs SEQC2 | 觀察值 | ⭐3 | 15.3% 不重疊 | HCC1395-only | [C13](cards/C13_LOH_Bed_SEQC2.md) |
| 14 | QS TO AUC=0.497 | NEGATIVE | ⭐4 | LOH penalty 反向 | 7/7 | [C14](cards/C14_QS_TO.md) |
| 15 | LOH 區域內甲基化全面失效 | ⭐3 | CramersV AUC~0.53 | 7/7 | [C15](cards/C15_LOH_Methylation_Failure.md) |
| 16 | HPFineNGroups somatic marker | POSITIVE REFINED | ⭐4 | NG=4+AF<0.4+NR≥80 TP rate 0.9281 | 7/7(5 ≥0.85) | [C16](cards/C16_HPFineNGroups.md) |
| 17 | LOH Subclone AF×Methylation | POSITIVE | ⭐4 | Δ_NG +0.705 (TO), 7/7 方向一致 | 7/7 | [C17](cards/C17_LOH_Subclone_AF.md) |
| 18 | HCC1954 CNV-driven reversal | CONFIRMED | ⭐4 | HER2+/chr5/8/17 雙路徑 | HCC1954-specific | [C18](cards/C18_HCC1954_Reversal.md) |
| 19 | Z3 amplicon blacklist | CONDITIONAL-NEG | ⭐5 | HCC1954 +0.0065 others -0.0044 | 7/7 | [C19](cards/C19_Z3_Blacklist.md) |
| 20 | Coverage_Multiple 非獨立 CN proxy | FALSIFIED | ⭐4 | z-norm 後 z_extreme 0.15% | 7/7 | [C20](cards/C20_CovM_Non_Independent.md) |
| 21 | LOH.bed 不受 self-phasing 汙染 | CONFIRMED | ⭐5 | Jaccard=1.0000 + C++ audit | 7/7 | [C21](cards/C21_LOH_Bed_Clean.md) |
| 22 | Zone-Aware F1 (H2) vs Characterization | NEGATIVE/CONFIRMED | ⭐4/⭐3 | Z1-Z5 characterization only | 7/7 | [C22](cards/C22_Zone_Aware.md) |

> 結論 22 實為「Zone-Aware F1 (H2) NEGATIVE + Zone Characterization CONFIRMED」雙重，原 00_INDEX 結論 20/21 合併呈現。

---

## 跨卡分析（cross_cutting/）

| 主題 | 檔案 | 涵蓋結論 |
|------|------|---------|
| CovM hardcoded 75.0 bug 影響矩陣 | [CovM_Bug_Impact.md](cross_cutting/CovM_Bug_Impact.md) | 17, 18, 19, 20, 21, 22（Z1-Z5）, 16 B.2-2 |
| Pooled OLS residualization trap 審查 | [Pooled_OLS_Audit.md](cross_cutting/Pooled_OLS_Audit.md) | 4, 5, 6, 15, 16, 17 |
| 全 AUC 掃描 FDR 校正框架 | [Multiple_Testing_Correction.md](cross_cutting/Multiple_Testing_Correction.md) | 3, 4, 5, 6, 7, 11, 16 |
| 主軸定位下的功能清單 | [Characterization_Functions.md](cross_cutting/Characterization_Functions.md) | 12, 16, 17, 18, 21, 22 |

---

## 跨卡 D1–D6 健康度統計

| 結論 | D1 內部一致性 | D2 方法論 | D3 證據鏈 | D4 數據信任度 | D5 統計嚴謹度 | D6 知識庫 | 整體 |
|------|--------------|----------|----------|--------------|--------------|----------|------|
| C01 Paired/TO | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | 🟢 ⭐5 |
| C02 PON | ✅ | ⚠️ | ✅ | ⚠️ | ⚠️ | ✅ | 🟡 ⭐3 |
| C03 TO AUC ceiling | ✅ | ⚠️ | ✅ | ⚠️ | ❌ FDR | ✅ | 🟡 ⭐4 |
| C04 O11 | ✅ | ✅ | ✅ | ✅ | ⚠️ | ✅ | 🟢 ⭐5 |
| C05 O12 | ✅ | ✅ | ✅ | ✅ | ⚠️ | ✅ | 🟢 ⭐4 |
| C06 O13 | ✅ | ✅ | ✅ | ✅ | ⚠️ | ✅ | 🟢 ⭐5 |
| C07 G1-G7 | ✅ | ⚠️ | ✅ | ✅ | ⚠️ Bonferroni | ⚠️ | 🟡 ⭐4 |
| C08 Read-level FP | ✅ | ⚠️ | ✅ | ⚠️ | ⚠️ | ⚠️ | 🟡 ⭐3 |
| C09 Self-phasing | ⚠️ LOH 層次 | ✅ | ✅ | ✅ | ✅ | ✅ | 🟢 ⭐4 |
| C10 PON-only | ✅ | ✅ | ✅ | ✅ | ⚠️ 單樣本 | ✅ | 🟢 ⭐4 |
| C11 Phase 1A F1 | ✅ | ✅ | ✅ | ✅ | ⚠️ effect size | ✅ | 🟡 ⭐3 |
| C12 ASM | ✅ | ✅ | ✅ | ✅ | ✅ | ⚠️ germline/somatic | 🟢 ⭐4 |
| C13 LOH.bed SEQC2 | ❌ r=0.997 誤標 | ✅ | ✅ | ✅ | ⚠️ 單樣本 | ✅ | 🟡 ⭐3 |
| C14 QS TO | ✅ | ✅ | ✅ | ✅ | ⚠️ | ✅ | 🟢 ⭐4 |
| C15 LOH methylation | ⚠️ | ⚠️ pooled OLS | ✅ | ✅ | ⚠️ | ✅ | 🟡 ⭐3 |
| C16 HPFineNGroups | ⚠️ filter 版本 | ⚠️ pooled OLS | ✅ | ⚠️ B.2-2 | ❌ FDR+CI | ⚠️ | 🟡 ⭐4 |
| C17 LOH Subclone AF | ❌ r=0.997 誤標 | ⚠️ pooled OLS | ✅ | 🔴 CovM bug | ⚠️ | ⚠️ | 🔴 ⭐4 |
| C18 HCC1954 reversal | ✅ | ⚠️ 單樣本泛化 | ⚠️ | 🟡 CovM bug | ⚠️ 單樣本 | ✅ | 🟡 ⭐4 |
| C19 Z3 blacklist | ✅ | ✅ | ✅ | 🟡 CovM bug | ⚠️ | ✅ | 🟢 ⭐5 |
| C20 CovM non-indep | ❌ r=0.997 誤標 | ⚠️ | ⚠️ | 🔴 Bug 直接 | ❌ CI | ✅ | 🔴 ⭐4 |
| C21 LOH.bed clean | ✅ | ✅ | ✅ | ✅ | N/A | ✅ | 🟢 ⭐5 |
| C22 Zone-Aware | ❌ r=0.997 誤標 | ⚠️ | ⚠️ Bug 影響 | 🔴 CovM bug | ❌ FDR+CI | ✅ | 🔴 ⭐4/⭐3 |

**圖例**：🟢 結論穩固 / 🟡 需補強但方向穩 / 🔴 bug 或陷阱影響下待驗證 / ✅ 通過 / ⚠️ 有瑕疵 / ❌ 明確問題

### 按維度聚合

| 維度 | ✅ | ⚠️ | ❌ | 主要問題 |
|------|----|----|----|---------|
| D1 內部一致性 | 17 | 3 | 4 | r=0.997 誤標跨 4 卡（C13/C17/C20/C22）；C09 LOH 層次；C15/C16 版本 |
| D2 方法論健康度 | 13 | 8 | 0 | pooled OLS 潛在影響 3 卡（C15/C16/C17）；C07 FDR；C18 單樣本 |
| D3 證據鏈 | 19 | 3 | 0 | CovM bug 影響 C18/C20/C22 鏈完整度 |
| D4 數據信任度 | 15 | 3 | 4 | CovM bug 影響 C17/C20/C22 🔴；C18/C19 🟡 |
| D5 統計嚴謹度 | 4 | 14 | 4 | FDR/CI 系統性缺失（C03/C16/C20/C22） |
| D6 知識庫對照 | 17 | 5 | 0 | germline/somatic（C12）、HER2+ cohort（C18）、臨床對接（C08） |

### 整體分佈

- **🟢 結論穩固**：9 / 22 (40.9%) — C01, C04, C05, C06, C09, C10, C12, C14, C19, C21
- **🟡 需補強**：10 / 22 (45.5%) — C02, C03, C07, C08, C11, C13, C15, C16, C18
- **🔴 bug 影響待驗證**：3 / 22 (13.6%) — C17, C20, C22

### 緊急修正矩陣（CovM bug 修正後必走）

| Card | Bug 等級 | 必修動作 | 影響其他卡 |
|------|---------|---------|-----------|
| C20 🔴 | Bug 直接 | 重算 r / z_extreme | C17 依賴 |
| C22 🔴 | Zone 邊界 | Z1-Z5 重定義 | C19 / C18 依賴 |
| C17 🔴 | step3 CN1 | CovM 重算 + within-group OLS | F2 功能成熟度 |
| C18 🟡 | z-score 重算 | per-sample 特徵穩定性驗證 | - |
| C19 🟡 | Zone 邊界 | Z3 blacklist pilot 重跑 | - |
| C16 🟡 | B.2-2 間接 | 關聯強度重估 | - |

---

## P0–P3 修正動作 Checklist（彙整自 22 卡）

### P0（結論反轉風險 / 必先解決）

- [ ] **P0-A** CovM hardcoded 75.0 infra bug ｜ 影響: C17, C18, C19, C20, C22 ｜ 動作: 另立 `/cpp-change` 修 KDE baseline + 重跑
- [ ] **P0-B** pooled OLS 在 C15/C16/C17 可能保留分組信號 ｜ 動作: within-group OLS 重算 Δ_NG
- [ ] **P0-C** P0-2 haplotag ReadParser 修正未完成 → C03/C05/C08 汙染 ｜ 動作: 依 `/cpp-change` 6 步流程修正

### P1（數值錯誤 / 文獻錯標）

- [ ] **P1-A** `08_Zone_Aware.md:49,73` 誤標 r=0.997（實 r=0.831）｜ 對應 C17/C20
- [ ] **P1-B** 全 AUC 掃描缺 FDR 校正 ｜ 影響 C03/C04/C05/C06/C07/C11/C16
- [ ] **P1-C** HPFineNGroups residualized AUC 缺 bootstrap CI ｜ 對應 C16

### P2（用詞精確化 / 方法論補強）

- [ ] **P2-A** C09 「62% TO LOH 消失」LOH.bed vs HP_Ratio 層次混用 ｜ CURRENT_FOCUS.md:53 / 07_LOH_CN_AF.md:102 補註
- [ ] **P2-B** AUC 0.58 門檻定義跨結論不一 ｜ 對應所有 AUC 聲明結論
- [ ] **P2-C** C12 ASM 32-66% 未分 germline vs somatic ｜ 引用文獻 + 註明類型
- [ ] **P2-D** C16 HPFineNGroups NG=3 非單調現象無解釋 ｜ 補充「生物假說 vs artifact」
- [ ] **P2-E** C09 Self-phasing 演算法層面根因缺 ｜ 補 LongPhase-TO 機制層分析

### P3（文件同步）

- [ ] **P3-A** HPFineNGroups filter 版本：CURRENT_FOCUS.md:98-103 更新為 NG=4+AF<0.4+NR≥80
- [ ] **P3-B** 22 結論總表有 18-22 未納入 06_結論穩定性審查.md（實際：17-21 已納入，19-21 為新增）→ 需補齊 Zone-Aware F1 拆分
- [ ] **P3-C** P0 TODO 編號：CURRENT_FOCUS vs 06_審查 衝突 → 統一
- [ ] **P3-D** 早期 archive TP=30,490/FP=4,842 舊數據殘留 → grep 全檔修正
- [ ] **P3-E** bip8_output_archive PRE-FIX 結論版本警告

---

## 待用戶決策（R-id）

| R-id | 決策點 | 候選 | 建議 | 狀態 |
|------|--------|------|------|------|
| R-01 | CovM bug 時機 | (a) 立即修+重跑 / (b) Phase 2 後 / (c) 僅標註 | (a) | **✅ 已決定: (a)** |
| R-02 | 論文定位 | (a) characterization / (b) subclone marker / (c) self-phasing diag / (d) negative paper | 主軸 (a) + Subclone Marker 為功能 | **✅ 已決定** |
| R-03 | 期刊目標 | (a) Genome Research / (b) Genome Biology / (c) Bioinformatics | (a) | **待決定** |
| R-04 | 13 NO-GO 結論重審 | (a) 接受 / (b) 統一 n_reads confound 檢查 | (b) | **待決定** |
| R-05 | Pooled OLS 全面回溯 | (a) 僅 C15/C17 / (b) 全 residualized 結論 | (b) | **待決定** |
| R-06 | 知識庫交叉驗證範圍 | (a) 僅現有 / (b) 全搜 / (c) 5 個主結論深搜 | (c) | **待決定** |

---

## 開放問題（Q-id）

**戰略層**
- **Q-01** E→A+D→B→C 優先序在 CovM bug 揭露後是否仍成立？
- **Q-02** 若 CovM 修正後 Zone-Aware Z1-Z5 失效，是否重新定義 zone？
- **Q-03** ISM 作為 read-level epigenetic characterization 的具體下游應用？
- **Q-04** TO 模式甲基化增益為負，是否接受為獨立 negative result 章節？
- **Q-05** Phase 2 A+D Normal BAM 成功標準是什麼？

**方法論層**
- **Q-06** 13 個 NO-GO 結論是否值得 ~2 週統一 within-group validation？
- **Q-07** HCC1395 單樣本結論（C02, C10, C13）是否需其他 5 樣本驗證？
- **Q-08** ASM germline vs somatic 區分是否需要 Normal BAM（回到 Phase 2 A+D）？

**生物學層**
- **Q-09** C16 HPFineNGroups NG=3 非單調：生物機制還是 artifact？
- **Q-10** C17 LOH Subclone AF×Methylation 臨床 downstream 可對接哪個 cohort？

---

## 驗收狀態

| 項目 | 狀態 |
|------|------|
| 22 張卡片齊備 | ✅ C01-C22 |
| 每卡欄位完整 | ✅ D1-D6 + 修正建議 + 整體評分 |
| 4 份 cross_cutting 齊備 | ✅ CovM_Bug_Impact / Pooled_OLS_Audit / Multiple_Testing_Correction / Characterization_Functions |
| 跨卡 D1-D6 統計回填 | ✅ 含按維度聚合 + 整體分佈 + 緊急修正矩陣 |
| P0/P1/P2/P3 checklist 完整 | ✅ |
| R-id 已標示待決定 | ✅（R-01/R-02 已決定；R-03/R-04/R-05/R-06 待定） |
| Q-id 開放問題完整 | ✅ Q-01 ~ Q-10 |

---

## 下一步派發（本審查結束後）

1. **立即啟動** `/cpp-change`：CovM KDE baseline 修正（對接 R-01 + CovM_Bug_Impact.md Phase A）
2. **同步啟動** P0-B within-group OLS 重算（C15/C16/C17，對接 Pooled_OLS_Audit.md）
3. **P1-A 立即修**：`08_Zone_Aware.md:49,73` r=0.997 誤標（5 分鐘操作，不需重跑）
4. **P3-A/B/C/D/E** 文件同步：獨立小任務，可平行進行
5. **R-03/R-04/R-05/R-06** 留待用戶最終決策
