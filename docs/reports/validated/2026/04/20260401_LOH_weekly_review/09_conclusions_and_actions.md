# 09. 整合結論與下週行動

<!--
建立時間: 2026-04-01 00:00
目標: 彙整本週 (2026-03-25 ~ 2026-03-31) 所有觀察結論，分類為確認/否決/待定，並提出分級行動建議
處理範圍: Section 01-08 全部結論的整合分析
關聯檔案:
  - docs/experiments/INDEX.md
  - docs/CURRENT_FOCUS.md
  - docs/reports/validated/2026/04/20260401_LOH_weekly_review/01-08 各 section
-->

---

## 1. 本週關鍵結論清單

共 19 項結論，分為三類：

| 類別 | 數量 | 代表性結論 |
|------|------|-----------|
| 確認的結論 | 14 項 | 有充分數據支持，可直接驅動後續行動 |
| 否決的假說 | 4 項 | 明確 negative results，避免重複投入 |
| 待定的問題 | 1 項 | 需要更多數據才能判定（原 3 項中 2 項已解決） |

---

## 2. 確認的結論

| # | 結論 | 數據摘要 | 來源 |
|---|------|---------|------|
| C1 | HP integer tag bug 修正完成，TO regions 品質大幅提升 | 88% TO regions 在 Tier A/A+（修正前僅 ~50%） | HP fix report |
| C2 | TO LOH enrichment = TP-enriched，方向與 paired 翻轉 | TO LOH enrichment ratio = 0.805x（TP 方向），paired 為 FP-enriched | O3 |
| C3 | TO 系統性過判 LOH | 同位點 LOH 判定率：TO 為 paired 的 16-52 倍 | O3, 文獻 |
| C4 | TO 無單一有效特徵 | 所有 TO 特徵 AUC < 0.58；5/9 甲基化方向在 Paired/TO 之間反轉 | O1, O4, O5 |
| C5 | LOH penalty 是 TO QS 失效根因 | LOH penalty trigger: TP 44.5% vs FP 35.8%（反向）；移除後 AUC +0.049 | O2 |
| C6 | Paired/TO 必須分離模型 | 特徵方向反轉 5/9；QS 組件效果相反；LOH enrichment 方向相反 | O1, O2, O3, O5 |
| C7 | GQ 是 paired 最強單一特徵 | AUC = 0.811（paired），TO 無 GQ 可用 | O4 |
| C8 | 樣本異質性極大 | 7 個樣本之間最大 8.6 倍 metric 差異 | O8 |
| C9 | VerificationClass 無效 | Cramer's V = 0.023（閾值 < 0.1 視為無效） | O6 |
| C10 | QS Mode-Aware 修改已完成 | commit b9eaba7；TO: loh_penalty=0, verify_bonus=0；Paired 不變 | Section 07 |
| C11 | COLO829 outlier 根因 = 低測序深度 (~30x) | median num_reads=29（其他 65-103）；hp/total=0.862（最差）；enrichment=0.956x（無區分力）；Z-score=22.09 | Q1 深度分析 |
| C12 | Tier A vs A+ enrichment 差異 = Simpson's Paradox | 6/7 樣本 A+ 優於 A；pool 後反轉因 hp_ratio 離散化放大；甜蜜點 hp_reads 45-50 | Q2 深度分析 |
| C13 | TO LOH TP 富集機制 = TO phasing 系統性 HP 偏移 | TP LOH rate +51.6%（29.3%→44.5%），FP 不變；39,724 位點翻轉中 71.6% TO min(HP1,HP2)=0；7/7 樣本一致 | Q3 深度分析 |
| C14 | O10 methyl_low_fraction AUC 有效（膨脹極小） | Read-level 0.737 → Region-level 0.728 [0.686,0.771]；膨脹僅 0.009；ICC=0.805 但 AUC 未膨脹；跨樣本 0.600-0.840 | O10 bootstrap |

---

## 3. 否決的假說

| # | 假說 | 否決依據 | 教訓 |
|---|------|---------|------|
| N1 | **O11: Within-group heterogeneity 可區分 TP/FP** | epipolymorphism AUC 0.845 -> 0.530 after n_reads correction；高 AUC 源自 n_reads confound 而非生物學信號 | 任何與 read count 正相關的 metric 必須做 n_reads correction |
| N2 | **O12: TO LOH 三場景甲基化可區分** | AlleleDelta = AF confound；L2 collider bias（近常數特徵 residualize on AF 產生虛假信號）；所有 corrected AUC < 0.58 | 必須以 L3 AF-bin 交叉驗證任何看起來有效的 L2 特徵 |
| N3 | **HP0 filter 可改善 TO performance** | HP0 reads 在 TP 中的比例反而高於 FP；filter 方向相反 | TO phasing 品質差使 HP0 比例無法作為可靠指標 |
| N4 | **LOH 作為 binary FP filter** | 直接移除 LOH regions 的 F1 delta < 0（paired 損失 TP 過多；TO 方向反轉） | LOH 資訊應作為 continuous feature 而非 binary filter |

---

## 4. 待定的問題

| # | 問題 | 所需數據 | 預估影響 | 狀態 |
|---|------|---------|---------|------|
| P1 | **O9 FN 觀察**：FN 位點的 ISM 甲基化特徵分佈如何？是否有可 rescue 的子集？ | 需要 FN ISM output（目前僅有 TP/FP） | 若 FN rescue 可行，可能比 FP filter 更有效提升 F1 | 待數據 |
| ~~P2~~ | ~~**HCC1395 chr8 LOH+HPSig 特異性**~~ | ~~chromosome-stratified 深度分析~~ | ~~若為樣本特異性，不宜作為通用 feature~~ | 維持 P2，非阻塞 |
| ~~P3~~ | ~~**Low AF TP LOH rescue**~~ | ~~FN-level 數據 + AF 分層分析~~ | ~~潛力大但不確定~~ | 併入 P1 |

### 本輪新解決的原待定問題

| 原問題 | 結論 | 升級為 |
|--------|------|--------|
| COLO829 phasing coverage 根因 | 低測序深度 (~30x) 連鎖效應，非 melanoma 生物特異性 | C11 |
| Tier A vs A+ enrichment 差異 | Simpson's Paradox + hp_ratio 離散化放大；6/7 樣本 A+ 實際更優 | C12 |
| TO LOH TP 富集機制 | TO phasing 缺乏 normal reference → 系統性單 HP 偏移 → TP LOH rate +51.6% | C13 |
| O10 AUC 膨脹驗證 | Region-level bootstrap AUC=0.728 [0.686,0.771]，膨脹僅 0.009，確認有效 | C14 |
| LOH-like = clonal fraction proxy | AF=0.580 反映 somatic clonality 而非基因組 LOH | 併入 C13 |
| COLO829 concordance 25.1x 受低 eff_hp 放大 | 低深度 hp_ratio 離散化效應 | 併入 C11 |
| TO QS 不存在 paired-like 樣本受損 | TO 全域 AUC=0.497，LOH penalty 反向作用，移除只能改善 | 併入 C10 |

---

## 5. P0/P1/P2 行動建議

### 5.1 行動清單

| 優先級 | 行動 | 依據 | 預期效果 | 狀態 |
|--------|------|------|---------|------|
| **P0** | 移除 TO LOH penalty 和 verify bonus | O2, O3, Section 07 | QS AUC +0.049 (0.497->0.546) | **已完成** (b9eaba7) |
| **P0** | 建立 Paired/TO 分離策略框架 | O1, O5, O7 — 特徵方向反轉 5/9 | 避免 Paired/TO 互相抵消 | 待執行 |
| **P1** | Phase 1A ML 特徵集確定（Paired） | O4 (GQ AUC=0.811), O5 (AF features), O3 (LOH subtype) | Paired AUC ~0.85+ | 待設計 |
| **P1** | TO 專用特徵探索 | O4 (TO 全 AUC<0.58), Direction E (ROCIT 啟發) | 從 region-level 轉向 read-level 特徵 | 待設計 |
| **P2** | 執行 O9 FN 觀察 | 待 FN ISM 數據生成 | 量化 FN rescue 潛力 | 待數據 |
| **P2** | HCC1395 chr8 深度分析 | O3 LOH+HPSig 7.4x enrichment | 判定是否為通用或特異性現象 | 待資源 |

### 5.2 行動依賴關係

```
P0 [已完成] 移除 TO LOH penalty
     |
     v
P0 建立 Paired/TO 分離策略
     |
     +---> P1 Paired 特徵集 (GQ + AF + LOH subtype + effect size)
     |
     +---> P1 TO 特徵探索 (read-level features, Direction E)
              |
              v
         P2 O9 FN 觀察 (需 FN ISM 數據)
         P2 chr8 深度分析
```

---

## 6. 下週預計工作

### 6.1 Week 2026-04-01 ~ 2026-04-06 目標

| 日 | 工作項目 | 預期產出 |
|----|---------|---------|
| W1-W2 | Paired/TO 分離策略文件撰寫 | 策略文件：features 分別列表、model 架構選擇 |
| W2-W3 | Paired ML 特徵集設計 | Feature matrix 定義（GQ, AF, LOH subtype, coverage, effect size） |
| W3-W4 | FN ISM 數據生成 | 執行 ISM 於 FN 位點，產出 FN methylation metrics |
| W4-W5 | O9 FN 觀察初步分析 | FN vs TP/FP 甲基化特徵分佈比較 |

### 6.2 風險與依賴

| 風險 | 影響 | 緩解策略 |
|------|------|---------|
| FN ISM 數據生成耗時 | 延遲 O9 分析 | 先以 HCC1395 單樣本 pilot |
| Paired ML 過擬合 | 泛化性差 | 預留 COLO829 做 cross-validation |
| TO read-level 特徵設計困難 | 延遲 TO 改善 | 先以 ROCIT-style pilot 驗證可行性 |

---

## 7. 研究限制：Cell Line 資料的外推性

**本研究所有結論均基於 7 個癌症細胞系**（HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437, COLO829），不是臨床患者組織樣本。細胞系具有以下特殊性：

- **Purity ≈ 100%**：臨床樣本 purity 通常 20-80%，低 purity 下特徵表現可能不同
- **無微環境**：缺乏 tumor microenvironment、immune infiltration 等臨床因素
- **有限的 LOH 多樣性**：7 個細胞系的 LOH 模式不代表所有腫瘤類型
- **樣本量限制**：所有統計基於 7 組比較，per-sample 推論（如 COLO829 的低 enrichment）需謹慎

這些結論（特別是 enrichment 幅度、特徵 AUC 閾值）推廣到臨床樣本前，需要在不同 purity 和腫瘤類型上進行校正驗證。

---

## 8. 研究方向定位

根據本週結論和 2026-Q2 研究策略（Direction E -> A+D -> B -> C）：

| 方向 | 本週進展 | 下步 |
|------|---------|------|
| **E: Read-level epigenetic context** | 文獻調查完成（ROCIT, MethylBERT 外部驗證） | P1: TO read-level feature pilot |
| **A: Paired model improvement** | GQ AUC=0.811 確認為 anchor feature | P1: ML feature set design |
| **D: Mode-aware QS** | QS mode-aware 重構完成 (b9eaba7) | P0: Paired/TO 分離策略 |
| **B: FN rescue** | 尚未啟動（待 O9 FN 數據） | P2: FN ISM 數據生成 |
| **C: LOH biology integration** | O12 negative; 文獻調查完成 | 暫緩（三場景不可區分） |

---

## 待驗證問題

1. Paired ML 特徵集（GQ + AF + LOH subtype + effect size）是否能達到 AUC > 0.85？
2. TO read-level 特徵是否比 region-level aggregate 有顯著提升？
3. FN 位點的 ISM 甲基化特徵是否存在可 rescue 的子集？
4. QS mode-aware 修改後的實際 benchmark 結果是否與預期一致（TO AUC ~0.546）？
5. 樣本異質性 8.6x 在 cross-sample ML model 中如何處理？（sample-level normalization vs sample-aware features）
6. ROCIT-style read-level classification 在 ONT InterSubMod 資料上的可行性如何？
7. LOH 判別應採用 hp_reads-dependent 閾值（甜蜜點 45-50）而非固定 tier 切分？（C12 啟示）
8. 跨樣本分析是否應排除 COLO829 或以 hp_reads-matched 控制？（C11 建議）

## 認知門檻補充建議

- **AUC 解讀**：0.5 = 隨機猜測，0.7 = 可參考，0.8+ = 實用水準。TO 全特徵 AUC < 0.58 意味著目前沒有任何可靠的 TO-specific filter。
- **F1 delta**：直接移除 LOH regions 的效果以 F1 score 變化量衡量。Delta < 0 表示移除反而降低整體性能。
- **L2 / L3 驗證層級**：L2 為 regression-based（特徵 vs outcome），L3 為 AF-bin stratified cross-validation。L2 的 significant result 必須通過 L3 確認，否則可能是 collider bias。
- **Collider bias**：當一個近常數的特徵被 residualize on AF，可能產生與 AF 完全相關的虛假殘差，導致看似有效但實為 confound 的特徵。O12 的 AlleleDelta 即為此例。
- **Direction E**：研究策略中的 "Read-level epigenetic context" 方向。核心思想是從 region-level 甲基化 aggregate statistics（如 mean methylation, distance matrix）轉向 per-read 甲基化 pattern classification，受 ROCIT/MethylBERT 啟發。
- **FN rescue vs FP filter**：當前研究的一個關鍵轉折——O2-O12 的系統性觀察顯示 ISM 的 FP filter 能力有限（特別是 TO），因此研究重心可能需要從「過濾 FP」轉向「rescue FN」，即利用 ISM 甲基化信號恢復被 caller 錯誤排除的 TP 位點。
