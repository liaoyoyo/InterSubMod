<!--
建立時間: 2026-04-04 21:00
更新時間: 2026-04-18
目標: InterSubMod 研究全景索引 — 多檔案文件系統入口
處理範圍: 2026-03-11 至 2026-04-18 全部研究方向的整理、驗證與路線圖
關聯檔案:
  - docs/reports/research_landscape/01_TO_FP問題全貌.md
  - docs/reports/research_landscape/02_Self_Phasing根因.md
  - docs/reports/research_landscape/03_ISM分析價值界定.md
  - docs/reports/research_landscape/04_暫停判定與重評估.md
  - docs/reports/research_landscape/05_證據鏈總覽.md
  - docs/reports/research_landscape/06_結論穩定性審查.md
  - docs/reports/research_landscape/07_LOH_CN_AF_研究總整理.md
  - docs/reports/research_landscape/08_Zone_Aware.md
  - docs/reports/research_landscape/09_Part_B.md
-->

# InterSubMod 研究全景文件系統

> **版本**: v1.2 (2026-04-18)
> **涵蓋範圍**: 2026-03-11 ~ 2026-04-18 全部研究方向
> **目的**: 輔助研究者與合作者快速理解目前研究狀態、推論邏輯與待驗證事項

### 關鍵圖表速覽

| 圖表 | 說明 | 位置 |
|------|------|------|
| ![穩定度+HP分類](figures/01_stability_and_hp_dependency.png) | 14 結論穩定度 + 特徵 HP 依賴分類 | 本文件 / 06 |
| ![FP過濾漏斗](figures/02_fp_filtering_funnel.png) | TO FP 三層過濾：PON 99.48% → ISM 0% | 01 |
| ![Self-Phasing影響](figures/03_self_phasing_impact.png) | HP bias 17.3:1 + LOH 組成 + 修正預測 | 02 |
| ![證據鏈依賴](figures/04_evidence_chain_dependencies.png) | 8 條證據鏈層次與依賴關係 | 05 |
| ![決策樹](figures/05_decision_tree.png) | 修正後研究路徑分支 | 04 |
| 🗺️ Zone-Aware Framework | Z1-Z5 zone 定義 + TO Z3 最低 TP rate 0.608 | 08 |
| 🔬 Part B 驗證矩陣 | 8 項質疑全通過 + F pilot 新 canonical filter | 09 |

---

## 一句話摘要

**InterSubMod 嘗試用甲基化特徵區分 TO 模式的 TP/FP，歷經 60+ 特徵 × 3 層次 × 7 樣本的系統探索後，發現 (1) self-phasing 循環依賴汙染了約 38% 的特徵，(2) 不依賴 HP 的 55% 特徵全部 AUC < 0.64，(3) 根因是 TO 的 31% FP 為 germline leak，甲基化上與 TP 不可區分。研究重心已從「variant filter」轉向「read-level epigenetic characterization」。**

---

## 文件導覽圖

```mermaid
graph TD
    INDEX["00_INDEX<br/>📌 本文件 — 全局入口"]

    F01["01_TO_FP問題全貌<br/>📊 TO 30.6% FP 從何而來？<br/>PON 過濾了什麼？剩什麼？"]
    F02["02_Self_Phasing根因<br/>🔬 為什麼 TO 的 HP tag 不可信？<br/>完整因果鏈 + 量化影響"]
    F03["03_ISM分析價值界定<br/>💡 ISM 能做什麼？不能做什麼？<br/>HP依賴 vs 非依賴特徵分類"]
    F04["04_暫停判定與重評估<br/>⏸️ 哪些結論需要修正後重測？<br/>修正後的預期改善幅度"]
    F05["05_證據鏈總覽<br/>🔗 8 條完整證據鏈<br/>前提→設計→結果→依賴"]
    F06["06_結論穩定性審查<br/>⚖️ 14 個結論的穩定度評分<br/>隱含假設 + 驗證需求"]
    F07["07_LOH_CN_AF_研究總整理<br/>🧬 LOH/CN/AF 三維度統合<br/>已關閉方向 + 未探索方向"]
    F08["08_Zone_Aware<br/>🗺️ Zone-Aware Framework<br/>5 zone 定義 + H1-H5 驗證<br/>F1 ❌ / Characterization ✅"]
    F09["09_Part_B<br/>🔬 Part B 質疑驗證<br/>8 項反挑戰 + F pilot<br/>⭐3 → ⭐4 升級"]

    INDEX --> F01
    INDEX --> F02
    INDEX --> F03
    INDEX --> F04
    INDEX --> F05
    INDEX --> F06
    INDEX --> F07
    INDEX --> F08
    INDEX --> F09

    F01 -->|"FP 來源確認後"| F02
    F02 -->|"self-phasing 影響範圍"| F03
    F03 -->|"哪些結論受影響"| F04
    F05 -->|"每條鏈的穩定性"| F06
    F01 & F02 & F03 -->|"LOH/CN/AF 維度統合"| F07
    F07 -->|"LOH/CN/AF 可行動化"| F08
    F07 -->|"HPFineNGroups 挑戰"| F09

    style INDEX fill:#4a90d9,color:#fff
    style F01 fill:#e8f4fd
    style F02 fill:#fde8e8
    style F03 fill:#e8fde8
    style F04 fill:#fdf8e8
    style F05 fill:#f0e8fd
    style F06 fill:#fde8f0
    style F07 fill:#e8eaf6
    style F08 fill:#e1f5fe
    style F09 fill:#fff3e0
```

---

## 建議閱讀順序

| 閱讀目的 | 建議路徑 | 預估時間 |
|----------|----------|----------|
| **快速掌握全貌** | 本文件 → 01 (前半) → 02 (機制圖) | 15 min |
| **理解推論邏輯** | 01 → 02 → 03 → 04 | 40 min |
| **審查結論可靠性** | 05 → 06 | 30 min |
| **規劃下一步實驗** | 04 → 06 (驗證需求欄) | 20 min |
| **LOH/CN/AF 維度統合** | 07（三維度交叉 + 未探索方向） | 25 min |
| **Zone-Aware 框架** | 07 → 08（從 LOH/CN/AF 到 zone confidence） | 20 min |
| **Part B 質疑驗證** | 09（8 項反挑戰 + F pilot 新 filter） | 30 min |

---

## 全局研究狀態速覽

```mermaid
pie title 14 個結論的穩定度分佈
    "穩定度 5（堅固）" : 3
    "穩定度 4（穩固）" : 6
    "穩定度 3（需注意）" : 5
```

### 研究方向結論總表

| # | 方向 | 結論 | 穩定度 | HP 依賴 | 文件 |
|---|------|------|--------|---------|------|
| 1 | Paired/TO 分離 | **CONFIRMED** | ⭐5 | 無 | 01, 05 |
| 2 | PON 移除 99.48% FP | CONFIRMED | ⭐3 | 無 | 01 |
| 3 | TO 無特徵 AUC>0.58 | **NEGATIVE** | ⭐4 | 混合 | 03, 05 |
| 4 | O11 heterogeneity | **NEGATIVE** | ⭐5 | 無 | 05, 06 |
| 5 | O12 LOH 場景 | **NEGATIVE** | ⭐4 | 部分 | 05, 06 |
| 6 | O13 跨區域 correlation | **NEGATIVE** | ⭐5 | 無 | 05, 06 |
| 7 | G1-G7 germline FP 識別 | **NO-GO** | ⭐4 | 無 | 03, 05 |
| 8 | Read-level FP 識別 | **CONDITIONAL NO-GO** | ⭐3 | 部分 | 04, 06 |
| 9 | Self-phasing 是 TO LOH 主因 | **CONFIRMED** | ⭐4 | — | 02, 05 |
| 10 | PON-only phasing | **PARTIAL SUCCESS** | ⭐4 | — | 02, 04 |
| 11 | Phase 1A paired F1+0.011 | **POSITIVE** (弱) | ⭐3 | 部分 | 04, 06 |
| 12 | ASM 32-66% | **POSITIVE** | ⭐4 | 部分 | 03, 05 |
| 13 | LOH.bed Jaccard=0.847 | 觀察值 | ⭐3 | 無 | 06 |
| 14 | QS TO AUC=0.497 | **NEGATIVE** | ⭐4 | 間接 | 03, 06 |
| 15 | LOH Subclone AF×Methylation | **POSITIVE** (雙模式確認) | ⭐4 | 部分 | 07 |
| 16 | CNV zone-aware filter | **CLOSED** | ⭐4 | 無 | 07 |
| 17 | Coverage_Multiple CN proxy | **CONFIRMED** (r=0.831) | ⭐4 | 無 | 07 |
| 18 | 甲基化 vs CN | **NEGATIVE** (CN-blind) | ⭐4 | 無 | 07 |
| 19 | GC bias 校正 | **NOT NEEDED** | ⭐4 | 無 | 07 |
| 20 | Zone-Aware F1 (H2) | **NEGATIVE** | ⭐4 | 無 | 08 |
| 21 | Zone Characterization | **CONFIRMED** | ⭐3 | 部分 | 08 |
| 22 | HPFineNGroups (Part B 精確化) | **POSITIVE REFINED** | ⭐4 | 間接 | 09 |

---

## 核心問題的三層結構

```mermaid
graph TB
    subgraph L1["Layer 1: Caller 層（ISM 無法控制）"]
        A1["ClairS-TO 輸出<br/>TP 11,598 + FP 14,244<br/>FP rate = 30.6%"]
        A2["PON 過濾<br/>移除 99.48% raw FP<br/>（2,717,339 / 2,731,541）"]
        A3["殘餘 FP = 14,202<br/>98.6% = raw_absent<br/>（paired 從未 call 出）"]
    end

    subgraph L2["Layer 2: Phasing 層（已發現根因）"]
        B1["LongPhase-TO<br/>somatic 自己 phase 自己<br/>= 循環依賴"]
        B2["HP tag 汙染<br/>94.6% somatic → HP1<br/>bias 17.3:1"]
        B3["TO LOH 62% 是 artifact<br/>HP_Ratio 跨模式 r=0.001"]
    end

    subgraph L3["Layer 3: ISM 特徵層（系統探索完成）"]
        C1["55% 特徵不依賴 HP<br/>~42 個特徵<br/>全部 AUC < 0.64"]
        C2["38% 特徵直接依賴 HP<br/>~29 個特徵<br/>TO 結果不可信"]
        C3["7% 間接依賴 HP<br/>~14 個特徵<br/>TO mode 已移除大部分影響"]
    end

    A1 --> A2 --> A3
    A3 -->|"FP 是 germline leak"| L3
    B1 --> B2 --> B3
    B3 -->|"HP 汙染"| C2

    style L1 fill:#fff3e0
    style L2 fill:#fce4ec
    style L3 fill:#e8eaf6
```

---

## 研究路線圖

```mermaid
gantt
    title InterSubMod 研究路線圖 (2026 Q1-Q2)
    dateFormat  YYYY-MM-DD

    section Phase 0 基礎設施
    HP tag bug 修正           :done, hp, 2026-03-30, 1d
    Self-phasing 因果鏈驗證   :done, sp, 2026-04-02, 2d
    PON-only phasing 測試     :done, pon, 2026-04-03, 1d
    HP 分類報告(A/B/C/D)      :done, cls, 2026-04-04, 1d

    section Phase 1 核心修正 (P0)
    ISM ReadParser 修正       :active, rp, 2026-04-07, 3d
    HP-dependent 特徵重測     :hp2, after rp, 5d
    within_dom_alt_frac 重評  :wda, after hp2, 3d

    section Phase 1A ML
    Paired methyl+context     :done, p1a, 2026-03-20, 10d
    TO 重測(待 phasing 修正)  :to1a, after hp2, 5d

    section Phase 1B LOH Evidence Panel
    LOH.bed 機制確認          :loh1, 2026-04-10, 3d
    Region-level LOH panel    :loh2, after loh1, 7d

    section Phase 2 Normal Methylation
    Baseline construction     :nm, 2026-05-01, 14d
```

---

## 待驗證事項清單 (按優先級排序)

### P0 — 阻塞性（必須先完成才能推進）

- [ ] **ISM ReadParser 修正**：忽略 somatic HP tags (HP:i:11/21/33)，只使用 germline HP:i:1/2
- [ ] **HP-dependent 特徵全量重測**：PON-only phasing + 修正 ReadParser 後，重跑 29 個 B 類特徵
- [x] **LOH.bed 生成機制確認**：已確認基於 VCF AF/VAF（`PhasingGraph.cpp:1817`）（2026-04-11 完成）

### P1 — 重要（影響結論可靠性）

- [ ] **多樣本 PON-only phasing**：在 COLO829, HCC1954 等 2-3 個額外樣本驗證
- [ ] **within_dom_alt_frac 重測**：self-phasing 修正後 AUC 變化（預期可能改善）
- [ ] **Purity estimation 修復**：PON-only 後 purity 從 0.927→0，需排查根因
- [ ] **H2009 Phase 1A 負向根因**：7 個樣本中唯一負向，需解釋

### P2 — 補強（提升論文說服力）

- [ ] PON 覆蓋率跨樣本穩定性（6 個非 HCC1395 樣本）
- [ ] ASM 效應量門檻分析（|delta|>0.2 後比例變化）
- [ ] GBM/RF 非線性模型測試（確認 LR LOSO 0.638 是天花板）
- [ ] 低純度臨床樣本驗證（purity 30-70%）
