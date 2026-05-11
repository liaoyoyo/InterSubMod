<!--
建立時間: 2026-04-23
對應週報: docs/reports/validated/2026/04/20260423_研究週報_20260416_20260423_NG2_LOH_constrained_phasing與TO_pivot_01.md
受眾: PI（熟悉 ONT、cancer genomics、somatic calling；對 InterSubMod 專案全景已掌握至 0421 週報）
故事弧: Three-Act（定錨+CN KDE方法學 → filter framework卡關+HP雜訊方案 → NG=2 真機制揭露+TO論文主軸 pivot）
工具: Phase 3.5 導演審查；下一步 Phase 4 三件套 + build_pptx
總 slide 數: 26
-->

# PPTX Storyboard — NG=2 LOH-constrained phasing 發現與 TO-層論文主軸 pivot

## 導演指引

### 三幕結構

```
Act I · 定錨與方法學（Slides 1-9, ~30%）
    封面 → 脈絡 roadmap → 上週回顧 → 本週提問 → CN 方法學升級（舊 75× → KDE）→ LOH×AF×CN 切分框架
    目的：建立「測量單位正確」+「切分框架初步 POSITIVE」的信任

Act II · 卡關與方案（Slides 10-17, ~30%）
    S4 ambiguous bucket 問題 → 引出 HP/ISM 特徵擴展需求 → TO self-phasing 雜訊問題 →
    --germline-hp-only 方案（PON 錨點）→ Phase 1 機制 ✓ filter FAIL
    目的：揭示「加特徵並未解決 S4」與「tag 方案機制對但 filter 無增益」的雙重卡關

Act III · 機制揭露與 pivot（Slides 18-24, ~30%）
    用戶提問「NG=2 與甲基有關係嗎」→ C++ 原始碼回查 → 4 bucket 純 phasing →
    obs18 跨 6 樣本驗證 → Inner same-hap 93-99% + TP gap +0.37 → LOH-constrained phasing 機制 →
    TO 論文主軸 pivot
    目的：呈現本週偵探時刻與論文主軸重定位

結語（Slides 25-26, ~10%）
    結論總表變動 → 下週優先行動
```

### 敘事層次分配原則

| 資訊類型 | 放置位置 |
|---------|---------|
| **核心結論數字**（2-3 個 key numbers） | Slide 主視覺 + 標題 |
| **一張關鍵圖**（讓 PI 一眼看出訊號）| Slide 視覺區 55-65% |
| **推論邏輯鏈**（為何這個結論合理） | Speaker notes |
| **方法學細節**（KDE 演算法、Fisher p、CI 計算） | Speaker notes |
| **已排除的替代解釋** | Speaker notes |
| **C++ 原始碼片段** | 僅 Thread D 機制揭露 slide 視覺，其餘 speaker notes |
| **統計檢定細節**（Wilson CI、Fisher odds ratio、bootstrap） | Speaker notes |
| **無 F1-relevant 的中間嘗試**（如 V2/V3 失敗迭代） | 完全**捨棄**（不進 PPT 也不進 notes） |
| **已關閉方向的歷史**（0406 LOH、Beyond-AUC 耗盡） | 僅在對應 slide 一行帶過，不開展 |

### 設計常數（承襲 PPTX_PROTOCOL）

- 色彩：`dark=#1E2A44`（標題/文字）、`bg=#F7F3EC`（背景）、`accent=#A85540`（強調）、`green=#009E73`（POSITIVE/PASS）、`red=#D55E00`（NEGATIVE/FAIL）、`blue=#0072B2`（對照）
- 字體：標題 Arial Bold 32-36pt、小標 Arial Bold 24pt、內文 Arial 16-18pt、caption 12-14pt、下限 10pt
- 版面：視覺 55-65%、留白 20-30%、文字 ≤15%、每 slide ≤4 bullet
- 雙語：中文主 + 英文副（EN 字號 60% + 縮排 0.25"）

---

## Act I · 定錨與方法學（Slides 1-9）

### Slide 1 · 封面

- **核心訊息**：週報定位為「方法學升級 + 機制更正」雙主線
- **視覺**（PPT 上放）：
  - 主標題「NG=2 LOH-constrained phasing 發現與 TO-層論文主軸 pivot」
  - 副標英文 "Mechanism correction: from methylation bimodality to phasing signatures in tumor-only sequencing"
  - 涵蓋區間 2026-04-16 ~ 2026-04-23
  - 報告人、日期
  - 背景圖：一張淡化的 obs18 堆疊圖作為視覺錨點（提示本週高潮）
- **口頭**（speaker notes）：
  - 本週最大發現：用戶 04-22 晚間一句「NG=2 與甲基有關係嗎」觸發 C++ 原始碼回查，推翻 HPFineNGroups 的 methylation bimodality 誤解
  - 週報涵蓋 8 天（與上週重疊 6 天）— 以 CN KDE 方法學為起點往上建，避免重複 0421 週報已展開的內容
- **捨棄**：無

---

### Slide 2 · 研究脈絡 Roadmap

- **核心訊息**：本週 4 Thread 在過去已關閉結論 × 下週優先行動之間的位置
- **視覺**：Mermaid 脈絡圖（與週報 Layer 0.1 同版，三欄：Past / ThisWeek / Future），4 Thread 顏色區分
- **口頭**：
  - 過去三大已關閉：LOH 雙定義 filter fail、Self-phasing 因果鏈、0421 基礎設施三軌
  - 本週 4 Thread 並非獨立 —— A 的新 CN tier 是 B 的測量單位；B 的 S4 bucket 卡關引出 C 的 HP 擴展；C 的 flag=on NG≥3=0 觀察是 D 的獨立 negative control
  - 下週 P0 只有 2 項：paired 對照（驗證 H-D3）+ Archive TO 6 樣本重跑
- **捨棄**：節點內的詳細 figure 數字（都在後續 slide 展開）

---

### Slide 3 · 上週結論回顧（橋接 0421 週報）

- **核心訊息**：進入本週時的三項「已解」結論 + 兩項 blocker
- **視覺**：2×3 格子，上排「已解」下排「未解」
  - 已解：LOH.bed Jaccard=1.0 ✓（獨立於 self-phasing）/ R1-Global Overfit Collapse 0.527 / F1-filter 路線終結（Verdict ΔF1=0）
  - 未解：stale-binary → master dataset 需重跑 / `--germline-hp-only` flag=on NG≥3=0 → marker 根基需重驗
- **口頭**：
  - 0421 週報核心訊息：「ISM 的論文定位從 variant filter → epigenetic characterization」已三重確認
  - 本週兩 blocker 進展：Blocker 1 → Thread A 方法學定案；Blocker 2 → Thread C/D 雙軌回應
- **捨棄**：
  - ClairS-TO Verdict 詳細數字（Germline FP 96.1% 等）
  - R1-Global 的 CI 細節
  - 0421 週報 Thread A/B/C 完整故事（PI 已看過）

---

### Slide 4 · 本週提問（定錨主軸）

- **核心訊息**：從 LOH 內 AF × CN 高 TP rate 切分出發 —— 其他特徵能否進一步切分 TP/FP？
- **視覺**：
  - 主題問句粗體置中：「除了 LOH、AF、CN，**加入 HP / ISM 甲基 / ISM 聚類特徵**是否能進一步區分 TP/FP？」
  - 副標：四個**子問題**小標（不展開，只提示）
    - Q1 · CN 怎麼切才準確？
    - Q2 · LOH × AF × CN 能切到什麼程度？
    - Q3 · HP 雜訊是否在拖累？
    - Q4 · NG 到底是什麼？
- **口頭**：
  - 用戶核心主軸：**「LOH 範圍、ALT 比例、還有這位置的 read 數量相對平均情況」是起點，但這 3 維已無法區分完整 TP/FP 結構**
  - 子問題 Q1 對應 Thread A；Q2 對應 Thread B；Q3 對應 Thread C；Q4 對應 Thread D
  - 四問題**有因果關係**：Q1 先解才能談 Q2；Q2 卡關（S4 bucket）才需 Q3；Q3 失敗的觀察直接成為 Q4 的線索
- **捨棄**：詳細子問題 bullet（下面 slide 各自展開）

---

### Slide 5 · CN 的老問題（Q1 引入）

- **核心訊息**：舊 pipeline 用 `--expected-coverage 75.0` 當所有樣本的 diploid baseline，這不準確
- **視覺**：左右對照
  - 左：「舊做法」— 7 樣本共用 75× hardcoded → CN tier 切分共用 ×0.5/×1.5/×2.5/×3.5
  - 右：「實際情況」— HCC1395 ≈ 54×（SEQC2 neutral median），舊值 bias +39%
  - 一句標題：**共用 75× 讓 CN tier 系統性偏移**
- **口頭**：
  - CN tier 是 Thread B 的測量單位；測量單位不準 → B 的所有 heatmap 比較都是 apples-to-oranges
  - 為何沒早發現：master dataset 是 2026-03-30 產出的 stale binary，至本週 KDE-fixed binary 已 commit 但尚未全量重跑
  - **per-sample KDE 必要性**：7 樣本實際 coverage 分佈差異大（HCC1395 54×、COLO829 ~47×、H2009 ~65×），單一預設物理上不可能涵蓋
- **捨棄**：
  - stale-binary 事件的詳細時間線
  - master dataset 的 40,096 rows / pilot 28,495 rows 差異（數字細節在 notes 說明是 TP-track filter 差異而非 KDE 效應）

---

### Slide 6 · KDE 雙 Pass 方法（Q1 解決方案）

- **核心訊息**：KDE 動態校準 per-sample 2× 高峰，偏差 <2%
- **視覺**：
  - 圖 A1（`fig3_two_pass_architecture.png`）架構圖放大展示
  - 三階段標註：Pass 1 並行 → Mid KDE peak detection → Pass 2 重算 CovM
  - 關鍵數字 inset：HCC1395 pilot **53.0×**（bias −1.9% vs SEQC2 54×）
- **口頭**：
  - Pass 1 並行計算 region metrics，seed=75.0 作為臨時 fallback
  - Mid：histogram → Gaussian smooth → 2nd-deriv peak detection，取 per-sample 2× diploid peak
  - Pass 2 重算 `Coverage_Multiple = NumReads / diploid_coverage_kde`
  - 5 種 fallback 路徑（valid regions<100、histogram 過窄、mode 超 sanity range [10×,300×]、user-specified、KDE 成功）都有 TSV audit column 紀錄
  - 2 個 commits：`374fad4`（audit column + DiploidEstimate struct）、`12d9b3e`（explicit fallback flag 取代 sentinel 誤判）；9+1 unit tests 全通過
- **捨棄**：
  - C++ 程式碼細節（僅 commit SHA 即可）
  - RegionProcessor.cpp 行號等實作細節

---

### Slide 7 · CN 修正量化驗收

- **核心訊息**：CovM median 0.880 → 1.245（×1.415 恰為 75/53），CN category 大幅重分類
- **視覺**：
  - 左半：圖 A2（`fig1_covm_distribution_shift.png`）兩組 CovM 密度對照
  - 右半：圖 A3（`fig2_category_reclassification.png`）Normal −5718 / CNV_Gain +2956 / High_Copy +2710 柱狀
  - 底標：**物理意義**：KDE 修正讓 diploid peak 真的落在 CovM=1，CN gain 不再被錯分為 Normal
- **口頭**：
  - ×1.415 不是巧合 —— 正好是 75/53，證明 CovM 只是被等比例拉伸
  - Normal → CNV_Gain 遷移代表：原本 NumReads ≈ 100 的 region 在 stale binary 下 CovM = 1.33（勉強 CNV_Gain）；KDE 後 CovM = 1.89（清楚 CNV_Gain）
  - 下游影響邊界：**影響**跨樣本 CovM 絕對值比較；**不影響** percentile-based filter（Variant A 數學證明 scale-invariant）、per-sample 內部排序、LOH.bed（Jaccard=1.0）
- **捨棄**：
  - 圖 A5/A6/A7（per-sample quantile drift / category migration / PCA）— 屬於跨樣本延伸，當前 pilot 僅 HCC1395 完成，其他樣本仍待 master rerun；放 speaker notes 提及即可

---

### Slide 8 · LOH × AF × CN 切分 · TP 熱圖（Q2 主結果）

- **核心訊息**：HCC1395 TO baseline TP 71.1% 被拆成**雙極分佈**
- **視覺**：
  - 主視覺：圖 B2（`fig_v2_1_to_tp_heatmap.png`）TP rate heatmap
  - 右側小字標註三個區域：
    - ⬛ 綠區：LOH 任何註解 × Extreme AF（TP 88-96%）
    - ⬛ 綠區：None × Near-half × CN∈T1/T2（TP 93-96%，canonical het）
    - ⬛ 紅區：None × Intermediate × CN≥3（TP 47-60%，FP hotspot）
- **口頭**：
  - 雙極分佈 = biology-informed filter 可行性的直接視覺證據
  - 綠區生物詮釋：LOH 區 Extreme AF = deletion/cnLOH 單 hap survive（FACETS/Battenberg 模型）；None + Near-half + CN~2 = Hardy-Weinberg 預期的 canonical somatic
  - 紅區（CN≥3 + Intermediate AF）有兩種典型 artifact：(a) repeat/segdup mapping chaos；(b) low-purity subclone leak
- **捨棄**：
  - FP 熱圖（B3）— 放下一 slide 或完全不放（與 TP 熱圖互補，視覺重複）
  - per-cell 詳細數字（PI 只需看色塊）

---

### Slide 9 · S1-S7 Filter Scheme 對照 + S3 勝出

- **核心訊息**：S3 Diploid Het **TP 95.5%**、S5 combo **FP reduction 99.37%**、S4 無辨別力
- **視覺**：
  - 主視覺：圖 B4（`fig_v2_3_filter_scheme_bar.png`）S1-S7 bar chart
  - 側邊小表格：S3 / S5 / S4 三列對照（TP rate / TP recall / FP reduction / fold-improvement）
  - 伏筆標註：**S4 含 75% TP + 76% FP → 需二級判別** ➡️（箭頭指向右下，提示 Act II）
- **口頭**：
  - S1 Deletion/cnLOH + Extreme (TP 90.1%, n=292) / S2 Subclonal LOH + Intermediate (87.4%, n=214) / S3 Diploid Het (95.5%, n=380) / S5 combo (91.8%, n=886, FP ↓99.37%)
  - S3 TP:FP fold-improvement 8.69× 是本週最強單一 filter 模組
  - S4 (LOH=None + Extreme AF) 含 75% TP + 76% FP 混雜 —— 不是「FP bucket」而是「無辨別力 bucket」；baseline TP 71.1% 無法用 5 維切分
  - S5 covers only 2.85% TP recall → 高 precision 白名單，**不是替代 caller，是 confidence tag**
- **捨棄**：
  - S6/S7（+NG≥3）— 邊際貢獻 <1pp，圖中次要，notes 一句帶過即可
  - paired_full 7 樣本 per-sample scheme 穩定性（圖 B5）— 僅確認 "框架不崩"，非主訊息；notes 提及

---

## Act II · 卡關與方案（Slides 10-17）

### Slide 10 · S4 卡關 → 特徵擴展需求（Q3 引入）

- **核心訊息**：S4 bucket 的 TP/FP 混雜無法用 LOH/AF/CN 解決，需加入 HP / ISM 特徵
- **視覺**：
  - 圖形化表達：一個大方塊「S4: n=30,432」裡面 TP (67%) 和 FP (33%) 交錯分佈，外圍 5 維座標（LOH×AF×CN×mode×sample）標 ❌
  - 右側箭頭指向三個候選特徵：🧬 HP 狀況 · 🧬 ISM 甲基 · 🧬 ISM 關聯/聚類
- **口頭**：
  - 用戶主軸清晰：「除了 LOH、AF、CN 之外，需要使用其他特徵...加入常用的 HP 狀況、ISM 甲基結果、ISM 關聯與聚類結論」
  - 注意：ISM 甲基 + 聚類在過去 Beyond-AUC 耗盡（2026-04-09）已全 AUC ≤ 0.58；但那是 pooled 測試，**biology-module 內是否仍有邊際增量尚未測**
  - 先查 HP → 因為 HP 特徵在 0418 F pilot 有 AUC ≥0.60 記錄（HPFineNGroups 0.617）
- **捨棄**：ISM 甲基 / 聚類的歷史細節（Beyond-AUC 耗盡的 25 特徵完整表）

---

### Slide 11 · HP 雜訊的根因（Self-Phasing 複習）

- **核心訊息**：TO self-phasing 讓 somatic HP tags 混入 HP1-1/HP2-1 分群，污染所有 HP-derived 特徵
- **視覺**：
  - 簡化因果圖：LongPhase-TO step03 → tumor_tagged BAM → HP1/HP1-1/HP2/HP2-1 都含 somatic 來源 → HPFineNGroups 的 fine-group 結構可能是 artifact
  - 關鍵數字引用：0421 週報「somatic bias 17.3:1 → 消除（PON-only phasing）」、「62% ISM LOH 消失」
- **口頭**：
  - 這是 0421 週報 Thread A 的結論延續：PON-only phasing 在上游已驗證（LOH.bed Jaccard=1.0 不變），但**下游 ReadParser 層仍保留 LongPhase 的 somatic HP tags** → 未解決
  - 為什麼本週要再動：0418 F pilot 的 HPFineN marker AUC 0.617 是在 somatic-tainted HP tag 上測得 → 機制可能是「真 subclone marker」，也可能是「self-phasing artifact」
  - 需要 orthogonal test：降低 somatic HP 參與 → 看 marker 是否消失
- **捨棄**：
  - Self-phasing 因果鏈的完整 62% LOH 消失細節（已在 0406 週報、CURRENT_FOCUS 詳載）
  - PON-only phasing 上游驗證的 N50 / phased rate 等指標

---

### Slide 12 · 方案：`--germline-hp-only` flag（PON 錨點）

- **核心訊息**：用戶提的設計 —— PON 當 germline block 錨點，只保留 PS ∈ germline block 的 HP tag；somatic 來源 demote
- **視覺**：
  - 左：原設計（flag=off）—— HP tag from {germline phasing, somatic phasing} 混合參與 bucket
  - 中：新設計（flag=on）—— 只信任 germline block PS tag；somatic HP demote 回 HP1/HP2/other
  - 右：預期效果 —— 4 bucket 分群純度提升，fine-group 結構應收斂
  - Commit SHA 標註：`775036c` (Phase 0 flag) + `a61779c` (Phase 1 validation docs) + `4dc2d73` + `2e2df22`
- **口頭**：
  - 用戶原話「使用 PON 當 tag 的錨點，somatic 用於分出 HP1-1 與 HP2-1 與 HP3」 — 這是**方法學設計**，實作位於 ReadParser 層（`src/core/ReadParser.cpp:123`）
  - 為何在 ReadParser 而非 LongPhase 上游：降低實作複雜度 + 不破壞既有 pipeline + default=off 不影響既有流程
  - 預期：filter AUC 改善（清乾淨 HP tag → 純度提升）、HP_Ratio 往 0.5 收斂、HPFineN 結構更接近真 biology
- **捨棄**：
  - flag 的 CLI 參數細節
  - Phase 0 smoke test HCC1395 chr19 subset 結果（僅作 sanity check）

---

### Slide 13 · Phase 1 機制驗證 ✓（Audit 獨立性）

- **核心訊息**：flag 機制正確 —— 僅 demote 不 remove，somatic tag 全基因 sum 兩 flag 下 identical
- **視覺**：
  - 主視覺：圖 C4（`fig_c4_audit_independence.png`）Audit 獨立性柱狀圖
  - 三個綠色 ✓ 標記 NHP_Somatic11 / 21 / 33 = identical
- **口頭**：
  - 這是**機制正確性的關鍵證據**：若 flag 有 side effect（誤 remove reads），下游失效時就分不清是「機制錯」還是「機制對但無用」
  - Audit identical → 排除「flag 實作錯」，任何下游觀察可安全歸因於機制效應
- **捨棄**：
  - NHP_Somatic11/21/33 的數值本身（347,213 / 308,762 / 124,096）— PI 只需看到「相同」即可
  - Somatic tag per-site 平均密度細節

---

### Slide 14 · Phase 1 AUC Gate FAIL

- **核心訊息**：18 個 HP-related 特徵無一達到 +0.02 Gate；**filter 方向 CONDITIONAL NEGATIVE**
- **視覺**：
  - 主視覺：圖 C1（`fig_c1_auc_before_after.png`）AUC before/after 雙 panel
  - 右側小條：4 個特徵降幅 ≥−0.025（HPFineNGroups、NHP3、HPMergedDelta、HPFine_NGroups_CF）
  - AlleleDelta 紅圈標註：HP-independent 對照 ΔAUC = 0.0000 ✓（預期一致）
- **口頭**：
  - Gate 設定 ≥+0.02 是 plan 預期的最低合理增益（此值低於 Beyond-AUC 耗盡 0.58 ceiling 的驗收門檻）
  - 實測：最大正向 ΔAUC = +0.0084（HP1FamilyN），遠不及 +0.02
  - 四個主要 HP-derived 特徵**反向下降**這個事實值得深思 —— 不只是「沒增益」，而是「去除 somatic HP tag 後，fine-group 結構本身消失」
  - Quality_Score 目標 ≥0.55 FAIL（實測 0.5148）→ 方法學獨立，見下一 slide
- **捨棄**：
  - 18 特徵完整表格（speaker notes 列全表，slide 只放圖）
  - 反向降幅的 p-value / bootstrap CI

---

### Slide 15 · HPFineNGroups 分佈崩潰（NG≥3 = 0）

- **核心訊息**：flag=on 下，原本 82% TP regions 的 HPFineN N≥3 完全消失
- **視覺**：
  - 主視覺：圖 C2（`fig_c2_ngroups_distribution.png`）NGroups 分佈堆疊對照
  - 關鍵數字：TP 22,732 regions + FP 8,148 regions 從 N=3/4 → N=2
  - 粗體警語：「0418 F pilot 的 **NG=4+NR≥80 subclone marker 訊號源至少部分依賴 somatic HP tag**」
- **口頭**：
  - 這是**雙重意義**的觀察：
    1. 證實 Phase 1 機制正確（flag 確實消除 somatic-derived HP bucket）
    2. **暴露了 HPFineN marker 的根基問題**
  - 0418 F pilot 結論的 89.1% TP rate 是在 master dataset（全 7 樣本 pooled、含 AF filter）上測的；本週 HCC1395 TO raw ClairS-TO split 無法重現（TP rate 0.6944 vs baseline 0.699，Fisher odds ratio 0.913 反向 p=3.5e-3）
  - 兩種解讀並存：(A) 資料集 dependency（master 結論仍成立，TO standalone 無富集是 pipeline 差異）；(B) 訊號部分來自 somatic HP 人工分組（flag=on NG≥3=0 orthogonal null）
  - **這個觀察是 Thread D 機制揭露的關鍵伏筆**
- **捨棄**：
  - Fisher odds ratio 計算細節
  - Marker re-audit 完整表（speaker notes 附）

---

### Slide 16 · HP_Ratio shift 與 R3 條件未觸發

- **核心訊息**：HP_Ratio 變化小，但方向正確；**Plan R3 條件（HP_Ratio 歸零）未觸發 → bug 不在上游**
- **視覺**：
  - 主視覺：圖 C3（`fig_c3_hpratio_shift.png`）HP_Ratio dumbbell shift
  - 四層數字：TP all / TP Non-LOH / TP LOH / FP all，每層 off vs on
  - 結論 bullet：「Plan 預期 0.836 → 0.55-0.65 不成立；baseline 本就已較低」
- **口頭**：
  - 最大 shift 僅 −0.023（Non-LOH TP）；LOH stratum 已極端不平衡（0.091）幾乎不動
  - Plan 引用的 0.836 baseline 可能來自舊版 haplotag 或不同 VCF 子集；V3-Fixed TO BAM 的 HP_Ratio 本就已較低
  - Plan R3 原條件：「修正後 HP_Ratio 仍偏 0 → bug 在 hp_tag 解析」未觸發
  - **定位確認**：bug 不在上游 hp_tag 解析；無需升級為 LongPhase-TO 上游 C++ 修復
- **捨棄**：
  - HP_Ratio 計算公式
  - 其他分層細節（flag 對 Potential_LOH=True/False 的非對稱效應）

---

### Slide 17 · 「Tag 可採用」的定位（用戶決策 A）

- **核心訊息**：機制 ✓ + filter 懸掛 + characterization 保留
- **視覺**：
  - 三欄狀態卡
    - ✅ 機制正確性（H-C1）：Audit independence 通過、HP_Ratio 方向正確、NG≥3 消失
    - ⏸ Filter 增益（H-C2）：HCC1395 TO Gate FAIL，需 Phase 2 全樣本才能判定
    - ⏸ Paired 一致性（H-C3）：本週未完成對照（需 Phase 2 master × 兩 flag）
  - 結論標語：**「技術可行但目前 filter 無增益；Phase 2 判定前結論懸掛」**
- **口頭**：
  - 用戶主軸原話「驗證效果合理可行，tag 可採用」 —— 依決策 A (b)，本週採「延後定論」定位，不宣告「可採用」
  - default=off 保留；flag=on 作為研究工具提供（乾淨 HP 分群用）
  - 下週 P0 行動：master dataset × 兩 flag × 7 樣本重跑 → 才能回答 H-C2/H-C3
- **捨棄**：
  - "tag 可採用" 肯定語（已明確避用）
  - Phase 1.5 within_dom_alt_frac downstream 重建（成本高、價值低，已判定 NOT DO）

---

## Act III · 機制揭露與 pivot（Slides 18-24）

### Slide 18 · 用戶提問觸發機制重查（Q4 引入 · 偵探時刻）

- **核心訊息**：一句提問「NG=2 與甲基有關係嗎」觸發 C++ 原始碼回查 → 發現誤解
- **視覺**：
  - 對話框視覺化（類似漫畫氣泡）：
    - **Past / AI V1 (2026-04-22 15:56)**：HPFineNGroups = "fine-grained methylation N groups" → 論文主軸 "Haplotype-loss-dependent methylation bimodality"
    - **User (2026-04-22 晚間)**：「NG=2 與甲基有關係嗎？」
    - **AI V2 (2026-04-22 20:40)**：C++ 原始碼回查 → **無 methylation 依賴** → 重新定性
  - 時間戳標註：**4 小時 44 分**的關鍵轉折
- **口頭**：
  - 這是本週**偵探時刻**，也是 memory 新增 `feedback_feature_name_vs_definition_rule` 的起源
  - Feature name 直覺解讀陷阱：`HPFineNGroups` 字面像 "HP-Fine-N-Groups" = "methylation N groups"，但實際是 HP × variant 的 4-bucket occupancy
  - 教訓：新 feature 分析前必查 C++ 原始碼，不可依字面意思推論生物語意
- **捨棄**：
  - V1 錯誤推論的完整邏輯鏈
  - 重構論文主軸時考慮過的候選（如 "allele-specific methylation skewing"）

---

### Slide 19 · HPFineNGroups 真機制（C++ 原始碼證據）

- **核心訊息**：NG 純粹是 `{HP1, HP1-1, HP2, HP2-1}` 4-bucket occupancy count，與 methylation 無直接關係
- **視覺**：
  - **C++ 程式碼片段**（少見地放視覺區，因為是關鍵證據）：
    ```cpp
    // src/core/LabelTest.cpp:265-302
    if (hp == "1")   group = 0;  // HP1    = hap 1, ref allele
    else if (hp == "1-1") group = 1;  // HP1-1 = hap 1, somatic allele
    else if (hp == "2")   group = 2;  // HP2    = hap 2, ref allele
    else if (hp == "2-1") group = 3;  // HP2-1 = hap 2, somatic allele
    ```
  - 右側註解引用：`include/core/Stats.hpp:324` 設計者註：`"Same haplotype, germline vs somatic"`
  - 底標：**Methylation 僅在 HPFineF/HPFineP 的 PERMANOVA 中作品質檢驗，NGroups 本身不依賴**
- **口頭**：
  - 這張 slide 的 C++ 程式碼是為**證據鏈的終點**放的；PI 看到程式碼比任何統計圖都直接
  - 設計者當初用 "FineN" 是指 fine granularity in phasing，不是 methylation；文件遺漏導致後續誤讀
  - Stats.hpp 的 "Same haplotype, germline vs somatic" 是當初 d_hp1_hp1s 距離指標的設計意圖
- **捨棄**：
  - 完整的 `hp_to_fine_labels()` 實作（238 行）
  - `HPFineF` / `HPFineP` PERMANOVA 的統計機制

---

### Slide 20 · NG=2 的四種組成可能

- **核心訊息**：NG=2 意指「4 bucket 中 2 個被 populate」；組成決定生物意義
- **視覺**：
  - 2×2 格子視覺化 4 bucket + 4 種 NG=2 組成
    - **same_HP1** (HP1+HP1-1)：單 hap 內部 ref+somatic → 物理必為 somatic
    - **same_HP2** (HP2+HP2-1)：同上（對稱）
    - **cross_het** (HP1+HP2-1)：canonical het phasing → germline het 與真 somatic 不可分
    - **cross_het_inv** (HP1-1+HP2)：對稱
  - 底標關鍵區分：**same-hap = somatic signature；cross-het = germline leak 風險區**
- **口頭**：
  - 組合拆分是 obs18 腳本的邏輯基礎；`research/tpfp_loh_af_kde_discrimination/scripts/obs18_TO_NG2_composition.py`
  - 物理推論：LOH 區只有單 haplotype → 必產 same-hap 分裂；非 LOH 區雙 hap 保留 → cross-het 是主要 pattern
- **捨棄**：
  - "other" 類（非以上四種的 NG=2）的細節
  - 各組合在 caller 算法中的生成機率

---

### Slide 21 · obs18 跨 6 樣本驗證（核心數據）

- **核心訊息**：**Inner × NG=2 在 6 TO 樣本全部 ≥93% 為 same-haplotype**（6/6 一致）
- **視覺**：
  - 主視覺：圖 D1（`obs18_NG2_composition_proportion.png`）6 樣本 Inner vs Outer 堆疊比例
  - 右側大字數字：`93.2% · 99.0% · 98.8% · 96.5% · 98.3% · 97.0%`（6 樣本 Inner same-hap%）
  - Median 97%（粗體）
- **口頭**：
  - 圖中每個樣本兩根柱（Inner / Outer）；Inner 柱綠色塊（same-hap）幾乎滿格
  - Outer 柱組成多變（cross-het 占 0.1-97.5%，median 44%）—— 但**主旨不在 Outer 組成，在 Inner 的跨樣本一致**
  - HCC1395_DORADO 的 Outer 也 97.5% same-hap 是 ONT DORADO basecaller 的特殊 phasing pattern，不影響主結論
- **捨棄**：
  - 每樣本 Inner / Outer 的 total n（在 speaker notes 和週報證據卡已記）
  - "other" 組成詳細比例

---

### Slide 22 · Inner vs Outer TP gap +0.37（6/6 正向）

- **核心訊息**：Inner same-hap TP rate median 0.93 vs Outer cross-het 0.55 → **gap +0.37 median, 6/6 正向**
- **視覺**：
  - 主視覺：圖 D2（`obs18_NG2_composition_heatmap.png`）TP rate heatmap
  - 左：Inner same_HP1 TP rates（0.96/0.94/0.76/0.43/0.93/0.92）
  - 右：Outer cross_het TP rates（0.50/0.55/0.24/0.08/0.88/0.69）
  - 中央 gap 箭頭：+0.46 / +0.39 / +0.52 / +0.35 / +0.05 / +0.23
- **口頭**：
  - 0 反向；H2009 gap +0.05 是 baseline 飽和（Inner 和 Outer TP 都很高，差距縮小）
  - HCC1954 outlier：Outer cross-het TP 0.08 極低 —— Potential_LOH 可靠性需下週 P1 專項分析（如果 HCC1954 的 "Outer" 其實含有被 LOH.bed 漏標的真 LOH，那 cross-het 的歸類本身有 bias）
  - Formal stats：Wilcoxon signed-rank on 6/6 gap 下週補齊
- **捨棄**：
  - Inner same_HP2（對稱版本）的數字
  - Wilson CI 的計算方式

---

### Slide 23 · LOH-constrained phasing 機制圖（論文主軸 pivot 核心）

- **核心訊息**：兩種物理場景 —— LOH 區 somatic 必產 same-hap；非 LOH 區 germline-somatic phasing 同構不可分
- **視覺**：
  - 雙面板機制圖（自繪，或用週報的 Mermaid 節點）
    - 左面板「LOH 區（Inner）」：單 hap → somatic SNV → HP1+HP1-1 → 綠框「物理必為 somatic TP」
    - 右面板「非 LOH 區（Outer）」：雙 hap → somatic het OR germline het → 都產 HP1+HP2-1 → 紅框「germline-somatic 不可分 → TO germline leak 物理根源」
  - 底標粗體：**新論文主軸（TO 層）：LOH-constrained phasing signatures distinguish somatic from germline-like variants in tumor-only sequencing**
- **口頭**：
  - 這是本週 TO 層論文主軸從 "methylation bimodality" 正式 pivot 的瞬間
  - 優勢（vs 舊主軸）：
    - 機制純 phasing → 無需 methylation 實驗驗證
    - 可直接從 LongPhase output + LOH bed 重現
    - 跨 basecall version（DORADO）、跨 pipeline（LongPhase-TO）穩定
    - 連接 FACETS/Battenberg/TITAN 既有文獻模型
  - **paired 層決策**：依用戶 C · 以 TO 為本週重點；paired mode AF × NGroups POSITIVE 結論保留不撤回，加註「需獨立 phasing-vs-methylation 驗證」
- **捨棄**：
  - 舊 "methylation bimodality" 的詳細誤解敘述（已在 slide 18）
  - 對 `project_loh_subclone_af_methylation_positive` memory 的具體改寫建議

---

### Slide 24 · Thread C 是 Thread D 的獨立 negative control

- **核心訊息**：flag=on 下 NG≥3=0 不只是 Thread C 的 filter 失敗，同時是 Thread D 機制的獨立佐證
- **視覺**：
  - 雙向箭頭視覺：左側 Thread C (`--germline-hp-only` flag=on) → NG≥3 = 0 (全 40k sites)
  - 右側 Thread D (obs18 flag=off) → Inner × NG=2 ≥93% same-hap
  - 中央關鍵句：**「somatic HP demote → NG≥3 消失 → NG 是 phasing 而非 methylation 的最強證據」**
- **口頭**：
  - Thread C 的 filter FAIL 原本只解釋為「tag 方案無增益」；加上 Thread D 機制後，重新解讀為：
    - NG 本身是 phasing bucket；demote somatic → 沒有 HP1-1/HP2-1 bucket → NG 只能落在 {0, 1, 2}
    - 這是**機制正確性的強驗證**，甚至比 obs18 的跨樣本一致更直接
  - 4 假說全 POSITIVE，結論穩定度 ⭐5（TO 層；paired 層需獨立驗證）
- **捨棄**：
  - 額外的 confound 討論（spatial autocorrelation、coverage effect）— 都已在 C++ 原始碼層被排除

---

## 結語（Slides 25-26）

### Slide 25 · 結論總表變動

- **核心訊息**：6 項結論變動 —— 3 新增 ⭐5 / ⭐4、1 降級 ⭐4→⭐3、1 加註保留、1 新增 ⭐3 懸掛
- **視覺**：表格（6 列 × 3 欄：結論 ID / 舊穩定度 / 新穩定度）
  - CL-LCP-001 (新 ⭐5)：**LOH-constrained phasing signatures** — 本週 pivot 核心
  - CL-S3-001 (新 ⭐4)：**S3 Diploid Het + S5 combo biology-informed filter**
  - CL-CN-KDE-001 (新 ⭐4)：CN KDE 雙 Pass 校準 bias <2%
  - CL-016 (⭐4 → ⭐3)：HPFineN subclone marker pipeline dependency 警告
  - CL-LAF-001 (⭐4 加註)：paired AF × NGroups 保留 + 註記
  - CL-HP-ONLY-001 (新 ⭐3)：`--germline-hp-only` Phase 1 mechanism ✓ / filter 懸掛
- **口頭**：
  - ⭐5 結論：C++ 原始碼 + 6/6 obs18 + Thread C 獨立 negative control 三重證實
  - ⭐4 結論：HCC1395 pilot 完整，待全量 7 樣本收斂升 ⭐5
  - ⭐3 結論：單樣本或懸掛
- **捨棄**：
  - 所有 ⭐1/⭐2 未變動結論
  - Tier 計算方法細節

---

### Slide 26 · 下週 P0-P2 優先行動 + 收斂圖

- **核心訊息**：4 項 P0/P1 行動下週完成，2 週內定案 LOH-constrained phasing 論文主軸
- **視覺**：
  - 左：里程碑時間軸（週報 Layer 4 圖），突出 2026-04-30 與 2026-05-07 兩節點
  - 右：P0/P1 行動卡片（4 項）
    - **P0** Paired normal-pilot obs18 對照 → 驗證 H-D3 gap 是否消失（1-2 天）
    - **P0** Archive TO 6 樣本重跑 ISM → 跨樣本 S1-S7 + marker × 兩 flag（~10 hr parallel）
    - **P1** HCC1954 outlier 根因分析（Potential_LOH 可靠性）（1 天）
    - **P1** Formal stats: Wilcoxon signed-rank on 6/6 gap（0.5 天）
- **口頭**：
  - 關鍵風險：(a) Paired 對照若 gap 未消失 → H-D3 修正 → Thread D 穩定度降為 ⭐4 TO-only；(b) master × 兩 flag 若 HPFineN marker TP rate 暴跌 → marker 生物學根基需重建
  - 緩解：已在風險表記錄；paired 層保留不撤回的決策提供 fallback
  - 2 週內兩個 milestone：下週 paired 對照 + Archive 重跑；再下週 Phase 2B master rerun × 兩 flag
  - 3 週後：S4 二級判別 + Wakhan pilot + ISM 論文 Figure 3/4 重寫
- **捨棄**：
  - P2/P3 延伸方向（Wakhan / S4 二級判別 / paired phasing-vs-methylation orthogonal test）放在 speaker notes

---

## Storyboard 自檢（導演審查）

### 敘事邏輯驗證

```
Q1 (CN 不準) → Slide 5-7（Act I 前段）
    ↓ 測量單位校準完成
Q2 (LOH×AF×CN 能切到什麼程度) → Slide 8-9（Act I 後段）
    ↓ S4 bucket 卡關
Q3 (HP 雜訊拖累？) → Slide 10-17（Act II）
    ↓ 機制 ✓ filter FAIL → 但 NG≥3=0 是偵探線索
Q4 (NG 到底是什麼) → Slide 18-24（Act III）
    ↓ C++ 原始碼回查 → obs18 6 樣本 → 機制更正 → 論文 pivot
結語：結論變動 + 下週（Slide 25-26）
```

**因果鏈無跳步** ✓。每一幕的轉折點都有明確驅動：Act I → Act II（S4 卡關）；Act II → Act III（NG≥3=0 觀察 + 用戶提問）。

### 重點平衡（按主軸 token 數）

| Thread | 主軸 tokens | Slides | 比例 |
|:-----:|------------|:------:|:---:|
| A (CN KDE) | ~120 字 | 3 slides (5-7) | 12% |
| B (LOH×AF×CN filter) | ~200 字 | 3 slides (8-10) | 12% |
| C (HP 雜訊 + germline-hp-only) | ~280 字 | 7 slides (11-17) | 27% |
| D (NG=2 LOH-constrained phasing) | ~350 字（主軸高潮）| 7 slides (18-24) | 27% |
| 定錨 + 結語 | — | 6 slides (1-4, 25-26) | 23% |

Act III (D) 佔比最高，符合本週**機制揭露為主旋律**的重點。Act II (C) 篇幅與 Act III 相當，因為 C 內部分兩層：(1) Phase 1 驗證邏輯（5 slides）；(2) tag 定位決策（slide 17）。

### 視覺/口頭/捨棄三層分配檢核

- **視覺層**（每 slide ≤4 bullet，視覺占 55-65%）：核心數字、關鍵圖、機制示意；避免完整 AUC 表、CI 計算、原始碼全文
- **口頭層**（speaker notes）：統計檢定細節、replaced 的替代假說、pipeline 版本關係、memory 教訓
- **捨棄層**：
  - 所有 ⭐1/⭐2 未變動結論（已關閉不重訪）
  - ClairS-TO Verdict 細節（0421 已述）
  - R1-Global 完整 CI
  - Wilson CI / Fisher odds 的計算方式
  - 0406 LOH 雙定義關閉的完整歷史
  - LongPhase-TO 與 LongPhase-S 的差異（Knowledge Base 已述）

---

## 後續產出（Phase 4 需同步生成）

1. `01_full_narrative_report.md` — 完整敘事報告（~2500 字，逐 slide 對應段落）
2. `02_ppt_slide_outline.md` — 每 slide 的 bullet / 圖 ID / notes 的結構化 JSON-like 大綱
3. `03_slide_layout_and_script.md` — 每 slide 的佈局座標 + 雙語文案
4. `pptx_config.json` — build_pptx.py 的配置（色彩、字體、aspect ratio）
5. 渲染 26 slides 為 PNG → 截圖驗證 aspect ratio / 文字溢出 / 色彩

---

**文件結束** · Phase 3.5 導演 Storyboard 定稿候選 · 待用戶審閱後進 Phase 4
