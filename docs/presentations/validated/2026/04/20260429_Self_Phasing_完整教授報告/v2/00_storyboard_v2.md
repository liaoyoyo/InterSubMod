# Self-Phasing PPT v2 · Storyboard（**24-slide v3**, 2026-04-30）

## 0. 報告定位

| 項目 | 設定 |
|------|------|
| 受眾 | 教授（NGS / long-read 基本背景，不熟 LongPhase-TO 與 InterSubMod 切分） |
| 風格 | 演講者模式：slide 文字最少（標題＋視覺）、細節入 speaker note |
| 主軸 | 用戶 6 段流程：高層次重點 → 觀察問題 → 為何重要 → 解釋與原因 → 改動驗證 → 未來規劃 |
| Slide 數 | **24**（v3 結構定案） |
| 演講節奏 | 30 min × 75 sec/slide |
| 比例 | 16:9，Latin Arial + CJK Droid Sans Fallback per-char |

## 0.1 修訂歷程（v1 → v2-32 → v2-26 → v3-24）

| 階段 | slide 數 | 主要差異 |
|------|:-:|---------|
| v1（2026-04-29 上午）| 34 | 線性敘事 |
| v2-32（2026-04-29 下午第一稿）| 32 | 6 段流程；補業界、五大目標 |
| v2-26（2026-04-29 reviewer 後）| 26 | 結論先行單頁、合併冗餘、補 prerequisite |
| **v3-24（2026-04-30 用戶 6 Section 逐項確認 + 根因敘事整合）** | **24** | 合併 5+6/11+12/22+24；拆 8/17；S2 業界吸收進 S20；**S11 加觸發條件主者**；Q1/Q2 整合 speaker notes |

## 0.2 v3 關鍵新增（最新根因敘事）

> **純樣本被 baseline purity calculator 估為 0.927（四捨五入 0.93），低於純樣本標準 0.95，因此「沒有」啟動 Two-Pass（兩條路）方法，走 baseline 標準三條路流程；因樣本實為純無 normal contamination 而流程假設有 → 暴露 tag 層 somatic-first 投票 bias → self-phasing 17.3:1。**

驗證來源：`v5_audit_suite/18_purity_calculator_failure_root_cause.md` + `PhasingProcess.cpp:197` 硬編碼 `purity > 0.95` 閾值 + `HaplotagProcess.cpp:512-563` getVote() somatic-first 順序。

---

## 1. 大架構（6 段、24 slides）

```
Section 1 · 高層次重點（2 slides · S1-S2）
Section 2 · 觀察到問題（5 slides · S3-S7）
Section 3 · 為何重要（2 slides · S8-S9）
Section 4 · 解釋問題與發生原因（4 slides · S10-S13）
Section 5 · 改動驗證與結論（7 slides · S14-S20）
Section 6 · 未來目標與規劃 + 結語（4 slides · S21-S24）
```

關鍵設計：**S1 衝擊式 TL;DR + 動機小字**；**S11 baseline 根因樹含觸發條件主者**；**S17 V5max1 climax + caveat**；**S20 業界家族樹 + 2x2 合併**；**S24 倒過來 take-home + Q&A 感謝**。

---

## 2. 逐 slide storyboard

### Section 1 · 高層次重點（2 slides）

| # | 標題 | 視覺 | speaker note 重點 |
|---|------|-----|-----|
| **S1** | 17.3 : 1 → 1 : 1（重點先講）| 大字「17.3:1 → 1:1」+ 綠字「+8.3 pp（clean PS, 全基因組）」+ 灰字「self-phasing 分層處理：V2b 解 phase scaffold；V3F/V5 解 tag 層」+ **底部小字「修補 self-phasing 是解鎖 ISM 五大研究目標的前提」** | 開場 60 秒鎖定核心結論；V5 不修 self-phasing 本身；InterSubMod 無 C++ 改動 |
| **S2** | 做了哪些事 — Pipeline + 4 階段 | 上半 Pipeline schematic（tumor BAM → ClairS-TO → longphase-to-mod[V5紅標] → tagged BAM → InterSubMod）+ 下半 4 階段（機制定位 → 4-commit 修補 → 驗證 → 影響評估）| 強調 InterSubMod ≠ 修補位置；本 repo 無 C++ 改動 |

### Section 2 · 觀察到問題（5 slides）

| # | 標題 | 視覺 | speaker note 重點 |
|---|------|-----|-----|
| **S3** | 基本定義 — HP tag 五整數值 + PS / LOH 兩層 | `figures/fig17_hp_tag_5versions.png` + slide 底 1 行 legend | HP1/HP2/11/21/33/0；PS = 同 phase block reads ID；LOH **兩層**先預告 |
| **S4** | 證據三層：理論 1:1 + 全基因組 17.3:1 + 個別位點 SP1 113:0 | 三欄並列：左理論隨機分散 schematic / 中 fig01d Panel D 全基因組 614K vs 35.5K / 右 SP1 IGV 縮圖 | 為何 17.3:1 是 artifact；94.6% 集中於 HP1；SP1 個別位點達 113:0 完全失衡 |
| **S5** | Paired 對照確實 ~1 : 1 | 自繪 bar：Paired (HP_Ratio 0.5) vs TO (HP_Ratio 0.94) | 同位點跨模式 r=0.001、Cohen's d=−1.20；預先回答「那不就用 paired 就好」→ TO 是必要研究方向 |
| **S6** | 拆鎖：phasing 沒問題，是 tag 的問題 | 兩欄表（phasing 層 LOH.bed Jaccard=1.0 ✓ / tag 層 HP_Ratio 17.3:1 ✗）| LOH.bed 由 VCF AD 產生、HP_Ratio 由 BAM HP tag 產生；本 PPT 焦點在 tag 層 |
| **S7** | LOH 兩層精確化 | 兩欄表：左 ISM HP_Ratio LOH（62% artifact）/ 右 LOH.bed region-level LOH（Jaccard=1.0 不變）| kappa=0.670 不一致性由「不同數據源」解釋；本工作只影響 ISM HP_Ratio LOH 層 |

### Section 3 · 為何重要（2 slides）

| # | 標題 | 視覺 | speaker note 重點 |
|---|------|-----|-----|
| **S8** | 為何重要 — TO 根基 + ISM 強依賴 HP tag | 左：示意 TO 純數據是 TP/FP 觀察根基；右：4-bucket 分群（HP1 / HP1-1 / HP2 / HP2-1）+ NGroups/HPSig/HP_Ratio | 整個 ISM 研究都建立在 TO 模式；4-bucket 是 ISM 核心 |
| **S9** | 這些 tag bug 影響哪些 ISM 特徵（共 85）| `figures/03_self_phasing_impact.png`（ISM 3-tier）+ slide 標 count 「29 / 14 / 42」 | 🔴 嚴重 **29 個** 必重跑；🟡 中度 **14 個**；🟢 不受影響 **42 個**；**口徑校準**：來源報告誤寫為 38%/7%/55%，slide 上只列 count |

### Section 4 · 解釋問題與發生原因（4 slides）

| # | 標題 | 視覺 | speaker note 重點 |
|---|------|-----|-----|
| **S10** | Prerequisite — 函數架構 + PON 概念（合併） | 左 4-row 表：getVote / judgeHaplotype / countSNPHaplotype / countINDELHaplotype；右 schematic：tumor-only 無 normal → PON 取代 → V2b `--pon-only-phasing` 啟用 | 教授不熟 longphase 內部架構必補；speaker note 補：**Purity 0.95 閾值** 概念（baseline 內建 Two-Pass 觸發鎖，`PhasingProcess.cpp:197`）|
| **S11** | ★ **Baseline 根因樹（含觸發條件 + 三 bug）** | 自繪 root-cause tree：**根 = 觸發條件「Purity 0.927 ≤ 0.95 → 未觸發 Two-Pass → 走三條路」** → 三葉 bug（① getVote somatic-first priority ② enum vs int literal mismatch ③ Phase scaffold self-phasing[V2b 已解]）+ 紅標代碼摘錄 | **新增觸發條件主者**：純樣本被估為 0.927，低於 0.95 標準，因此**沒**啟動 Two-Pass → 走 baseline 三條路標準流程 → 因樣本實為純無 normal contamination 而流程假設有 → 暴露 tag 層 somatic-first 投票 bias；三 bug 用顏色區分（紅 priority / 橘 enum / 灰 scaffold）|
| **S12** | ★ V5 三層投票流程（細化定位） | 自繪流程圖（V5 三層 with 紅圈標出對應 S13 diff 段落）| Layer 1 **germline-first 投票順序**（baseline 為 somatic-first 的關鍵翻轉）→ Layer 1.5 somatic fallback（**confidence 0.6 是 mid-low purity 防禦層**，針對所有 purity 範圍特別 0.6-0.93）→ Layer 2 encode 11/21/33；**Q2 答案整合**：V5 在 mid-low purity 比 baseline 更好（09 報告 0.6 sample HP33 比例 12.4% vs Baseline 2%，保守 +10pp 避免錯誤分配） |
| **S13** | **程式碼 diff slide**（baseline vs V5）| 兩欄 monospace 24 pt（每欄 ≤ 6 行）+ baseline 欄上方紅標 caveat **「somatic-first 投票順序」** + 紅標 priority bug + enum bug 修正點 | 用戶指定 1 張 diff slide；左：baseline `for(variantKeys){if any non-zero break}`（**somatic-first**）；右：V5 `if germline > 0 ... else if somatic > 0 ...`（**germline-first**）；標 enum bug 修為 integer literal 11/21；**baseline 並非「有 HP11/HP21 即 somatic」這個簡化敘述，而是 somatic-first 投票優先序設計缺陷** |

### Section 5 · 改動驗證與結論（7 slides）

| # | 標題 | 視覺 | speaker note 重點 |
|---|------|-----|-----|
| **S14** | ★ Sanity 4 項 PASS + 5 條獨立證據鏈 | 公式卡（守恆律 A/B + Layer 1.5 期望 1/2）+ 15/15 PASS 標章 + 右欄 5 條證據鏈摘要 | 守恆律 A/B PASS；**Q1 答案整合**：PON-only Two-Pass 在 0.93 sample（V5 88.2% vs BL 74.9%, +13.3pp）與 0.6 sample（09 報告：V5 HP tag 分布更接近真實）**都有效**；aggregate +6.65pp；clean PS +13.3pp |
| **S15** | 量化指標 — AMB / HP33 / Concordance | `figures/fig18_concordance_amb_f1.png` | AMB 17.5→8.0% (-9.5pp)；HP33 reads 239,679→110,197 (-54%)；clean PS paired GT V5 90.5% vs Baseline 82.2% (+8.3 pp **clean PS blocks 全基因組**)|
| **S16** | 4-commit 演進 timeline | `figures/fig01a_commit_evolution.png` | V2b 8b8c1fd PON-only → V3F 41ff147 兩層 + enum fix → INDEL guard 380e8d2 → V5 working tree Layer 1.5 |
| **S17** | 🎯 CLIMAX：V5max1 — 39 reads, 100% reassigned | `figures/C_V5max1_chr19_4639528.png`（圖佔 80% 寬度）+ 標題大字「39 reads, 100% reassigned」+ 底部小字 caveat **「但 V5 不修 self-phasing 本身（V2b 已處理 phase scaffold）」** | 演講高潮；**reviewer 警示防誤解**：V5 修的是 tag 層 directional reassignment，不是 phase 層 self-phasing |
| **S18** | HP 翻轉證據 — Bar chart + IGV 縮圖（climax 延伸）| 自繪 bar chart（baseline vs V5 vs paired，3 SP-extreme 位點）+ 右側小 IGV 縮圖（D_SP1）| 打破 S17 IGV 視覺單調；3/3 SP-extreme 一致 |
| **S19** | 矛盾或盲點 — F1 釐清 + cnLOH 邊界 | F1 流程圖（ClairS-TO raw F1=0.7166 三版本相同 → 經 ISM SuggestFilter → ΔF1=-0.0003 噪音）+ 4 條 caveat | F1 不能衡量本實作品質、真實價值在 read-level；cnLOH 區仍未解；AF>0.9 邊界；Confidence threshold 0.6 未直接驗證 |
| **S20** | 業界家族樹 + 2x2 比較表（合併原 S2 + S21）| 上半家族樹（LongPhase 2022 → LongPhase-S 2025 paired / **longphase-to-mod V5** TO+PON）+ 下半 2x2 matrix（Germline/Tumor × 公開工具/本實作）| 同實驗室相鄰工作（**非業界共識**）；本實作填補 tumor-only 在公開工具中的 gap；LongPhase-S 在 ClairS 上 SNV F1 +4.5% 為 paired 場景參考（**注意非直接可比，因本工作 F1 不變是預期行為**）|

### Section 6 · 未來目標與規劃 + 結語（4 slides）

| # | 標題 | 視覺 | speaker note 重點 |
|---|------|-----|-----|
| **S21** | 後續可動 / 已初步現（合併原 S22+S24）| 上半 4-row「可一起做」（Phase 2A normal methylation / HPFineNGroups marker 重驗 / archive haplotag_version / Thread D main axis）+ 下半 3 卡片「已初步現」（Thread D NG=2 6/6 POSITIVE / HPFineNGroups 重驗中 / LOH-constrained phasing 論文主軸候選）| 後續可推動方向 + 已初步發現 |
| **S22** | 五大目標銜接 + 研究發展樹 | 上半 5 卡片（per-CpG / clone / two-hit / TO normal / F1 panel）+ 下半 Mermaid 樹（self-phasing 修補 → 4-bucket 可信 → NG analysis 可信 → 五大目標解鎖）| 來源：`InterSubMod/docs/architecture/20260327_InterSubMod研究願景定錨_01.md`；目標 1/2/4 直接依賴 HP tag 正確 |
| **S23** | 短期 P0 行動清單 | 3 條 P0：F1 commit V5 working tree（高）/ F3 7 樣本 V5 全量重跑 / F4 master×flag 重驗 HPFineNGroups | F5-F8 入 speaker note |
| **S24** | Take-home + Next step + Q&A 感謝 | 上半大字 take-home：「TO 模式可用 · Tag 層問題已解 · 五大目標解鎖」+ 下半 3 條 next step 承諾 + **底部小字「Q&A 歡迎；詳細數據見 source_materials/」+ Thank you** | 不重複 S1；以承諾收尾比口號收尾更有力；強化結尾 |

---

## 3. v3 變動 vs v2-26 對照

| 變動類型 | 細節 |
|---------|------|
| **合併** | 原 S5+S6 → S4 證據三層；原 S11+S12 → S10 prerequisite；原 S22+S24 → S21；原 S2 業界 → 吸收進 S20 |
| **拆分** | 原 S8 → S6 拆鎖 + S7 LOH 兩層；原 S17 → S15 量化 + S16 timeline |
| **新增** | S11 root-cause tree **加觸發條件主者**「Purity 0.927 ≤ 0.95 未觸發 Two-Pass → 走三條路」 |
| **強化 speaker note** | S10 補 Purity 0.95 閾值；**S12 補 mid-low purity 防禦層 + Q2 答案**；**S14 補 Q1 PON 雙路徑有效性**；S13 baseline 欄紅標 somatic-first |
| **強化結尾** | S24 加底部小字「Q&A 歡迎」+ 感謝 |
| **數量** | 26 → **24 slides**（減 2）|

## 4. v3 截圖驗證檢查項

每張 slide 必須通過：

1. 文字無溢出（Latin + CJK 雙語不互相覆蓋；標題 ≤ 15 字；每頁 ≤ 5 視覺元素）
2. 圖片等比 fit-within（無強制縱橫比破壞）
3. CJK 字型正確（無方塊、無亂碼）
4. Speaker note 存在且 ≥ 350 字（**S11/S12/S14 因含 Q1/Q2 整合 ≥ 500 字**）
5. 圖片連結有效

驗證方式：`v2/scripts/screenshot_all.py` 必須輸出 `[ISSUES] none detected`。

## 5. v3 5 處口徑校準（reviewer 第二輪採納）

| # | 校準項 | 落實位置 |
|---|------|---------|
| 1 | Self-phasing 分層語法（V2b/V3F/V5）| S1/S17/S24/speaker script |
| 2 | 「同實驗室相鄰工作」非「業界共識」| S20/speaker script |
| 3 | ISM 3-tier 只列 count（29/14/42 共 85）| S9 標註「來源報告誤寫」 |
| 4 | enum 行號 `Util.h:20-25` | S3/S11/speaker script |
| 5 | concordance caveat「clean PS blocks，全基因組」| S1/S15/S24 |

## 6. v3 完成狀態

| 階段 | 狀態 |
|------|------|
| Storyboard 24 slides + v3 變動 | ✅ 完成（本檔）|
| Source materials（含 18 報告 04_purity_calculator_failure...md）| ✅ 完成 |
| Figures 33 張實體複製 | ✅ 完成 |
| 6 份 onboarding 文件（將更新 background/key_metrics/code_references/qa）| 🔄 進行中 |
| `scripts/build_pptx.py` 24-slide 重寫 | ⏭ 待 build 階段 |
| `scripts/screenshot_all.py` | ✅ 既有可重用 |
| `output/*.pptx` 24-slide 版 | ⏭ 待 build 階段 |
| `notes/speaker_script_v2.md` 24-slide 版 | ⏭ 待 build 階段 |

## 7. v3 build 必做清單

1. ⏭ build_pptx_v3 撰寫（重用 v2-26 helper + 新增/調整 schematic 適應 24-slide）
2. ⏭ speaker_script_v3 撰寫（每張 ≥ 350 字；S11/S12/S14 ≥ 500 字含 Q1/Q2 整合）
3. ⏭ screenshot 驗證 `[ISSUES] none detected`
4. ⏭ 抽查 5 張關鍵 slide：**S1（衝擊 TL;DR + 動機小字）**、**S11（根因樹含觸發條件）**、**S13（程式碼 diff baseline somatic-first 紅標）**、**S17（climax + caveat）**、**S24（take-home + Q&A 感謝）**

完整 build 指令：

```bash
cd InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v2
python3 scripts/build_pptx.py
python3 scripts/screenshot_all.py
```
