---
title: PPTX 18-slide Critical Issues Root Cause + 重新設計計畫
type: pptx_redesign_analysis
date: 2026-05-06
status: pending_user_decision
---

# 1. Critical Issues Root Cause（5-Why）

## 1.1 P0 #1：HP:i:33 數字陳述錯誤（比 Agent-D 報告的還深）

### 真相（來自 source 交叉驗證）

| metric 名稱 | 定義 | 數量級 | 出處 |
|------|------|------|------|
| **`HP:i:33`** | 4/29 audit / PI 報告 4 column header | 帶 SAM tag `HP:i:3` 的 read 計數（whole genome）| **239,679 → 110,197**（−54%）|
| **`HP_33` / `HP33`** | 3paths_audit / force_path2only 表格 column | 同樣是「reads with HP:i:3」計數 — **同一 metric！只是不同 BAM scope** | threshold_compare 全 BAM: 2,640 / 14,524 / 3,468；15-site cherry: 28 |
| **`HAPLOTYPE3` / phase block 3** | self-phasing 演算法層名稱 | longphase 內部對 phase block 3 的稱呼 | n/a |

**核心結論修正**（2026-05-07 釐清）：HP:i:33 與 HP_33/HP33 **是同一個 metric**（reads with HP:i:3），只是不同 audit 文件用不同命名 + 不同 BAM scope（whole genome vs threshold_compare 全 BAM vs 15-site cherry）。問題不是「兩個 metric」是「同一 metric 三個 scope 沒明示」。

### Why 鏈（5 層）

1. **Why PPTX 寫「14524 → 0」？** → master_draft 把 PI 報告 `HP:i:33 −54%` 跟 3paths_audit `HP_33 = 14524 → 0` 混為一談
2. **Why master_draft 會混？** → 兩個 metric 中文都簡稱「HP_33 / HP33 / HP:i:33」，名稱重疊
3. **Why 名稱會重疊？** → BAM 欄位 `HP:i:33` 是 SAM tag 格式，aggregate metric `HP_33` 是衍生計數變數 — 兩套命名不同層級但縮寫一致
4. **Why 沒在 master_draft / weekly-report 流程攔住？** → weekly-report W4 邏輯紅旗檢查只查「過度宣稱」「流水帳」，沒查「同名 metric 不同層」混用
5. **Why pptx-build 流程沒攔住？** → P3 section confirm 只看 focal point 是否服務 thesis，沒查每個數字 cite 是否 single-source。Agent-D Wave 2 才事後抓到，但連 Agent-D 自己也誤把 HP_33 當 HP:i:33

### Lesson（個人風格 candidate 規則）

> **Metric 引用必須 single-source verification**：每個數字進 PPTX 前必須 grep source 文件確認 metric 名稱、單位、層級、變化方向。Memory `feedback_feature_name_vs_definition_rule` 已強制 C++ feature；應擴展到 PPTX 數字引用。

---

## 1.2 P0 #2：noPath3 1.127 / HP_33=28 的 scope 錯誤

### 真相（來自 source）

force_path2only §3.4 標題明寫：「**BAM HP tag — 15-site cherry-picked sample**」
- noPath3：HP_33 = 28、HP1:HP2 = 1.127（**15 個 site 的小樣本**）
- §3.5 標「全 BAM 計數（待 monitor 完成）」，預期值 ≈ 0.735:1

### Why 鏈

1. **Why PPTX 把 1.127 / 28 呈現為與其他 4 版本同級的全 BAM 數字？** → 表格 column header 寫 HP_33（aggregate）但內容混了 15-site 與全 BAM
2. **Why 表格 mix scope？** → master_draft §0 Top Findings #4 與 Layer 0.1 把「15-site ≈ OLD V5 work tree」與「跨 23 染色體 17.3:1」並列寫成 evidence card
3. **Why 並列？** → 用戶當時想呈現「假說 PASS」結論，cherry-pick 是該假說 PASS 的最快證據；但沒標註 scope 差異
4. **Why 沒在 weekly-report W6 教授追問預測抓到？** → 預測題 7 個都是「結論層」追問，沒包含「scope coverage」追問
5. **Why pptx-build §20 雜訊紅旗沒抓到？** → §20.F 紅旗清單只有「過度宣稱」一條籠統項，沒細到「scope 標註缺失」

### Lesson（個人風格 candidate 規則）

> **任何數據引用須附 scope 標籤**：whole-genome / 15-site cherry / 1-sample / multi-sample 必明示在數字旁。

---

# 2. 受眾理解度評估（"相鄰專長" 標準）

## 2.1 受眾基線假設

依 `/big8_disk/liaoyoyo2001/Knowledge/README.md` 系統的「目標讀者：AI Agent 與新進研究人員」+ master_draft `audience: advisor`，本份簡報實際受眾是：

> **PI（領域專家但未親自操作程式碼）** — 來自 PI 報告 4/22 line 4 自我定義

也就是說：受眾**懂癌症基因組 + 長讀長 + ONT 大方向**，但**不熟 longphase-to 程式碼細節 / `HP:i:33` BAM tag / pphasing-graph 內部機制**。Knowledge base 為這類讀者寫，目標讀者為「AI Agent 與新進研究人員」。

## 2.2 18 張 slide 受眾理解度評估

| Slide | 主要術語 | 第一次看是否能懂？ | 待加強 |
|:-:|---|:-:|---|
| 01 Cover | V5 / 5 commits / Pass 1/2 / P0 | ⚠ Pass 1 vs Pass 2 需先解釋；P0 (priority 0) 是工程術語 | 加 1 行白話：「Pass 2 = second-round phasing；P0 = 最高優先」 |
| 02 Background | Thread D / Thread B / V5 audit chain | ❌ Thread D/B 是內部專案命名 PI 不一定記得 | 加 footnote「Thread D = LOH-constrained phasing 主軸」 |
| 03 Commit chain | SHA / cherry-pick / threshold | ⚠ cherry-pick 是 git 術語；threshold 0.95 vs 0.9 沒解釋指什麼 | 加 1 行：「threshold = highPurity 觸發 Pass 2 的純度門檻」 |
| 04 ★ Critical | ploidy bug / purity / highPurity / Pass 1 only | ❌ ploidy=0 → purity=0 → highPurity=false 因果鏈技術重 | 必補 1 行白話因果（slide 05 才有圖，slide 04 沒鋪墊） |
| 05 Mechanism | 5-stage flowchart | ✓ 5 box 圖示已視覺化 | OK |
| 06 ★ 三條路 | 路 1/2/3 / somaticCalling / second round phasing | ❌ 「somaticCalling 重跑」「只重 phase 不重 call」二者差別不解釋很難懂 | 補 1 行：「路 2 = 重新算誰是 somatic + 重新分配 haplotype；路 3 = 只重新分配，不重算」 |
| 07 5 版本表 | HP1:HP2 / HP_33 / 反轉 | ⚠ HP_33 沒定義 / 反轉指什麼沒解釋 | 補 1 行：「HP1:HP2 應約 1:1；偏離越遠 self-phasing 越嚴重；反轉=從偏 HP1 翻轉成偏 HP2」 |
| 08 Ablation | force / patch / 負控制 | ⚠ 「負控制實驗」是科學方法術語 PI 懂；hack 一行 patch 對 PI 是黑箱 | 不需大改，補 1 句即可 |
| 09 F1 | Caller F1 / FILTER / ISM SuggestFilter | ⚠ FILTER 欄位 vs F1 metric 兩者關係沒解釋 | 補：「F1 = caller 直接判定的精度；FILTER = VCF 中 PASS/LowQual 標記」 |
| 10 Caveat banner | 4 affected artifacts | ✓ 文字直白 | OK |
| 11 Decision tree | (a)(b)(c) 三選項 / 22.55 MB binary / upstream PR | ⚠ 「上游 PR」是 git 術語；22.55 MB 數字沒上下文 | 補：「上游 = longphase-to 原作者 zhenyu 的 repo」 |
| 12 Cross-sample | 7 個樣本表 / [U] | ✓ [U] 已加說明 | OK |
| 13 HPFineNGroups | marker / tier ⭐3 / phasing signature | ❌ marker / tier 升降是內部術語；phasing signature 是新詞 | 補：「marker = 候選生物標記；phasing signature = 因 phasing 演算法產生的副作用模式（非真實生物訊號）」 |
| 14 Thread B | Inter AF→NGroups / paired vs TO 層 | ❌ NGroups / paired vs TO 是內部術語 | 補：「paired = 含 normal 對照樣本；TO = tumor-only 無對照」 |
| 15 Future | P1-P5 / 25 hr 平行 | ✓ 工時 + 優先序直白 | OK |
| 16 Take-home | 3 件事 | ✓ 結構清晰 | OK |
| 17 Q&A | 必問 / 可能問 | ⚠ 整頁文字密 | 字級 + 排版改進（已在 issue list） |
| 18 References | 9 deliverable / 10 reference | ✓ 結尾頁 | OK |

**理解度評分**：可清楚理解 5 / 18（28%）；需補 1 行解釋 9 / 18（50%）；技術詞密度過高需重設計 4 / 18（22%）= **S04 / S06 / S13 / S14**

## 2.3 受眾必問 vs 我預測的差距

我在 master_draft §17 預測 7 個 Q&A，但**漏了至少 3 個受眾視角的基本問題**：

| 漏掉的問題 | 預期理由 |
|---|---|
| 「Pass 1 vs Pass 2 是什麼意思？」| Slide 04 直接用，沒鋪墊 |
| 「HP1:HP2 應該是多少？」| Slide 07 表格沒標 ground truth ratio（≈ 1:1）|
| 「為何不直接接受新 default？」| Slide 11 列三選項，但「為何維護 noPath3 變體」這個基本動機沒講 |

這三個問題 PI 在 30 sec 內就會問。當前 PPT 沒鋪墊。

---

# 3. 雜訊與必要性審計（去除 / 簡化清單）

## 3.1 應移除（不必要）

| # | 位置 | 原因 |
|:-:|---|---|
| N1 | S01 Cover「Report type: improvement_report (problem:progress 混合主線)」| 內部分類術語，PI 不需要看 |
| N2 | S15 「P1-P5 各對應 ⭐⭐⭐ / ⭐⭐ / ⭐ 標記」| 已用 P1-P5 + 紅 amber grey 配色傳達優先序，星號重複 |
| N3 | S03 timeline 中 commit SHA「8b8c1fd / 41ff147 / 380e8d2」| PI 不會記，可簡化為「Phase 0/1/2」+ 只保留 4/30 兩個 critical SHA |
| N4 | S06 三路 codeflow「(無 V5 flag)」「(PON 不重 phase)」括號內否定描述 | 干擾「描述路徑做什麼」主軸；改用肯定描述 |
| N5 | S07 全部行內「(bug)」「(forced)」標籤 | 表格已用紅綠 highlight 行，標籤多餘 |
| N6 | S17 Q4-Q7「⭐⭐」標記 | 此頁是 backup，重要性已隱含 |
| N7 | speaker note 多處 [ORAL-OPTIONAL] 提到 cpp 行號（如 PhasingProcess.cpp:142-220）| PI 不開 IDE，行號是給 follow-up audit 用的；移到 reference slide |

## 3.2 應簡化（必要但太重）

| # | 位置 | 簡化建議 |
|:-:|---|---|
| S1 | S04 「Evidence A：BAM provenance / Evidence B：機制因果鏈」雙欄 | 合併為 1 框「Evidence：BAM = 4/12 + log purity=0 + 無 'second round' 字串」 |
| S2 | S07 5 版本表 5 欄（版本/走的路/HP1:HP2/HP_33/反轉?）| 合併「走的路 + 反轉?」為 1 欄；HP_33 移 footnote |
| S3 | S14 Thread B「TO 層」「Paired 層」雙欄 | 合併為 1 段：「TO 層撤回，Paired 層保留 + Inter AF→NGroups +0.7」 |
| S4 | S18 Reference 10 項列表 | PI 不會即時讀；改 1 行「9 份 deliverable 路徑詳見 master_draft.md」 |

## 3.3 必加（受眾理解必需）

| # | 位置 | 該加 |
|:-:|---|---|
| A1 | S01 cover 下方 mini-glossary | 「Pass 2 = second-round phasing；P0 = 最高優先級」3 行 ≤ 30 字 |
| A2 | S04 因果鏈鋪墊 | slide 04 加 1 行白話：「ploidy 估值崩成 0 → purity 算出 0 → V5 跳過 second round」（slide 05 視覺化的文字版）|
| A3 | S06 路 2 vs 路 3 一行差異 | 「路 2 = 重算 somatic + 重分 hap；路 3 = 只重分 hap，不重算」（取代括號否定描述）|
| A4 | S07 ground truth ratio | 「HP1:HP2 應約 1:1；偏離 = self-phasing 嚴重」（標題副題）|
| A5 | S13 / S14 內部術語 footnote | marker / phasing signature / paired vs TO 一行解釋 |

---

# 4. 為何之前流程沒抓到 #1 / #2

## 4.1 weekly-report 流程的 gap

| 流程節點 | 當前查 | 沒查 | Gap |
|---|---|---|---|
| W4 邏輯紅旗 | 過度宣稱、流水帳 | metric 名稱單一性、scope 標註 | gap-1：metric ambiguity |
| W6 教授追問 | 結論層 5-7 題 | 受眾基線「Pass 1 是什麼」級的 0-th 問題 | gap-2：missing zero-th question |
| W7 4 桶分流 | tier S/A/B/C/D | 「相鄰專長受眾理解度」維度 | gap-3：no audience-comprehension axis |

## 4.2 pptx-build 流程的 gap

| 流程節點 | 當前查 | 沒查 | Gap |
|---|---|---|---|
| §20.E 6 問 audit | focal point / element ≤ 6 / 比例 / etc. | 每個數字是否 single-source verified | gap-4：no number-source verification |
| §20.F 雜訊紅旗 | 籠統「過度宣稱」 | scope label 缺失（whole vs cherry vs sample-1）| gap-5：scope label red flag missing |
| Wave 1 4 Agent | T/C/L/B（視覺）| 數字準確性（後來才有 Wave 2 D 抓到，但太晚）| gap-6：data accuracy 應在 Wave 1 跟視覺並行 |

## 4.3 個人風格 + 流程改進候選

依 feedback_classification 4 分類，建議升 active 規則的候選：

| 候選規則 | 通用 / 本次特定 |
|---|:-:|
| **R-G1**：metric 名稱混淆檢查（同名不同層 metric 必加 source path + scope 標籤）| **通用**（適用所有 PPT + master_draft）|
| **R-G2**：每張 slide ≤ 3 個 PI 不熟術語；多於 3 個須補 footnote | **通用** |
| **R-G3**：Wave 1 加 Agent-N（Number verification），與 T/C/L/B 並行 | **通用** |
| R-S1：HP:i:33 vs HP_33 vs HAPLOTYPE3 命名澄清條目 | 本次特定 |
| R-S2：self-phasing 主題 PPT 必含 ground-truth ratio (HP1:HP2 ≈ 1:1) | 本次特定 |

---

# 5. PPT 設計最佳實踐（網路研究整合）

## 5.1 Assertion-Evidence 結構（已採用，但可加強）

[Penn State Engineering / Alley 2013](https://writing.engr.psu.edu/ae_comprehension.pdf) + [iBiology](https://www.ibiology.org/professional-development/power-point-slide-design/)：
- 標題 = 完整 assertion 句（已做 ✓）
- Body = 1 主視覺支撐標題（已做 ✓ 但 S04/S07/S14 仍以文字為主，違反 "evidence = visual"）
- 字 ≤ 30 word per slide（PI 報告 60 字版較寬鬆，但 S02/S04/S08/S10 仍超 60 字）
- **2025 年 ScienceDirect study**：AE 結構降低 cognitive load + 提升 delayed recall（[Alzayed 2025](https://www.sciencedirect.com/science/article/pii/S2307187725001701)）

## 5.2 跨 expertise 受眾的雜訊降低（[Lucidchart](https://www.lucidchart.com/blog/how-to-explain-technical-ideas-to-a-non-technical-audience) / [Stanford Online](https://online.stanford.edu/10-tips-communicating-technical-ideas-non-technical-people)）：

關鍵 takeaway：
1. **避免行話 / acronym** — 必要時 1 行解釋於首次出現處
2. **「最重要的 10% 訊息重複 7 次」** — 本次 main thesis「V5 為 5 commits + Pass 1 only + Pass 2 重驗 P0」應在 S01 / S04 / S10 / S16 重複出現（目前 S01 + S04 + S10 + S16 已有但角度太分散；應有 1 致語）
3. **類比與譬喻** — Pass 1 vs Pass 2 可用「初稿 vs 校稿」類比；ploidy bug 可用「計算分母為 0」類比
4. **混合受眾平衡**：不簡化到專家覺得空，不複雜到非專家迷失 — 補 footnote 比改主敘述好

## 5.3 圖示繪製重點（為 #2 之後修正用）

依本研究，PPT figure 應「先告訴讀者要看什麼」：
- S03 timeline：✓ 已加紅色 highlight 4/30
- S05 因果鏈：⚠ 應在每個紅 box 加「結果」label（而非只有 metric 名稱）
- S06 三路 codeflow：⚠ 應在路 2 (綠) 加「success」、路 3 (紅) 加「fails」標題層 modifier
- S07 應改為 bar chart（替代純 table），讓「等價性」視覺化
- S15 priority gantt 比 5 cards 更傳達「P1+P2 平行、P3 後置」的時間關係

---

# 6. 修正計畫框架（待 ack 才動）

依 issue_list.md 既有 R1/R2/R3 + 本分析新增的 audit-driven 改進，整合為新 4 階段：

| Phase | 範圍 | 估時 | 必要性 |
|---|---|:-:|---|
| **Phase A — 事實層級修補（必修）** | P0 #1 (HP:i:33 真實值 239,679→110,197 + 加 metric scope 標) + P0 #2 (noPath3 加 15-site proxy 標) | 15 min | 不修則 PI 看到錯誤數字 |
| **Phase B — 受眾理解度補強** | A1-A5 必加項 + N1-N7 移除項 + S1-S4 簡化項 | 30 min | 不補則第一次聽 PI 抓不到主軸 |
| **Phase C — 視覺/字級/色彩** | issue_list R2-R3 全套（□ 方塊、字級 14pt 下限、紅綠雙編碼、focal chip 對齊）| 45 min | 視覺品質 |
| **Phase D — 流程改進規則寫入** | personal_style_log 加 R-G1 / R-G2 / R-G3（通用）+ R-S1 / R-S2（本次特定）| 15 min | 累積學習 |

**總估時：~1.75 hr**

# 7. 待你決策（3 個 confirm 點）

## Q1：修正範圍
- **建議：執行 Phase A + B + D（事實 + 受眾理解 + 規則寫入）**，估 1 hr
- 替代：Phase A only（極簡，10 min）/ A+B+C+D 完整（1.75 hr）

## Q2：HP:i:33 真實值寫法
我發現實際數字是 **239,679 → 110,197（−54%，4/29 技術報告 line 267）**，不是 master_draft 寫的「14524 → 0」也不是 Agent-D 推測的「0 → 14524」。
- **建議：改用 read-level「239,679 → 110,197」並標註「read 計數，PI 報告口徑」**
- 同時要在 S07 表格 caveat 加註：「HP_33（aggregate）≠ HP:i:33（read-level）兩者不同 metric」

## Q3：是否升 R-G1 / R-G2 / R-G3 為通用 active 規則？
- R-G1 metric 命名 single-source verification（**強烈建議升通用**，本次教訓深刻）
- R-G2 PI 不熟術語 ≤ 3 / slide（**建議升通用**）
- R-G3 Wave 1 加 Agent-N（建議升通用，需修 multi_agent_review.md）

請對 Q1-Q3 給答覆，我即按你的版本執行 Phase A → B → C（如選）→ D。

---

## Sources

- [Harvard Catalyst — Slides best practices](https://catalyst.harvard.edu/writing-communication-center/visualize-science/slides/)
- [Penn State / Alley — Assertion-Evidence comprehension PDF](https://writing.engr.psu.edu/ae_comprehension.pdf)
- [Alzayed 2025 — AE vs Traditional PPT cognitive load (ScienceDirect)](https://www.sciencedirect.com/science/article/pii/S2307187725001701)
- [iBiology — Rethinking Scientific Presentations](https://www.ibiology.org/professional-development/power-point-slide-design/)
- [assertion-evidence.com — Penn State Garner & Alley](https://www.assertion-evidence.com/)
- [Lucidchart — Explain Technical Ideas to Non-Technical Audience](https://www.lucidchart.com/blog/how-to-explain-technical-ideas-to-a-non-technical-audience)
- [Stanford Online — 10 Tips for Communicating Technical Ideas](https://online.stanford.edu/10-tips-communicating-technical-ideas-non-technical-people)
- [Mercedes Bernard — Explaining technical concepts](https://mercedesbernard.com/blog/how-to-talk-technical/)
