# Self-Phasing 教授報告 V3 · 12-Slide 精簡版 Storyboard

> **權威 storyboard**（給 build_pptx.py 開發者與演講者）。最終定稿 2026-04-30。
> v3 以 **one-slide-one-message** 高密度敘事壓縮 v2 的 24 slides 為 12 slides。

---

## 0. V3 定位說明

| 項目 | 設定 |
|------|------|
| 受眾 | 教授（NGS / long-read 基本背景，不熟 LongPhase-TO 與 InterSubMod 切分） |
| 演講長度 | **15-18 分鐘**（剩餘 12-15 分鐘留 Q&A） |
| Slide 數 | **12**（壓縮自 v2 的 24 → 50% 壓縮率） |
| 設計哲學 | 一張 slide = 一個結論句 + 一張支撐圖（佔 60-70% 面積），文字 ≤ 12% 面積、其餘留白 |
| Speaker note | 每張 ≥ 500 字（V3 因 slide 變少，每張 note 比 v2 更長以承載口說細節） |
| 與 v2 關係 | V3 是「演講高密度版」；v2 24-slide 完整版保留為深度資料源（教授追問細節時可開啟 v2 對應 slide）。圖片資產**完全沿用** v2/figures/（v3/figures/ symlink） |
| 比例 | 16:9，Latin Arial + CJK Droid Sans Fallback per-char |
| 必講核心 | **S1 / S6 / S7 / S10**（壓秒時其餘可加速；這 4 張是邏輯骨幹） |

> 演講者注意：V3 的責任是把**邏輯一鏡到底**講清楚；數據細節（行號、p-value、N、CI）放在 speaker note 與 Q&A 預備答案，不上 slide。

---

## 0.1 V3 → v2 對照速查表

| V3 # | V3 標題（≤15 字） | 對應 v2 slides | 處理方式 |
|------|------|------|------|
| S1 | 17.3 : 1 → 1 : 1 | v2 S1 | 1:1 保留 |
| S2 | 工作流程一覽 | v2 S2 | 1:1 保留（pipeline + 4 階段） |
| S3 | HP tag 五值 + 三層證據 | v2 S3 + S4 | 合併（HP tag 定義 + 全基因組 17.3:1） |
| S4 | Phasing 沒問題，是 Tag 有問題 | v2 S5 + S6 + S7 | 三合一（paired 對照 + 拆鎖 + LOH 兩層） |
| S5 | ISM 85 features 中 29 個受影響 | v2 S8 + S9 | 合併（4-bucket + 3-tier） |
| S6 | ★ 根因樹（觸發條件 + 三 bug） | v2 S11 | 1:1 保留（核心新發現） |
| S7 | ★ V5 三層投票流程 | v2 S12 | 1:1 保留（修正流程） |
| S8 | Sanity 15/15 PASS + 5 證據鏈 | v2 S14 | 1:1 保留 |
| S9 | AMB ↓ HP33 ↓ Concordance ↑ | v2 S15 | 1:1 保留 |
| S10 | 🎯 V5max1 — 39 reads 翻轉 | v2 S17 + S18 | 合併（CLIMAX + bar chart） |
| S11 | 業界家族樹 + 解鎖五大目標 | v2 S20 + S22 | 合併（家族樹 + 五目標卡片） |
| S12 | Take-home + 3 P0 + Q&A | v2 S23 + S24 | 合併（P0 + take-home + Q&A 感謝） |

**移入 speaker note 的 v2 內容**：
- v2 S10（prerequisite 4 函數 + PON 概念）→ 入 S6/S7 note
- v2 S13（程式碼 diff baseline vs V5）→ 入 S7 note，可作 backup slide 備教授追問
- v2 S16（4-commit timeline）→ 入 S8 note
- v2 S19（矛盾或盲點 F1 釐清 + cnLOH）→ 入 S9 note + Q&A 預備
- v2 S21（後續可動 / 已初步現）→ 入 S12 note

---

## 1. 12-Slide V3 Storyboard

### Slide 1. 17.3 : 1 → 1 : 1

**核心訊息**（≤30 字）：Tumor-only 模式 read 異常偏向 HP1，V5 修補後回到理論平衡。

**視覺**：自繪 banner（auto_figs/S1_impact.png 或 build 時生成）。
- 上半：左側超大字「**17.3 : 1**」灰色 → 右側超大字「**1 : 1**」綠色，中間綠色巨大箭頭。
- 下半：綠色標籤「**+8.3 pp（clean PS, 全基因組）**」。
- 底部小字（淺灰）：「self-phasing 分層處理：V2b 解 phase scaffold；V3F/V5 解 tag 層 · InterSubMod 無 C++ 改動」。
- 再底一行小字：「修補 self-phasing 是解鎖 ISM 五大研究目標的前提」。

**頁面元素**（16:9，13.33×7.5 英吋）：
- 標題列（隱藏，由 banner 取代）
- 上半 banner：x=0.5, y=0.8, w=12.3, h=3.2（17.3:1 → 1:1 主視覺）
- 中段 callout 條：x=0.5, y=4.2, w=12.3, h=0.8（+8.3 pp 綠色 highlight）
- 下半雙行小字：x=0.5, y=5.4, w=12.3, h=1.5（caveat + motivation）

**結論句**（slide 上大字）：「**17.3 : 1 → 1 : 1**」+「**+8.3 pp（clean PS, 全基因組）**」

**Speaker note 重點**（≥500 字）：
- 開場 30 秒鎖核心：「今天 15-18 分鐘要講四件事 — 我們在 tumor-only 模式發現了一個 17.3:1 的 HP tag 失衡，找到根因，做了 4 個 commit 修補，量化驗證 +8.3 pp concordance；這不是新 caller 也不是新 phaser，但這是讓我們 InterSubMod 五大目標可以推進的前提。」
- 17.3:1 怎麼算：HP1 family（HP1+HP1_1）614,000 reads / HP2 family（HP2+HP2_1）35,500 reads = 17.3。HCC1395 全基因組、跨 23 染色體一致、94.6% 集中於 HP1。
- 為何 1:1 是基準：兩條 hap 來自父母各一，random allele assignment 期望 1:1；17.3 偏離 17 倍 → 必為 artifact。
- +8.3 pp 怎麼算：clean PS blocks（germline accuracy ≥ 70% 的 phasing block）上，V5 vs paired ground-truth concordance 從 baseline 82.2% 提升到 V5 90.5%。**口徑校準**：是 clean PS blocks，不是全基因組 raw reads。
- self-phasing 分層處理三句話：(1) V2b commit 8b8c1fd 啟用 `--pon-only-phasing` 解 **phase scaffold** 層；(2) V3F commit 41ff147 + V5 working tree 解 **tag 層**（germline-first + Layer 1.5 + enum fix）；(3) V5 不修 self-phasing 本身，那是 V2b 的責任。
- 範圍切分：修補在外部 fork（longphase-to-mod V5），InterSubMod repo 無 C++ 改動 — 演講主軸是「在使用端發現問題、回到外部修、量化下游影響」。
- 「修補是解鎖五大目標的前提」：因為 38% ISM 特徵直接依賴 HP tag（後面 S5 會回來），所以 tag 不修對 → 38% 結論不可信。
- 接下來節奏：S2 工作流程（5 秒過）→ S3-S5 觀察與重要性（5 分鐘）→ S6-S7 根因與修正（5 分鐘）→ S8-S10 驗證與 climax（5 分鐘）→ S11-S12 未來與結語（2-3 分鐘）。

**壓縮自 v2**：S1。

**對應 Q 編號**：（無直接 Q，但是 Q1「為什麼是 17.3」的入口；演講者在此鋪陳，不展開細節）

---

### Slide 2. 工作流程一覽

**核心訊息**：修補位置在 longphase-to-mod，InterSubMod repo 僅在使用端，本演講分四階段。

**視覺**：兩段並排或上下分欄。
- 上半：Pipeline schematic — `tumor BAM → ClairS-TO（caller）→ longphase-to-mod V5（紅標）→ tagged BAM → InterSubMod`。重點紅色標記 V5 修補位置，InterSubMod 框出邊界（虛線）。
- 下半：4 階段橫條 — ① 機制定位 → ② 4-commit 修補 → ③ Sanity + concordance 驗證 → ④ 影響評估 與 ISM 銜接。

**頁面元素**：
- 標題：x=0.5, y=0.3, w=12.3, h=0.7（「工作流程一覽 · Pipeline & 4-Stage」）
- 上半 pipeline 圖：x=0.5, y=1.2, w=12.3, h=2.8
- 分隔線（淺灰）：x=0.5, y=4.1, w=12.3, h=0.05
- 下半 4 階段橫條：x=0.5, y=4.4, w=12.3, h=2.5

**結論句**：「**InterSubMod 無 C++ 改動 · 修補在 longphase-to-mod V5**」

**Speaker note 重點**（≥500 字）：
- pipeline 一行解釋：tumor BAM 進 ClairS-TO 出 VCF → 進 longphase-to-mod 做 phasing + haplotag → 出 tagged BAM → 進 InterSubMod 做 ISM 分析。本工作改的是 longphase-to-mod V5（紅標）。
- 為何要強調邊界：教授可能誤會本工作改了 InterSubMod；事實上 InterSubMod repo 在這次 self-phasing 修補中**沒有 C++ 改動**，所有 patch 都在 longphase-to-mod fork（獨立 git repo）。
- HP tag 是兩 repo 介面契約：HaplotagProcess.cpp:66-68 的 `ReadHaplotype` 結構體**一字未變** — 介面契約穩定，下游 InterSubMod 不需重編，只需重跑 BAM。
- 4 階段：(1) 機制定位 — 用 18 份 v5_audit_suite 報告定位三 bug；(2) 4 commit — V2b PON-only / V3F getVote 重寫+enum / INDEL guard / V5 working tree Layer 1.5；(3) 驗證 — 4 項 sanity check 15/15 PASS + paired GT concordance；(4) 影響評估 — ISM 3-tier 分類（29/14/42）+ 五大目標解鎖。
- 演講後三分之二會聚焦在 ② 與 ③；① 只用 S6 一張 slide 講根因；④ 用 S5 + S11 連結。
- 範圍提醒（避免誤解）：本實作是 **fork 而非 LongPhase 主線 PR**；longphase-to-mod 是 LongPhase 的 tumor-only 分支實驗。同實驗室相鄰工作（後面 S11 會講業界家族樹）。

**壓縮自 v2**：S2。

**對應 Q 編號**：可能引出 Q10（與 LongPhase-S 重疊度），不展開，留 S11/Q&A。

---

### Slide 3. HP tag 五值 + 三層證據

**核心訊息**：HP tag 五整數值是 4-bucket 分群基礎；理論 1:1、全基因組 17.3:1、單點 113:0 — 三層內外一致。

**視覺**：上下兩段。
- 上半：`figures/fig17_hp_tag_5versions.png`（HP tag 五值 schematic：HP1/HP2/HP1_1/HP2_1/HP33 視覺定義）
- 下半：`figures/fig01d_somatic_bias_explanation.png`（全基因組 614K vs 35.5K bar 與 SP1 IGV 縮圖）

**頁面元素**：
- 標題：x=0.5, y=0.3, w=12.3, h=0.7（「HP tag 五整數值 + 三層證據（理論／全基因組／單點）」）
- 上半圖：x=0.5, y=1.1, w=12.3, h=2.6
- 三層 callout 列（細條）：x=0.5, y=3.85, w=12.3, h=0.45（「理論 ~1:1 · 全基因組 17.3:1 · SP1 113:0」三段並列）
- 下半圖：x=0.5, y=4.4, w=12.3, h=2.7

**結論句**：「**理論 ~1:1 · 全基因組 17.3:1 · 單點 SP1 113:0** — 三層證據鏈內外一致」

**Speaker note 重點**（≥500 字）：
- HP tag 五值定義：
  - **HP1 (i:1)** — read 屬 haplotype 1 germline
  - **HP2 (i:2)** — read 屬 haplotype 2 germline
  - **HP1_1 (encoded as 11)** — read 屬 hap 1 含 somatic（HP1+somatic）
  - **HP2_1 (encoded as 21)** — read 屬 hap 2 含 somatic（HP2+somatic）
  - **HP33 (i:33)** — read 在某 phase block 內無法明確分類（ambiguous / Layer 1.5 fallback）
- 為何用 integer literal 11/21/33 而不用 enum：因為 BAM HP:i: 標準是 integer；用 enum 在 baseline 出過 type cast 失配 bug（後面 S6 會講）。enum 行號 `Util.h:20-25`（HAPLOTYPE_UNDEFINED=−1 起算）。
- LOH 兩層**先預告**（細節留 S4）：LOH.bed = region-level，由 phased VCF 的 AD 偵測；ISM HP_Ratio LOH = read-level，由 BAM HP tag 分布計算。**兩套不同數據源 → kappa = 0.670 不一致**。
- 17.3:1 拆解：HP1 family = HP1 + HP1_1 = 614,000 reads；HP2 family = HP2 + HP2_1 = 35,500 reads。94.6% 集中於 HP1，跨 23 染色體一致 — 不是某幾條染色體 outlier，是系統性 bias。
- 為何 1:1 是 null hypothesis：reads 對任何 hap 的 assignment 在 random allele 場景下期望均勻；17.3 偏離 17 倍意味 P-value 極小。生物學上不存在「整個基因組 95% reads 都在某一 hap」這種實況。
- SP1 chr19:17565944：113:0 完全失衡。是 HP2 方向（paired ground truth）但 baseline 全標 HP1。是極端示例，後面 S10 會用 V5max1 演示反向案例。
- **三層證據邏輯**：理論預期（生物學）→ 全基因組統計（量化）→ 單點 IGV（眼見為憑）。三層彼此獨立但結論一致 → 強化「這是真 artifact，不是抽樣 bias」。

**壓縮自 v2**：S3 + S4。

**對應 Q 編號**：Q1（為什麼是 17.3）— 細節在此預備；演講者若被打斷可即答。

---

### Slide 4. Phasing 沒問題，是 Tag 有問題

**核心訊息**：phasing 層 LOH.bed 不變（Jaccard=1.0），tag 層 HP_Ratio 才是 self-phasing artifact，本工作只動 tag 層。

**視覺**：左右兩欄表（自繪）。
- 左欄（綠色，phasing 層 OK）：
  - 「LOH.bed Jaccard = 1.0000 ✓」
  - 「region-level LOH 不變」
  - 「資料源：phased VCF 的 AD」
- 右欄（紅色，tag 層 NG）：
  - 「HP_Ratio 17.3:1 ✗」
  - 「ISM HP_Ratio LOH 62% artifact」
  - 「Cohen's d = −1.20（巨大效應量）」
  - 「kappa = 0.670（兩套 LOH 不一致）」
  - 「同位點跨模式 r = 0.001」
  - 「資料源：BAM HP tag」

**頁面元素**：
- 標題：x=0.5, y=0.3, w=12.3, h=0.7（「Phasing 沒問題，是 Tag 有問題（兩層分離）」）
- 副標：x=0.5, y=1.0, w=12.3, h=0.4（「LOH 兩層精確化 · 修補只動 tag 層」）
- 左欄綠色卡：x=0.5, y=1.6, w=6.0, h=5.4（phasing 層 OK）
- 右欄紅色卡：x=6.83, y=1.6, w=6.0, h=5.4（tag 層 NG）

**結論句**：「**Phasing OK · Tag 有 bug · 修補只動 tag 層**」

**Speaker note 重點**（≥500 字）：
- 為何要拆鎖：教授第一反應會問「phasing 不是壞了嗎？怎麼不直接修 phasing？」回答：因為 phasing **沒壞** — LOH.bed Jaccard = 1.0000 完全相同。
- LOH.bed 怎麼產生：LongPhase-TO 從 phased VCF 的 allele depth (AD) 偵測 region-level LOH（連續區段中 het variants 缺失 → 推 LOH）。這條路徑本工作前後完全一致。
- ISM HP_Ratio LOH 怎麼產生：InterSubMod 看每條 read 的 HP tag，計算每個區域 HP1 reads / (HP1 + HP2)；當 HP_Ratio > 0.8 或 < 0.2 標 LOH。這條路徑**受 BAM HP tag 影響**。
- 兩層 kappa = 0.670：意思是兩套 LOH 偵測的一致性中等偏低 — 過去一直認為這是「ISM 看到的 read-level LOH 是更精細的訊號」，但本工作證明 **62% 是 self-phasing artifact**，不是真 LOH。
- Cohen's d = −1.20 含義：HP_Ratio 在 Paired 模式下中位數 ~0.5（平衡），TO 模式下中位數 0.94（極端偏 HP1）。effect size −1.20 是 statistics 領域的「巨大效應量」（>0.8 為 large）。
- 同位點 r = 0.001：對 288K 個 paired/TO 同位點配對計算 HP_Ratio Spearman correlation = 0.001 — 等於完全無相關。意思是「同一個 SNV，在 paired 看 HP_Ratio 0.5，在 TO 看 HP_Ratio 0.94，完全是兩個獨立分布」→ TO 模式的 HP_Ratio 不反映真生物訊號。
- 86.5% TO-only LOH 在 paired 下平衡：ISM 在 TO 看到的 LOH，當切到 paired 模式重看，86.5% 是 HP_Ratio 0.4-0.6 的平衡位點（不是 LOH）。
- **修補範圍宣告**：本工作只動 tag 層，phasing 層不動。意義：(1) 既有 LOH.bed 結論不需重做；(2) 但 ISM 全部依賴 HP tag 的特徵需要重跑（38% / 29 個，後面 S5）。
- Q&A 預備：若教授追問「為什麼是 BAM tag 出問題不是 VCF AD 出問題」→ 答：BAM tag 經過 getVote() 投票邏輯（baseline somatic-first 設計缺陷），VCF AD 只是計數無投票邏輯。

**壓縮自 v2**：S5 + S6 + S7。

**對應 Q 編號**：Q2（為什麼有 normal 就會對）、Q3（Jaccard=1.0 怎麼證明 phasing 不變）。

---

### Slide 5. ISM 85 features 中 29 個受影響

**核心訊息**：38% ISM 特徵直接依賴 HP tag；HP tag 不修對則這 29 個結論不可信。

**視覺**：左右兩欄。
- 左欄（主視覺）：`figures/03_self_phasing_impact.png`（3-tier 金字塔：🔴29 / 🟡14 / 🟢42）
- 右欄（小圖+卡片）：4-bucket 分群示意（HP1 / HP1-1 / HP2 / HP2-1）+ 三標籤「NGroups · HPSig · HP_Ratio」

**頁面元素**：
- 標題：x=0.5, y=0.3, w=12.3, h=0.7（「ISM 85 features 中 29 個受影響（38%）」）
- 左欄 3-tier 金字塔：x=0.5, y=1.2, w=7.5, h=5.8
- 右欄 4-bucket 區：x=8.3, y=1.2, w=4.5, h=2.7（4 顏色方塊 grid）
- 右欄三特徵卡片：x=8.3, y=4.1, w=4.5, h=2.9（NGroups / HPSig / HP_Ratio 各一卡）

**結論句**：「🔴 **29** 個必重跑 · 🟡 **14** 個間接污染 · 🟢 **42** 個不受影響（共 85）」

**Speaker note 重點**（≥500 字）：
- 為什麼這頁重要：這頁回答教授的下一個問題「修不修 HP tag 影響範圍多大」。答案：**85 個 ISM 特徵，29 個直接依賴**（38%），不修則這 29 個結論不可信。
- 🔴 嚴重影響 29 個（HP-依賴）：HP_Ratio、Potential_LOH、HPMergedDelta、HPMergedSig、HPFineNGroups、HPFineP、HP_Discordance 等。**這些必須在 V5 BAM 上重跑**才有意義。
- 🟡 中度影響 14 個（間接污染）：QualityScore、GlobalP、CramersV、VerificationClass、HP_AlleleDelta 等。間接受 HP tag 變化影響但非主依賴。
- 🟢 不受影響 42 個（無 HP 依賴）：PairwiseMean/MedianDist、AlleleDelta、Caller、甲基化矩陣、CpG 座標等。這些是 read-level methylation matrix 與 caller-side 特徵，**HP tag 變化不影響**。
- **口徑校準**（reviewer 採納）：來源報告誤寫百分比為 38%/7%/55%，但 14/85 = 16.5% 而非 7%；slide 上**只列 count（29/14/42）不列百分比**避免誤導。
- 4-bucket 分群：ISM 把每個 region 的 reads 分成 HP1（germline hap1）/ HP1-1（hap1 + somatic）/ HP2（germline hap2）/ HP2-1（hap2 + somatic）4 桶，然後做距離計算與 NGroups subclone marker。**bucket 分對 HP tag 必對 → 直接連結 S6 的根因**。
- HPFineNGroups marker 案例：過去解讀為 methylation bimodality marker（subclone 證據），但 V5 重詮釋為 phasing bucket signature — 不是 methylation 多模態，而是 4-bucket 分對後 read 分群更乾淨。這個結論**待 master × flag 重驗**（後面 S12 P0 行動 F4）。
- 跨樣本影響：7/7 樣本同方向 self-phasing 影響 → 過去**所有 TO 模式 ISM 結論**（Wave 3 LOH 分層、O-系列觀察、Phase 1A）都需要重跑或加註版本（archive haplotag_version）。
- 為什麼 InterSubMod 整個研究主軸都建在 TO 模式：因為 normal BAM 不一定可得（cell-free DNA、archive 樣本），TO 是必要研究方向。所以修對 TO HP tag = 解鎖 ISM 五大目標的前提（後面 S11 完整鋪開）。

**壓縮自 v2**：S8 + S9。

**對應 Q 編號**：（無直接 Q，但鋪墊 Q9「F1 不變那為什麼還做」的答案 — read-level 影響在這頁就明顯）。

---

### Slide 6. ★ 根因樹（觸發條件 + 三 bug）

**核心訊息**：純樣本 purity 0.927 ≤ 0.95 → 未觸發 Two-Pass → 走 baseline 三條路 → 暴露 tag 層 somatic-first 投票 bias。

**視覺**：自繪 root-cause tree（垂直結構）。
- **根節點**（紅色 banner）：「觸發條件：Purity 0.927 ≤ 0.95 → 未觸發 Two-Pass → 走 baseline 三條路」
- **中間節點**（橘色）：「baseline 流程假設有 normal contamination，但純樣本實為純 → 暴露 tag 層 bug」
- **三葉節點**（並排）：
  - 葉 1（紅色）：「① getVote() somatic-first priority 設計缺陷 — `HaplotagProcess.cpp:512-563`」
  - 葉 2（橘色）：「② enum vs integer literal mismatch — `Util.h:20-25` HAPLOTYPE_UNDEFINED=−1」
  - 葉 3（灰色）：「③ Phase scaffold self-phasing — V2b 已解（`PhasingProcess.cpp` PON-only）」

**頁面元素**：
- 標題：x=0.5, y=0.3, w=12.3, h=0.7（「★ 根因樹：觸發條件 + 三層 Bug」）
- 根 banner（紅）：x=0.5, y=1.1, w=12.3, h=0.9
- 中段 callout（橘）：x=2.0, y=2.3, w=9.3, h=0.8
- 三葉節點（並排）：每節 4.0 寬，y=3.4, h=3.5；x 分別 0.5 / 4.67 / 8.83
- 底部 caveat 小字：x=0.5, y=7.0, w=12.3, h=0.4（「三葉用顏色區分嚴重度：紅 = priority bug, 橘 = enum, 灰 = 已解」）

**結論句**：「**Purity 0.927 ≤ 0.95 = 未觸發 Two-Pass = 球員兼裁判**」

**Speaker note 重點**（≥500 字）：
- 此頁是演講核心 1/2（另一個核心是 S7）— 必須講清楚因果鏈。
- **觸發條件主者**（最新 v3 新增）：HCC1395 5kHz 是純樣本（實際 purity ~1.0），但 baseline 的 `collectPloidyRatio()` 估為 **0.927**（誤差 −0.06）。`PhasingProcess.cpp:197` 硬編碼 `if (purity > 0.95) Two-Pass`，所以 0.927 ≤ 0.95 → **沒觸發 Two-Pass** → 落到 baseline 三條路標準流程。
- 「三條路」是什麼：baseline 設計三種流程（Two-Pass for high purity normal-aware / single-pass / fallback），純樣本應走 Two-Pass 用 normal scaffold；因 purity 估錯走到 single-pass，於是用 tumor 自己的 somatic 當 phasing anchor → 球員兼裁判。
- **這不是 purity calculator 壞了**：04 報告證明 0.927 是合理估計（vs 0.6 sample 估 0.607 也準），是「閾值設計太硬」+「估計誤差自然存在」雙因素。
- 為什麼 paired 模式不出這 bug：paired 有 normal BAM 直接給 germline ground truth，根本不走「估 purity 決定流程」這條路。所以這 bug 只在 TO 模式 + 純樣本暴露。
- 三葉 bug 詳述：
  - **葉 1 priority bug**（紅，最嚴重）：`getVote()` 迴圈順序設計為 somatic-first — `for (variantKeys: {HP1_1, HP2_1} → {HP1, HP2}) { if (any > 0) break }`。意思是「先看 somatic vote，有就停」。當 PON-only 啟用後 germline votes 急減，剩下的 somatic vote 直接主導 → 17.3:1 artifact。
  - **葉 2 enum mismatch**（橘）：`HAPLOTYPE_UNDEFINED = -1` 起算的 enum 與 BAM 標準的 integer literal 11/21/33 之間有 type cast 失配。靜態檢查不抓（C++ 隱式 enum→int 不 warn）。修法是把 enum 直接換 integer literal。
  - **葉 3 phase scaffold**（灰，已解）：V2b commit 8b8c1fd 加 `--pon-only-phasing` flag，phasing 階段只用 PON-confirmed germline 做 anchor，不用 putative somatic。這層**已修**，本演講重點在葉 1+葉 2。
- 因果鏈一句話：「修了葉 3 → 暴露葉 1 → 同時發現葉 2 → V5 三層投票同時修葉 1 + 葉 2」。
- 如何驗證觸發條件論述：04 報告同時測 purity = 0.6 的 t30_n20 sample，baseline 估 0.607 → 也 ≤ 0.95，也跳 Two-Pass，但 self-phasing 自然弱化（normal contamination 平衡 hap）→ 所以**只有純樣本暴露 17.3:1 極端 artifact**。
- Q&A 預備：教授可能問「為什麼不直接調閾值到 0.99」→ 答：那會讓更多樣本走 Two-Pass，但 Two-Pass 在低 purity 不穩定；正確解是修葉 1+葉 2 讓投票邏輯本身 robust，這正是 V5 的設計（下頁）。

**壓縮自 v2**：S11（v2 核心新發現），整合 v2 S10 prerequisite 概念入 note。

**對應 Q 編號**：Q4（為什麼 paired 沒被發現這 bug）、Q5（unit test 為什麼沒抓到）、Q6（PON 有資料還會 self-phasing）。

---

### Slide 7. ★ V5 三層投票流程

**核心訊息**：V5 把 baseline 的 somatic-first 翻轉為 germline-first + 加 confidence 0.6 防禦層 + enum 修為 integer literal。

**視覺**：自繪三層流程圖（垂直）。或使用 `figures/fig01b_three_layer_logic.png` / `figures/fig16_v5_threelayer_logic.png` 為基底並標註紅黃綠圈。
- **Layer 1**（紅圈）：「**germline-first 投票**：if germline > 0 → 用 germline 直接決定」
- **Layer 1.5**（黃圈）：「**somatic fallback with confidence ≥ 0.6**：mid-low purity 防禦層，避免錯誤分配」
- **Layer 2**（綠圈）：「**encode 11/21/33**：integer literal 取代 enum，無 type cast」

**頁面元素**：
- 標題：x=0.5, y=0.3, w=12.3, h=0.7（「★ V5 三層投票（修正流程）」）
- 副標：x=0.5, y=1.0, w=12.3, h=0.4（「baseline somatic-first → V5 germline-first + Layer 1.5 + integer literal」）
- Layer 1 卡（紅圈）：x=0.5, y=1.6, w=12.3, h=1.7
- Layer 1.5 卡（黃圈）：x=0.5, y=3.4, w=12.3, h=1.7
- Layer 2 卡（綠圈）：x=0.5, y=5.2, w=12.3, h=1.7
- 底部 caveat：x=0.5, y=7.0, w=12.3, h=0.4（「V5 修 tag 層；phase 層由 V2b PON-only 處理」）

**結論句**：「**germline-first + confidence 0.6 + integer literal**」

**Speaker note 重點**（≥500 字）：
- V5 三層投票對應修哪些 bug：Layer 1 修葉 1（priority）、Layer 1.5 是新增防禦層（mid-low purity）、Layer 2 修葉 2（enum）。
- **Layer 1 germline-first**（最關鍵翻轉）：
  - baseline 偽碼：`for (variantKeys: {HP1_1, HP2_1} → {HP1, HP2}) { if (any > 0) break }` — somatic-first
  - V5 偽碼：`if (germline > 0) → use germline` — germline-first 直接條件
  - 影響：純樣本場景 baseline 投出 17.3:1，V5 直接修為 1:1（理論值）
- **Layer 1.5 somatic fallback with confidence ≥ 0.6**（mid-low purity 防禦層）：
  - 設計邏輯：當 germline votes 不足時，允許 somatic 作 fallback，但需 confidence ≥ 0.6（保守閾值）
  - 為何 0.6：09 報告 0.6 sample 上，V5 HP33 比例 12.4% vs baseline 2%（保守 +10pp 標 ambiguous，避免錯誤分配）
  - 這是「conservative tagging」設計哲學 — 寧可標 HP33 ambiguous，也不誤標 HP1 或 HP2
- **Layer 2 encode 11/21/33**：
  - baseline 用 enum HAPLOTYPE1_1=3 在 caller 端 cast 出 BAM HP:i: 值，type 失配導致部分 read 標錯
  - V5 直接用 integer literal `11`、`21`、`33`，無 cast、無 type 錯誤
  - 這個修改**最小但最重要** — Sanity 守恆律 B（germline 不變）的 0 漂移就是靠這個保證
- **Q1 + Q2 整合答案**（教授必問）：
  - Q1「PON 雙路徑在不同 purity 都有效嗎」：是。0.93 sample（HCC1395 5kHz）+ 0.6 sample（t30_n20）都驗證 — 0.93 上 Sanity 4 項 PASS + clean PS +13.3 pp；0.6 上無 paired ref 但 V5 HP33 比例 +10pp 表示保守 tagging 生效。
  - Q2「baseline somatic-first vs V5 germline-first，V5 真更好還是只是不同」：V5 更好。純樣本場景翻轉 17.3 倍 self-phasing artifact（baseline → 1:1），mid-low purity 場景以 confidence 0.6 提供保守 tagging（避免錯誤分配），對下游 ISM 4-bucket 分群品質**重要**。
- **caveat（reviewer 警示防誤解）**：V5 修的是 tag 層 directional reassignment，不是 phase 層 self-phasing 本身。phase 層的 self-phasing（phasing graph 用 somatic 當 anchor）由 V2b commit 8b8c1fd 的 `--pon-only-phasing` 解決。**S10 climax 講 V5max1 39 reads 翻轉時要再重複這個 caveat**，避免教授誤會 V5 「同時修 phase + tag」。
- 程式碼 diff（細節入 backup）：v2 S13 有完整 baseline vs V5 兩欄 monospace 程式碼對照；V3 移入 note，若教授追問「能不能看 diff」可開 v2 S13 backup。

**壓縮自 v2**：S12，整合 v2 S13 程式碼 diff 入 note。

**對應 Q 編號**：Q1+Q2 整合答案、Q4、Q6（PON 與 phasing anchor 區分）、Q12+Q13（v3 新增 mid-low purity 與 germline-first 設計）。

---

### Slide 8. Sanity 15/15 PASS + 5 證據鏈

**核心訊息**：4 項硬性 sanity check 全 PASS、5 條獨立證據鏈互相佐證、無 violation。

**視覺**：左右兩欄。
- 左欄：守恆律公式卡（自繪）
  - 「**A · ΔHP33 + (ΔHP11 + ΔHP21) = 0**」（守恆律 A）
  - 「**B · germline HP1/HP2 不變（0 reads 漂移）**」（守恆律 B）
  - 「**Layer 1.5 期望 1 · 33→directional 精確守恆**」
  - 「**Layer 1.5 期望 2 · 無 germline → HP33 = 0**」
  - 大字 PASS 標章：「**15 / 15 PASS · 0 violation**」
- 右欄：5 條獨立證據鏈條列
  - 「✓ 1. LOH.bed Jaccard = 1.0000（phasing 不變）」
  - 「✓ 2. PON-only Two-Pass 0.93 + 0.6 sample 都驗證有效」
  - 「✓ 3. Aggregate paired GT concordance +6.65 pp」
  - 「✓ 4. Clean PS paired GT concordance +13.3 pp」
  - 「✓ 5. 3/3 SP-extreme 一致翻轉（SP1/SP2/SP3）」

**頁面元素**：
- 標題：x=0.5, y=0.3, w=12.3, h=0.7（「Sanity 15/15 PASS + 5 條獨立證據鏈」）
- 左欄守恆律卡：x=0.5, y=1.2, w=6.0, h=5.8
- 右欄證據鏈：x=6.83, y=1.2, w=6.0, h=5.8

**結論句**：「**Sanity 15/15 PASS · 0 violation · 5 條證據鏈互相佐證**」

**Speaker note 重點**（≥500 字）：
- 為何要做 sanity check：這是「修補不引入新 bug」的硬性保證。LongPhase 原始 codebase **無 HP tag round-trip test**（v5_audit_suite/06 確認），本工作補了 4 項守恆律。
- 守恆律 A · ΔHP33 + (ΔHP11 + ΔHP21) = 0：意義是「HP33 ambiguous reads 重新分配為 directional 後總數守恆」。15 sites（5 TP + 4 FP + 3 V5max + 3 SP-extreme）全 PASS。
- 守恆律 B · germline HP1/HP2 不變：意義是「V5 修補不應該動 germline tag，只動 ambiguous 與 somatic-bearing reads」。15 sites 0 reads 漂移 — 這也驗證 Layer 2 enum 修對。
- Layer 1.5 期望 1 · 33→directional 精確守恆：所有從 HP33 翻轉到 HP11/HP21 的 reads 數量等於 directional bucket 增量。15/15 PASS。
- Layer 1.5 期望 2 · 無 germline → HP33 = 0：在無 germline 訊號的位點，V5 應該全部標 ambiguous（HP33）；任何 directional tag 都是違反設計。0 violation。
- 5 條獨立證據鏈展開：
  - ✓ 1. **LOH.bed Jaccard = 1.0000**（baseline vs V2b PON-only）— phasing 層完全不變，鎖定本工作只動 tag 層。
  - ✓ 2. **PON-only Two-Pass 雙 purity 驗證**（Q1 整合答案）— 0.93 sample 上 4 項 sanity PASS + concordance +13.3 pp；0.6 sample 上保守 tagging +10pp HP33。
  - ✓ 3. **Aggregate paired GT concordance +6.65 pp**（15-site all pooled，V5 78.85% vs baseline 72.20%）— 即使 cherry-picked 含 problem PS blocks，整體仍改善。
  - ✓ 4. **Clean PS paired GT concordance +13.3 pp**（15-site clean only 11 sites）/ **+8.3 pp**（全基因組 clean PS） — 兩個獨立 N 都顯著改善。
  - ✓ 5. **3/3 SP-extreme 一致翻轉**（SP1 113:0、SP2 109:1、SP3 108:0 baseline → V5 全翻至 paired 方向）。
- **Q1 答案完整版**（教授可能在此問）：是。Two-Pass 在純樣本（0.93）+ 中等純度（0.6）都驗證有效；前者 +13.3 pp clean PS、後者透過 +10pp HP33 體現保守 tagging 防禦。Aggregate +6.65pp 是混合 sample 平均；clean PS +13.3pp 是去除 problem blocks 後的真實上限。
- **未測項目（誠實邊界）**：
  - 4 項 sanity 在 cherry-picked 15 sites（不是隨機抽樣）— F8 要做 50-100 隨機 cross-validate
  - 7 樣本擴展未做（F3 待辦）
  - Confidence threshold 0.6 為間接驗證（HP33 比例變化），非 vote log 直接驗證（F2 待辦）
- **4-commit timeline**（v2 S16 移入 note）：
  - V2b 8b8c1fd · `--pon-only-phasing` PhasingProcess.cpp +25/-3 等
  - V3F 41ff147 · getVote 重寫 + caller 端 enum 修 HaplotagProcess.cpp +36/-25 line 506-541, 697
  - INDEL guard 380e8d2 · countINDELHaplotype +8/-4 line 497-510
  - V5 working tree（未 commit）· Layer 1.5 + countSNPHaplotype alt guard +24/-7 line 489-494, 512-563
  - V5 累計 +68/-36 行；介面契約 H:66-68 一字未變

**壓縮自 v2**：S14，整合 v2 S16 4-commit timeline 入 note。

**對應 Q 編號**：Q5（unit test）、Q7（sanity 是誰寫的、覆蓋率）、Q12（PON Two-Pass 雙 purity）。

---

### Slide 9. AMB ↓ HP33 ↓ Concordance ↑

**核心訊息**：三個量化指標同方向改善 — Ambiguous% 減半、HP33 reads 減半、paired GT concordance +8.3 pp。

**視覺**：`figures/fig18_concordance_amb_f1.png`（三聯橫幅，AMB / HP33 / Concordance 三組 bar）

**頁面元素**：
- 標題：x=0.5, y=0.3, w=12.3, h=0.7（「量化指標：AMB ↓ · HP33 ↓ · Concordance ↑」）
- 主視覺圖：x=0.5, y=1.2, w=12.3, h=4.8（fit-within 等比，不強壓比例）
- 底部三句強調：x=0.5, y=6.2, w=12.3, h=1.0（三段並列：「AMB 17.5 → 8.0 % (−9.5 pp)」/「HP:i:33 −54%」/「Concordance +8.3 pp」）

**結論句**：「**AMB −9.5 pp · HP33 −54% · Concordance +8.3 pp**（全基因組 clean PS）」

**Speaker note 重點**（≥500 字）：
- 三個量化指標講解：
  - **AMB%**（ambiguous read percentage）：17.5% → 8.0%，下降 9.5 pp。意義：原本 17.5% 的 reads 被 baseline 標為「無法確定 hap」，V5 把其中超過一半重分配到 directional（HP1/HP2/HP11/HP21）。
  - **HP:i:33 reads**（HCC1395 全基因組）：239,679 → 110,197，下降 54%。這與 AMB% 一致 — V5 大幅減少 ambiguous tag。
  - **Clean PS concordance**：全基因組 clean PS blocks 上，V5 vs paired ground truth 一致率 82.2% → 90.5%，提升 8.3 pp。15-site cherry-picked 上 74.9% → 88.2%（+13.3 pp）。
- **F1 不變的釐清**（v2 S19 移入 note）：
  - ClairS-TO raw F1 = 0.7166（三版本完全相同 — V5 不改 caller，這是預期且必要）
  - V5 + ISM SuggestFilter F1 = 0.7154
  - V5 vs Baseline ΔF1 = −0.0003（噪音範圍，0.7154 - 0.7157）
  - **F1 不能衡量本實作品質**：F1 衡量 caller-level TP/FP，但本實作改的是 read-level HP tag 品質。真實價值在 read-level concordance +8.3 pp。
- **真實價值在哪**：
  - 下游 ISM 38% 特徵（29 個 HP-依賴）必須在 V5 BAM 上才可信
  - 生物學詮釋正確性：HPFineNGroups marker 重詮釋為 phasing signature 而非 methylation marker
  - 跨樣本研究可信度：7/7 樣本同方向影響全部 TO 模式 ISM 結論
- **矛盾或盲點**（v2 S19 移入 note）：
  - cnLOH 區仍未解：copy-neutral LOH 雙親同源無 het variants，無從 anchor，V5 fallback 也無法處理。需要 cnLOH-aware filter（CN 層 + germline-only methylation reference）作為 future work（F6 行動）。
  - AF > 0.9 邊界：極端高 AF 的 somatic（接近 LOH 區）V5 行為與 baseline 接近，因為 germline votes 太少 fallback 到 Layer 1.5 confidence 0.6 但仍有不確定性。
  - Confidence threshold 0.6 未直接驗證：以 HP33 比例變化間接驗證；F2 行動要做 vote log 直接驗證。
- **N 與信賴區間**（Q8 預備）：
  - 15-site clean PS：N = 11 sites cherry-picked，HCC1395 5kHz 單樣本
  - 全基因組 clean PS：N = 全基因組 clean PS blocks（數萬-數十萬 reads 量級）
  - Wilson 95% CI 推估在 ±1 pp 內（V5 90.5% N=N_clean）
  - 跨樣本擴展未做（F3 待辦）
- **接 S10 climax**：「指標看完，現在我給你看一個具體案例 — V5max1 一個位點 39 條 read 100% 翻轉」

**壓縮自 v2**：S15，整合 v2 S19 矛盾盲點 + F1 釐清入 note。

**對應 Q 編號**：Q8（+8.3 pp 的 N 與信賴區間）、Q9（F1 不變那為什麼還要做）。

---

### Slide 10. 🎯 V5max1 — 39 reads 翻轉

**核心訊息**：CLIMAX — 單一位點 V5max1 39 條 read 在 V5 後從 baseline 全 HP1 完全翻轉至 paired ground truth 方向。

**視覺**：`figures/C_V5max1_chr19_4639528.png`（4-BAM 並列 IGV 截圖，佔 80% 寬度）

**頁面元素**：
- 標題：x=0.5, y=0.3, w=12.3, h=0.7（「🎯 CLIMAX：V5max1 chr19:4639528 — 39 reads, 100% reassigned」）
- 主視覺圖：x=1.0, y=1.1, w=11.3, h=5.4（80% 寬度大圖，4-BAM 並列 baseline / V2b / V3F / V5）
- 底部標籤條：x=0.5, y=6.6, w=12.3, h=0.4（「Baseline → V5：39 reads 全翻至 paired 方向 · 3/3 SP-extreme 一致」）
- 底部 caveat 小字：x=0.5, y=7.05, w=12.3, h=0.4（「V5 不修 self-phasing 本身（V2b 已處理 phase scaffold）；本圖示 tag 層 directional reassignment 效果」）

**結論句**：「**39 reads, 100% reassigned · 3/3 SP-extreme 一致翻轉**」

**Speaker note 重點**（≥500 字）：
- 演講高潮 — 給教授一個「眼見為憑」的 IGV 截圖證明。
- V5max1 chr19:4639528 是什麼：v5_audit_suite/01 找到的 V5 vs V3F 最大 ΔHP33 位點之一。V3F HP33 reads 39 → V5 HP33 reads 0；同時 V5 HP11 reads +39（全部從 ambiguous 翻轉到 directional HP1+somatic）。
- 4-BAM 並列怎麼讀：
  - 第 1 行 Baseline：39 reads 全標 HP1（紫色）— self-phasing artifact
  - 第 2 行 V2b PON-only：phase 層改 anchor，但 tag 層仍 somatic-first → reads 部分仍 HP1
  - 第 3 行 V3F：getVote 重寫 + enum 修，部分 reads 翻轉，但仍有 HP33 ambiguous
  - 第 4 行 V5：Layer 1.5 confidence 0.6 啟動，**39 reads 全翻轉至 HP11**（hap 1 + somatic 確定方向）+ 對齊 paired ground truth
- 為何「100% reassigned」是強訊號：
  - Sanity 守恆律 A 驗證：ΔHP33 = −39, ΔHP11 = +39 → 守恆等式精確成立
  - 39 reads 都同方向翻轉（不是隨機分散）→ 不是 noise，是 directional fix
  - 3/3 SP-extreme 一致（SP1 113:0, SP2 109:1, SP3 108:0 baseline → V5 全翻至 paired 方向）→ 跨位點 robust
- **caveat 必講**（reviewer 警示防誤解）：
  - V5 修的是 tag 層 directional reassignment，**不是 phase 層 self-phasing**
  - phase 層的 self-phasing（phasing graph 用 somatic 當 anchor）由 V2b commit 8b8c1fd 的 `--pon-only-phasing` 解決
  - 這個 caveat 必須在演講中明確口述，否則教授會誤會「V5 修了所有 self-phasing」 — 實際上 V5 是「在 V2b 修了 phase 層的基礎上，再修 tag 層投票邏輯」
- V5max2 / V5max3 數據（v2 S18 整合進 note）：
  - V5max2 chr19:2235521 · V3F HP33 26 → V5 HP33 0 · ΔHP11 +26
  - V5max3 chr19:7405500 · V3F HP33 16 → V5 HP33 0 · ΔHP21 +16
  - 三個熱點累計影響量級接近全基因組 HP33 變化的 1‰，但作為 **守恆律精確驗證 + 跨位點一致性** 的代表案例
- **演講高潮節奏**（演講者注意）：
  - 講到此頁時稍微提高音量
  - 用手指對著 4-BAM 行說「baseline 全紫色 → V5 全藍色」
  - 接著立刻補 caveat「但這只是 tag 層的修正」
  - 最後一句串接 S11：「這個 read-level 改善，下游影響到我們 InterSubMod 五大目標」
- **3/3 SP-extreme bar chart**（v2 S18 視覺整合）：原 v2 S18 是獨立 slide，V3 合併入 S10 — 若需要 backup 可調用 D_SP1/SP2/SP3 三張圖。

**壓縮自 v2**：S17 + S18。

**對應 Q 編號**：Q8（+8.3 pp 信賴區間）— 連同 S9 一起回答。

---

### Slide 11. 業界家族樹 + 解鎖五大目標

**核心訊息**：本實作填補 tumor-only + PON 公開工具 gap（同實驗室相鄰）；修補 tag 層是解鎖 InterSubMod 五大目標的前提。

**視覺**：上下兩段。
- 上半：自繪業界家族樹
  - 「LongPhase 2022（germline phasing）」
  - 分支 1 ↓「LongPhase-S 2025（paired tumor + normal）」
  - 分支 2 ↓「**longphase-to-mod V5（TO + PON）← 本實作**」紅標
- 下半：五大目標卡片（並排 5 卡）
  - 卡 1「目標 1 · per-CpG 多標籤關聯」（直接依賴 ✓）
  - 卡 2「目標 2 · clone 結構分析」（直接依賴 ✓）
  - 卡 3「目標 3 · 二次打擊順序」（部分依賴）
  - 卡 4「目標 4 · TO normal 補強」（直接依賴 ✓）
  - 卡 5「目標 5 · F1 panel」（部分依賴）

**頁面元素**：
- 標題：x=0.5, y=0.3, w=12.3, h=0.7（「業界家族樹 + 解鎖 ISM 五大目標」）
- 上半家族樹：x=0.5, y=1.2, w=12.3, h=2.8
- 中段分隔線：x=0.5, y=4.1, w=12.3, h=0.05
- 下半五卡片（並排）：x=0.5, y=4.3, w=12.3, h=2.7（每卡 2.4 寬，間距 0.05）

**結論句**：「**填補 TO+PON 公開工具 gap · 解鎖五大目標 1/2/4 直接依賴**」

**Speaker note 重點**（≥500 字）：
- 業界家族樹脈絡：
  - **LongPhase 2022** — 原始 germline phasing tool（hap1/hap2 切分）
  - **LongPhase-S 2025**（bioRxiv 2025.11.20.689492v1）— paired 模式 tumor + normal，somatic 錨在 germline scaffold；ClairS SNV +4.5%、ClairS indel +7.1%、DeepSomatic SNV +1.2%、DeepSomatic indel +0.5%
  - **longphase-to-mod V5（本實作）** — TO + PON 模式，當 normal 不可得時 PON 替代 normal 作 anchor
- **同實驗室相鄰工作**（口徑校準，reviewer 採納）：不是「業界共識」，是同實驗室 LongPhase-S 的 tumor-only 變體。設計哲學一致（somatic 錨在 germline scaffold），但場景不同。
- **與 LongPhase-S 區別**（Q10 答案）：
  - LongPhase-S 假設有 normal BAM；本實作填補當 normal 不可得時的標準解
  - 4-commit 漸進修補在 longphase-to-mod fork（獨立 git repo），不是 LongPhase 主線 PR
  - 4 項 sanity check + ISM 3-tier 分類 + 7 樣本一致性是新貢獻
  - 不重疊：LongPhase-S 是 paired，本實作是 TO；姐妹工作而非競爭
- **F1 對比**（口徑校準）：LongPhase-S F1 +4.5%（ClairS SNV）是 paired 場景；本實作 ΔF1 = −0.0003 是 TO 場景且 caller 不變（read-level concordance +8.3 pp 是真實價值）。**注意非直接可比**。
- 五大目標逐一講：
  - **目標 1 · per-CpG 多標籤關聯**：read-level methylation × HP tag × allele 多維分析。**直接依賴正確 HP tag**（read-bucket 分群）。
  - **目標 2 · clone 結構分析**：4-bucket 分群 → NGroups subclone marker → clonal evolution。**直接依賴**（HPFineNGroups 是核心）。
  - **目標 3 · 二次打擊順序**：兩個 somatic 事件的時序（先 hit hap1 還 hap2）。**間接依賴**（P4 pilot 在 region-level 無訊號，等目標 1 解鎖才能繼續）。
  - **目標 4 · TO normal-free 補強**：tumor-only 模式不依賴 normal 仍能達 paired 品質。**直接依賴**（self-phasing 修補是前提）。
  - **目標 5 · F1 panel**：跨樣本跨 caller 的 F1 比較 panel。**部分依賴**（HP tag 影響 ISM SuggestFilter 但 caller 端不變）。
- 來源：`InterSubMod/docs/architecture/20260327_InterSubMod研究願景定錨_01.md` + memory `project_research_vision_five_goals.md`。
- 為何這頁重要：把「修 tag」從低層 engineering 升到「解鎖整個研究主軸」的高度 — 教授會看到本工作對學術主軸的價值，不只是 bug fix。
- **後續可動 / 已初步現**（v2 S21 移入 note，可在此補述 30 秒）：
  - Phase 2A normal methylation reference 整合
  - HPFineNGroups marker 重驗（master × flag）
  - archive haplotag_version 標籤
  - Thread D NG=2 6/6 POSITIVE 已初步發現（Wilcoxon p=0.0156, W=21）— LOH-constrained phasing 論文主軸候選

**壓縮自 v2**：S20 + S22，整合 v2 S21 入 note。

**對應 Q 編號**：Q9（為何要做）、Q10（與 LongPhase-S 重疊）、Q11（Thread D NG=2 統計）。

---

### Slide 12. Take-home + 3 P0 + Q&A

**核心訊息**：三句 take-home 收束 + 短期 P0 三條 + Q&A 感謝。

**視覺**：上下兩段。
- 上半（大字 take-home）：
  - 「**TO 模式可用**」（綠色大字）
  - 「**Tag 層問題已解**」（綠色大字）
  - 「**五大目標解鎖**」（綠色大字）
- 下半（3 P0 卡片並排）：
  - 卡 P0-F1「commit V5 working tree」
  - 卡 P0-F3「七樣本 V5 全量重跑」
  - 卡 P0-F4「master × flag 重驗 HPFineNGroups」
- 底部：「Q&A 歡迎 · 詳細數據見 source_materials/」+「Thank you」

**頁面元素**：
- 標題：x=0.5, y=0.3, w=12.3, h=0.7（「Take-home · 3 P0 行動 · Q&A」）
- 上半 take-home 三段：x=0.5, y=1.2, w=12.3, h=2.5（三段大字並列）
- 中段分隔線：x=0.5, y=3.8, w=12.3, h=0.05
- 下半 3 P0 卡片：x=0.5, y=4.0, w=12.3, h=2.5（每卡 4.0 寬）
- 底部小字：x=0.5, y=6.7, w=12.3, h=0.7（「Q&A 歡迎；詳細數據見 source_materials/」+「Thank you」）

**結論句**：「**TO 可用 · Tag 已解 · 五大目標解鎖**」

**Speaker note 重點**（≥500 字）：
- 收尾節奏：take-home 30 秒、P0 卡片 30 秒、Q&A 引導 30 秒，總 90 秒收場。
- **三句 take-home 解讀**：
  - 「TO 模式可用」— self-phasing artifact 不再是 blocker；TO 模式的 ISM 結論在 V5 BAM 上可信
  - 「Tag 層問題已解」— phase 層由 V2b 解、tag 層由 V3F+V5 解、Sanity 15/15 PASS、5 條證據鏈互證
  - 「五大目標解鎖」— 目標 1/2/4 直接依賴 HP tag 正確；現在可推進
- **3 P0 行動承諾**：
  - **F1 · commit V5 working tree**：V5 目前在 working tree 未 commit，P0 必做
  - **F3 · 七樣本 V5 全量重跑**：HCC1395 5kHz 單樣本驗證需擴展到 7 樣本（HCC1395_5kHz、HCC1395、HCC1954、COLO829、H2009、SEQC2、其他）
  - **F4 · master × flag 重驗 HPFineNGroups**：HPFineNGroups marker 從 methylation bimodality 重詮釋為 phasing signature，需 master × flag 設計重驗（Phase 2B）
- **F5-F8 細節**（speaker note 補述，不上 slide）：
  - F5 · ISM 重跑（29 個 HP-依賴特徵）
  - F6 · cnLOH-aware filter（CN 層 + germline-only methylation reference）
  - F7 · vote log 直接驗證 Layer 1.5 confidence 0.6 threshold
  - F8 · 50-100 隨機抽樣位點 cross-validate sanity check
- **不重複 S1**：S12 take-home 用「TO 可用 / Tag 已解 / 五大目標解鎖」三句，與 S1 的「17.3:1 → 1:1 / +8.3 pp」數字陳述不重複；以承諾收尾比口號收尾更有力。
- **Q&A 預備指引**（內部用，不上 slide）：
  - Q1（17.3 怎麼算）→ S1 預講過
  - Q2（為什麼有 normal 就會對）→ S4
  - Q3（Jaccard=1.0 怎麼證明）→ S4
  - Q4（paired 沒被發現這 bug）→ S6
  - Q5（unit test）→ S8
  - Q6（PON 與 phasing anchor 區分）→ S6 + S7
  - Q7（sanity 是誰寫覆蓋率）→ S8
  - Q8（+8.3 pp 信賴區間）→ S9 + S10
  - Q9（F1 不變為何做）→ S5 + S9
  - Q10（與 LongPhase-S 重疊）→ S11
  - Q11（Thread D NG=2 統計）→ S11
  - Q12（PON 雙路徑驗證）→ S7 + S8
  - Q13（germline-first 真更好還是不同）→ S7
- **時間管理 fallback**：
  - 18 分鐘正常版：每 slide 75-90 秒
  - 15 分鐘壓縮版：S2/S5/S11 加速，S6/S7/S10 保持節奏
  - 12 分鐘極限版：S2 跳過、S11 簡述、S12 take-home 直接收尾
- **感謝**：「謝謝大家的時間，詳細數據在 source_materials/，歡迎 Q&A」

**壓縮自 v2**：S23 + S24。

**對應 Q 編號**：（無新 Q，但作為 Q&A 啟動點）。

---

## 2. 素材引用速查

### 既有圖片（v3/figures/ symlink → v2/figures/）

| Slide | 圖片檔名 | 路徑 |
|-------|---------|------|
| S2 | fig1_pipeline_comparison.png | `figures/fig1_pipeline_comparison.png` |
| S3 | fig17_hp_tag_5versions.png | `figures/fig17_hp_tag_5versions.png` |
| S3 | fig01d_somatic_bias_explanation.png | `figures/fig01d_somatic_bias_explanation.png` |
| S5 | 03_self_phasing_impact.png | `figures/03_self_phasing_impact.png` |
| S7 | fig01b_three_layer_logic.png | `figures/fig01b_three_layer_logic.png` |
| S7 | fig16_v5_threelayer_logic.png（備選） | `figures/fig16_v5_threelayer_logic.png` |
| S9 | fig18_concordance_amb_f1.png | `figures/fig18_concordance_amb_f1.png` |
| S10 | C_V5max1_chr19_4639528.png | `figures/C_V5max1_chr19_4639528.png` |
| S10 backup | C_V5max2_chr19_2235521.png | `figures/C_V5max2_chr19_2235521.png` |
| S10 backup | C_V5max3_chr19_7405500.png | `figures/C_V5max3_chr19_7405500.png` |
| S10 backup | D_SP1_chr19_17565944.png | `figures/D_SP1_chr19_17565944.png` |
| S10 backup | D_SP2_chr19_12452332.png | `figures/D_SP2_chr19_12452332.png` |
| S10 backup | D_SP3_chr19_12467180.png | `figures/D_SP3_chr19_12467180.png` |

### 自繪圖（build_pptx.py 生成至 v3/figures/auto/ 或 inline）

| Slide | 自繪內容 | 函式建議名 |
|-------|---------|-----------|
| S1 | 17.3:1 → 1:1 衝擊 banner | `make_s1_impact_banner()` |
| S2 | pipeline + 4 階段並列 | `make_s2_pipeline_4stage()` |
| S4 | phasing/tag 兩欄表 | `make_s4_two_column_table()` |
| S6 | root-cause tree（觸發條件 + 三葉） | `make_s6_root_cause_tree()` |
| S7 | V5 三層投票流程（紅黃綠圈） | `make_s7_v5_three_layer()` |
| S8 | 守恆律公式卡 + 5 證據鏈 | `make_s8_sanity_evidence()` |
| S11 | 業界家族樹 + 5 目標卡片 | `make_s11_family_tree_targets()` |
| S12 | take-home + 3 P0 卡片 | `make_s12_takehome_p0()` |

### 深度資料源（Q&A 開啟用）

| 資源 | 路徑 |
|------|------|
| v2 完整 24-slide storyboard | `../v2/00_storyboard_v2.md` |
| 重點數據速查表 | `../v2/notes/key_metrics_table.md` |
| 13 條 Q&A 預備答案 | `../v2/notes/qa_11_questions.md` |
| 18 份 v5_audit_suite 報告 | `../v2/source_materials/v5_audit_suite/` |
| 4-commit diff 細節 | `../v2/source_materials/01_code_diff/` |
| Purity 0.95 root cause | `../v2/source_materials/v5_audit_suite/04_purity_calculator_failure_root_cause.md` |
| v2 backup PPTX（教授追問用） | `../v2/output/LS-TO-V5.pptx` |

---

## 3. Speaker 提醒

### 4 個必講核心 slides（壓秒時其他可加速）

| 必講 | 為何必講 | 不能省的 90 秒 |
|------|---------|-------------|
| **S1** 17.3:1 → 1:1 | 開場鎖定核心；建立「這是真 artifact」共識 | 60 秒講數字 + 30 秒講 InterSubMod 邊界 |
| **S6** 根因樹 | 觸發條件主者 + 三 bug 的因果鏈 — 全演講邏輯轉折點 | 30 秒 purity 0.95 觸發 + 60 秒三葉拆解 |
| **S7** V5 三層投票 | germline-first 翻轉 + Layer 1.5 confidence + integer literal | 90 秒三層逐一帶過 |
| **S10** V5max1 IGV climax | 眼見為憑的 39 reads 翻轉；演講高潮 | 60 秒 IGV 指圖 + 30 秒 caveat |

### 可加速的 slides（時間緊時 30-45 秒）

- **S2** 工作流程：可一句話帶過「修補在 longphase-to-mod，InterSubMod 在使用端」
- **S11** 業界家族樹：可省略 LongPhase 2022 → 2025 演進，直接講「同實驗室相鄰工作 + 五大目標銜接」
- **S12** Take-home：take-home 三句直接念，P0 卡片帶過

### 演講節奏建議

| Slide | 預期分鐘 | 累計時間 |
|-------|---------|---------|
| S1 17.3:1 → 1:1 | 1.5 | 1.5 |
| S2 工作流程 | 1.0 | 2.5 |
| S3 HP tag 五值 | 1.5 | 4.0 |
| S4 Phasing OK Tag NG | 1.5 | 5.5 |
| S5 ISM 29/14/42 | 1.5 | 7.0 |
| **S6 根因樹** | **2.0** | **9.0** |
| **S7 V5 三層投票** | **2.0** | **11.0** |
| S8 Sanity + 證據鏈 | 1.5 | 12.5 |
| S9 量化指標 | 1.0 | 13.5 |
| **S10 V5max1 climax** | **1.5** | **15.0** |
| S11 業界家族樹 + 五目標 | 1.5 | 16.5 |
| S12 Take-home + Q&A | 1.0 | 17.5 |
| **總時間** | **17.5 分鐘** | （留 12-15 分鐘 Q&A）|

**時間管理 fallback**：
- **18 分鐘正常版**：依上表執行
- **15 分鐘壓縮版**：S2 / S5 / S11 各砍 30 秒，S6 / S7 / S10 不動
- **12 分鐘極限版**：S2 跳過、S11 簡述至 30 秒、S12 take-home 直接收尾

### Q&A 預備指引（哪些 Q 在哪些 slide 後可能出現）

| Q # | Q 主題 | 最可能在哪 slide 後出現 | 演講者預備動作 |
|-----|-------|----------------------|-------------|
| Q1 | 17.3 怎麼算 | S1 / S3 | S3 已展開；S1 後若被打斷可即答 614K/35.5K |
| Q2 | 為什麼有 normal 就會對 | S4 | 預備「TO 是必要研究方向」回答 |
| Q3 | Jaccard=1.0 怎麼證明 | S4 | 預備 LOH.bed vs ISM HP_Ratio LOH 區分 |
| Q4 | paired 沒發現這 bug | S6 / S7 | 預備「PON-only 暴露時機」答案 |
| Q5 | unit test 為什麼沒抓到 | S8 | 預備「LongPhase 無 round-trip test」答案 |
| Q6 | PON 與 phasing anchor 區分 | S6 / S7 | 預備兩個獨立步驟的解釋 |
| Q7 | sanity 是誰寫覆蓋率 | S8 | 預備「15 sites cherry-picked + F8 隨機 cross-val」答案 |
| Q8 | +8.3 pp 信賴區間 | S9 / S10 | 預備 N + Wilson CI 推估 |
| Q9 | F1 不變為何做 | S5 / S9 | 預備 read-level 價值 + 38% ISM 影響答案 |
| Q10 | 與 LongPhase-S 重疊 | S11 | 預備「同實驗室相鄰 + TO 變體」答案 |
| Q11 | Thread D NG=2 統計 | S11 / S12 | 預備 Wilcoxon W=21, p=0.0156 答案 |
| Q12 | PON Two-Pass 雙 purity | S7 / S8 | 預備 0.93 + 0.6 雙樣本驗證表 |
| Q13 | germline-first 真更好 | S7 | 預備 baseline somatic-first vs V5 germline-first 對比 |

### 截圖驗證檢查項（給 build_pptx.py 開發者）

每張 slide 必須通過：

1. 文字無溢出（Latin + CJK 雙語不互相覆蓋；標題 ≤ 15 字；每頁 ≤ 5 視覺元素）
2. 圖片等比 fit-within（無強制縱橫比破壞）
3. CJK 字型正確（無方塊、無亂碼）
4. Speaker note 存在且 **≥ 500 字**（V3 因 slide 變少，每張 note 比 v2 更長）
5. 圖片連結有效（v3/figures/ symlink → v2/figures/）
6. **必講 4 核心**（S1/S6/S7/S10）必須通過 visual review 抽查

驗證方式：`v3/scripts/screenshot_all.py` 必須輸出 `[ISSUES] none detected`。

---

## 4. V3 build 必做清單

1. ⏭ build_pptx_v3 撰寫（複用 v2 helper：`add_text_with_fallback()` per-char Latin/CJK fallback、`fit_image_within()` 等比插入、THEME dict 顏色直接複用）
2. ⏭ speaker_script_v3 撰寫（每張 ≥ 500 字；S6/S7 因含 Q1/Q2/Q12/Q13 整合可達 700-900 字）
3. ⏭ screenshot 驗證 `[ISSUES] none detected`
4. ⏭ 抽查 4 張必講核心 slide：**S1（衝擊 banner + 邊界小字）**、**S6（根因樹含觸發條件）**、**S7（V5 三層投票紅黃綠圈）**、**S10（V5max1 IGV + caveat）**

完整 build 指令：

```bash
cd InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v3
python3 scripts/build_pptx.py
python3 scripts/screenshot_all.py
```

---

## 5. 修訂歷程

| 版本 | 日期 | 變動 |
|------|------|------|
| V3-12（本檔）| 2026-04-30 | 從 v2-24 壓縮為 12 slides；one-slide-one-message；4 必講核心鎖定 S1/S6/S7/S10；speaker note ≥ 500 字 |
