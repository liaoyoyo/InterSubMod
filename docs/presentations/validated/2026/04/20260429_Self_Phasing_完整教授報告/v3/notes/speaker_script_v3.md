# Self-Phasing 教授報告 V3 精簡版 — 演講稿

> 產出日期：2026-04-30｜對應 PPTX：v3 12-slide 精簡版｜目標時長：15-18 分鐘 + Q&A
> 來源整合：v2/notes/speaker_script_v2.md（26 slides 完整版）、v2/notes/key_metrics_table.md、v2/notes/qa_11_questions.md（Q1-Q13）
> 嚴守 5 處口徑校準：(1) V5 不修 self-phasing 本身（V2b 已處理）；(2) 「同實驗室相鄰工作」非「業界共識」；(3) ISM 影響只列 count（29/14/42 共 85），來源報告誤寫 7%；(4) `Util.h:20-25`（不是 21-25）；(5) +8.3 pp（clean PS blocks，全基因組）。

---

## 演講前準備

| 項目 | 說明 |
|------|------|
| **受眾** | 教授（NGS / long-read 基本背景，不熟 LongPhase-TO 內部、不熟 InterSubMod 細節） |
| **目標時長** | 15-18 分鐘正講 + Q&A（保留 12-15 分鐘） |
| **與 v2 關係** | V3 是 v2 的精簡演講版；深度資料完整保留在 v2/source_materials 與 v2/notes |
| **4 個必講核心** | **S1**（開場衝擊）、**S6**（根因樹）、**S7**（V5 三層投票）、**S10**（V5max1 climax） |
| **預備 backup** | v2 PPTX 在 `v2/output/`；遇深度提問可隨時切換到 v2 對應 slide |
| **5 張最低核心**（時間嚴重不足時保留） | S1、S6、S7、S10、S12 |
| **緊急加速** | S2 / S5 可壓 30 秒；S8 / S9 可只報核心數字；不可壓：S1 / S6 / S7 / S10 |

---

## 開場分層分工說明（30 秒，講在 S1 之前）

「各位老師好。今天用 15 到 18 分鐘把 self-phasing 這個故事完整講一遍。先說明一下今天這份簡報的定位：這是 V3 精簡版，從原本 v2 的 26 張壓到 12 張，每張都是 one-slide-one-message。**所有深度數據、完整證據鏈與 13 條 Q&A 預備答案都保留在 v2 的 source_materials 跟 notes，需要查證隨時可以切回去**。今天演講分四段：開場衝擊一張、問題與根因六張、修補與驗證四張、收束承諾一張。Q&A 隨時打斷，謝謝。」

---

## Slide 1. 17.3 : 1 → 1 : 1（1.0 min · 必講）

### 開場（10 秒）
「先看這個對照數字 — 17 比 1 對 1 比 1。這就是今天整場演講要解的問題。」

### 主敘事（45 秒，必講核心，≥ 250 字）

「一句話定位：原始 LongPhase-TO 在 ClairS-TO Tumor-Only 模式下出現 self-phasing。我們在 HCC1395 5kHz 全基因組 baseline 觀察到，HP1 family 對 HP2 family 的 read 比是 614,000 對 35,500，等於 17.3 比 1，94.6% 的 reads 集中在 HP1，而且這個現象**跨 23 條染色體都一致**。生物學上這比應該是 1 比 1，因為 germline het 在染色體上是隨機分散的，所以 17.3 比 1 必然是 phasing 階段引入的 artifact。我們在 longphase-to-mod fork 用 4 個 commit 漸進修補（V2b → V3F → INDEL guard → V5），把它降回大約 1 比 1，**read-level concordance 對 paired ground truth 提升 +8.3 pp，這是 clean PS blocks 在全基因組的數字**。

兩個一定要先建立的認知：第一，**InterSubMod 本 repo 沒有任何 C++ 改動**，修補全部在外部 fork 的 longphase-to-mod，ISM 是被動受惠者。第二，self-phasing 問題鏈是被分層處理的 — phase 層的 self-phasing scaffold 由 V2b PON-only phasing 解掉；V3F 跟 V5 解的是 tag 層的 priority bug 跟 enum mismatch；**V5 不修 self-phasing 本身，這個分層後面 S10 climax 我會再強調一次**。

本場演講分四段：問題（S2-S5）、根因與修正（S6-S7）、驗證（S8-S10）、影響與收束（S11-S12）。」

### 過場（5 秒）
「先進工作流程，定位修補的位置。」

### 預期 Q&A
- Q：「17.3 怎麼算？」 → A：HP1 family（HP1 + HP1_1）/ HP2 family（HP2 + HP2_1）= 614K / 35.5K = 17.3。詳見 qa_11_questions.md Q1。
- Q：「Paired 模式是多少？」 → A：~1:1，HP_Ratio 落在 0.49-0.51，全 7 樣本一致。後面 S4 會展開。
- Q：「跨樣本一致嗎？」 → A：CV-2 同方向 7/7 樣本，Cohen's d = -1.20。

---

## Slide 2. 工作流程一覽（含 InterSubMod 邊界）（1.5 min）

### 開場（10 秒）
「先看 pipeline。哪個 box 是 caller、哪個是 phasing、哪個是我們本 repo，三件事必須分清楚。」

### 主敘事（70 秒，≥ 250 字）

「Pipeline 從左到右是這樣的：tumor BAM 進來，先送 ClairS-TO 這個外部 caller 產生 VCF；VCF 跟 BAM 一起進 longphase-to-mod 做 phasing 跟 haplotagging，輸出兩個東西：phased VCF 跟一個帶 HP:i tag 的 tumor_tagged.bam；最後 InterSubMod 讀 tagged BAM，計算 region-level 的 ISM 特徵。**今天 V5 修補的位置是中間的 longphase-to-mod，不是 ClairS-TO，也不是 InterSubMod**。再強調一次：HaplotagProcess.cpp 在 longphase-to-mod 這個獨立的 git repo，本 InterSubMod repo 完全沒有 C++ 改動。

為什麼不直接在 ISM 修？這違反單一職責 — phasing 跟 haplotag 該在 longphase-to-mod 處理。而且 ISM 的 TO F1 增益本來只有 0.0124，對 caller 端的 germline FP 沒有修復力，加再多甲基化特徵都沒辦法取代上游 tag 品質的修補。

四階段工作的 framing：(1) 機制定位 — 從理論預期 1:1 對照觀察到的 17.3:1，推導 self-phasing 機制；(2) 4-commit 修補 — V2b → V3F → INDEL guard → V5，每個 commit 解一層 bug；(3) 驗證 — 4 項硬性 sanity check 全 PASS 加 5 條獨立證據鏈；(4) 影響評估 — 對 ISM 85 個特徵分 3-tier 並接到五大研究目標。」

### 過場（10 秒）
「下一張定義 HP tag 跟三層證據鏈，這是後面所有討論的基礎詞彙。」

### 預期 Q&A
- Q：「為什麼不在 ISM 修？」 → A：違反單一職責；ISM TO F1 增益 0.0124 本來就無修復力。詳見 v2 S3 speaker note。
- Q：「longphase-to-mod 跟原版 LongPhase 差多少？」 → A：4 commits, +68/-36 行，集中在 3 函數，介面契約一字未變（drop-in replace）。

---

## Slide 3. HP tag 五值 + 三層證據鏈（2.0 min · 合併原 v2 S3+S4）

### 開場（10 秒）
「這張一張說兩件事：HP tag 五值定義，跟為什麼我們敢說 17.3:1 是 artifact 的三層證據。」

### 主敘事（90 秒，≥ 250 字）

「先講 HP tag 編碼。BAM 裡每條 read 帶一個整數 HP tag，longphase-to-mod 用 6 種值：HP:i:0 是 unphased；HP:i:1 跟 HP:i:2 是 germline 的兩個 hap；HP:i:11 跟 HP:i:21 是 somatic on HP1 跟 HP2；HP:i:33 是 somatic ambiguous。**這個編碼是 longphase-to-mod 自訂，不在 LongPhase 官方 spec**，跟 C++ enum（在 Util.h 第 20 到 25 行）的 HAPLOTYPE1_1 = 3 是不同數字 — 這個型別語意失配是後面 S6 enum bug 的伏筆，先記住「3 跟 11 是兩個不同數字但長得很像」這件事。

接下來是 LOH 兩層必須區分。**ISM 的 HP_Ratio LOH 是從 BAM HP tag 統計的（read-level），會受 self-phasing 影響，62% 是 artifact**；**LOH.bed 是 longphase-to-mod 直接從 VCF allele depth (AD) 偵測的（region-level），不受 self-phasing 影響，V5 前後 Jaccard = 1.0 完全不變**。兩套用不同數據源，過去把 ISM HP_Ratio LOH 當作 LOH.bed 是常見誤解。

最後是三層證據鏈，內外一致：(1) 理論層 — germline het 隨機分散，預期 1:1；(2) 全基因組層 — 觀察到 17.3:1 跨 23 染色體一致；(3) 個別位點層 — 後面 S4 會看到 SP1 chr19:17565944 baseline 113:0 全失衡。三條獨立證據都指向同一個結論：17.3:1 是 artifact。」

### 過場（10 秒）
「既然是 artifact，下一張就拆鎖：問題在 phasing 層還是 tag 層？」

### 預期 Q&A
- Q：「HP:i:33 是怎麼來的？」 → A：兩 hap 都不確定時的 fallback bucket，baseline enum bug 讓它幾乎不出現（0/15 sites），V5 修補後正確顯示。
- Q：「LOH.bed 跟 ISM HP_Ratio LOH 為何用不同數據源？」 → A：LOH.bed 走 VCF AD（區段層），ISM HP_Ratio 走 read-level HP tag，本身就是兩個維度。kappa = 0.670 不一致由此解釋。

---

## Slide 4. Phasing OK，Tag 有問題（1.5 min · 合併原 v2 S5+S6+S7）

### 開場（10 秒）
「拆鎖關鍵：self-phasing 影響哪一層？答案是 tag 層，不是 phasing 層。」

### 主敘事（70 秒，≥ 250 字）

「左欄 phasing 層：路徑是 VCF allele depth → LongPhase region 偵測，輸出 BED 格式 region-level LOH。**baseline LOH.bed 跟 V2b PON-only 啟用後的 LOH.bed 比對，Jaccard = 1.0 完全相同**，phasing 層的 LOH region 偵測不受 self-phasing 影響。V2b 的 PON-only phasing 已經把 phase scaffold 的 self-phasing 解掉了。

右欄 tag 層：路徑是 BAM HP:i tag → ISM per-region HP_Ratio 統計。Baseline 17.3:1 嚴重偏差，**62% 的 ISM HP_Ratio LOH 在 paired 比對下其實是 artifact**，V5 修補對象就是這層。

兩個關鍵對照證據佐證 tag 層才是問題：第一，**同位點跨模式 paired vs TO 的 HP_Ratio 相關係數 r = 0.001，n = 288K pairs，幾乎完全無相關** — 表示 TO 的 HP_Ratio 完全不反映真實 haplotype。第二，**TO-only 標記為 LOH 的位點，在 paired 模式下 86.5% 是平衡的**，TO 的 LOH 標記大部分是 self-phasing artifact。同樣本同位點，paired HP_Ratio 落在 0.49-0.51，TO 落在 0.91-0.95，Cohen's d = -1.20，巨大效應量。這個對照排除了「TO 樣本本身就 LOH 多」的解釋。」

### 過場（10 秒）
「Tag 層問題確認後，下一張說明這個問題對 ISM 研究造成多大影響。」

### 預期 Q&A
- Q：「Jaccard = 1.0 是什麼比對？」 → A：baseline LOH.bed 跟 V2b PON-only LOH.bed 的 region-level Jaccard 相似度。詳見 qa_11_questions.md Q3。
- Q：「為什麼不直接用 paired 就好？」 → A：normal BAM 不一定可得（臨床、archive、cell-free DNA），TO 模式是必要研究方向。

---

## Slide 5. 29 / 14 / 42 個 features 受影響（共 85）（1.5 min）

### 開場（10 秒）
「為什麼 tag 層 bug 重要？因為 ISM 的 38% 特徵直接依賴 HP tag。」

### 主敘事（70 秒，≥ 250 字）

「先說 ISM 的核心架構：ISM 把 reads 分成 4 個 bucket — HP1（germline HP1）、HP1-1（somatic on HP1）、HP2、HP2-1，計算 region-level 特徵。如果 HP tag 偏成 17.3:1，這個 4-bucket 分群就是錯的，所有依賴它的特徵都不可信。

V5 修補對下游 ISM 特徵的影響量化為 3 個 tier：

**嚴重影響 29 個（HP-依賴）**：HP_Ratio（62% 假 LOH）、Potential_LOH、HPMergedDelta、HPMergedSig（方向反轉）、HPFineNGroups（含 NG=2 LOH-constrained phasing marker）。這 29 個必須在 V5 BAM 上重跑才能信任，舊結論需加註版本（haplotag_version = baseline 或 V5）。

**中度影響 14 個（間接污染）**：QualityScore（AUC 0.497，已移除 LOH penalty）、GlobalP（取 HP/Allele 最小值）、CramersV（取 HP/Allele 最大值）、VerificationClass（label_sig 含 HP 成分）。重評多數影響微弱。

**不受影響 42 個（無 HP 依賴）**：PairwiseMeanDist、AlleleDelta、AlleleP、Caller 特徵（AF / GQ / DP / SB）、甲基化矩陣本身（BAM MM/ML tag 不依賴 HP）、CpG 座標、region_methyl_mean。這 42 個結論穩固不需重測。

口徑校準：**slide 上只列 count 不列百分比** — 來源報告寫 14/85 = 7%，但 14/85 實際是 16.5%，來源百分比有誤。為避免傳遞錯誤數字，本頁全部用 count 表達。重跑工作可聚焦在嚴重 + 中度的 43 個，不需全 85 重跑。」

### 過場（10 秒）
「下一張 climax：這個 17.3:1 的根因到底是什麼？」

### 預期 Q&A
- Q：「HPFineNGroups 重詮釋是什麼意思？」 → A：過去評為 subclone marker（甲基化 bimodality），V5 後重詮釋為 phasing signature（LOH-constrained phasing），論文主軸候選。詳見 v2 S22 + Q11。
- Q：「29 個都要重跑嗎？」 → A：是，加上 14 個中度，共 43 個進 V5 重跑名單。

---

## Slide 6. ★ 根因樹（觸發條件 + 三 bug）（2.0 min · 必講核心）

### 開場（10 秒）
「核心新發現：self-phasing 不是單一 bug，而是一個觸發條件加上三層獨立故障。」

### 主敘事（90 秒，必講核心，≥ 250 字）

「根（藍色）是觸發條件：**Purity 0.927 ≤ 0.95**。LongPhase-to-mod 在 PhasingProcess.cpp 第 197 行硬編碼 `if (purity > 0.95)` 才會觸發 Two-Pass 路徑。HCC1395 5kHz 真實 purity > 0.95（純 tumor），但 baseline calculator 估計成 0.927（四捨五入 0.93）— 這個估計**正確、不是誤判**，但因為 ≤ 0.95 沒跨閾值，所以走「三條路」標準流程，這個流程**假設樣本含 normal contamination 但實為純樣本**，於是 somatic 反客為主進 phasing graph，暴露 tag 層 somatic-first 投票 bias，最後產生 17.3:1。這是「特定條件觸發」而非單純設計缺陷，重點是純樣本剛好落在 0.927 卡 0.95 邊緣。

接下來是三葉，三個獨立 bug：

**Bug 1（紅色，getVote priority）**：`HaplotagProcess.cpp:506-541` 裡 variantKeys 順序是 `{HP1_1, HP2_1}` 排第一，任一 somatic count > 0 立刻 break，germline 完全跳過。為什麼 paired 模式沒事？paired germline votes 多但 priority bug 仍存在，paired 有 normal 校正最終結果；TO PON-only 後 germline votes 急減，bug 立刻顯形。V3F commit 41ff147 修。

**Bug 2（橙色，enum vs int literal mismatch）**：C++ enum 在 **Util.h 第 20 到 25 行**（注意是 20-25 不是 21-25），HAPLOTYPE1_1 定義成 3。但 BAM HP tag integer 是 11。caller 端寫 `if(hpResult != HAPLOTYPE1_1)` 用 enum=3 比較 hpResult=11/21/33，fallback 死分支永不執行，HP:i:33 永不出現（baseline 0/15 sites）。C++ 把 enum 隱式轉 int 不會 warn，靜態檢查不易抓。V3F 同 commit 改為 integer literal。

**Bug 3（灰色，scaffold，已由 V2b 解）**：`PhasingProcess.cpp:154-157` 的 `convertNonGermlineToSomatic()` 不被觸發。已由 V2b commit 8b8c1fd 啟用 PON-only phasing 解決，本 PPT 不展開，視覺只佔 20%。

焦點重點：**今天主要講的是 V3F + V5 修補的 Bug 1 跟 Bug 2，V2b 的 scaffold 已經處理掉了**。」

### 過場（10 秒）
「下一張看 V5 怎麼修補這個三層 bug — 答案是三層投票流程。」

### 預期 Q&A
- Q：「為什麼 paired 沒被發現這 bug？」 → A：paired germline votes 多掩蓋 priority bug，paired 有 normal 校正最終結果；PON-only 後 germline votes 急減立刻顯形。詳見 qa_11_questions.md Q4。
- Q：「PON 有資料還會 self-phasing？」 → A：要區分 PON classification（VCF 層）vs phasing anchor（read 連結層）。baseline 預設 `--pon-only-phasing=false` 把 somatic 拿去當 anchor。詳見 Q6。
- Q：「purity 估錯了嗎？」 → A：沒估錯，0.927 vs 真實 >0.95 在誤差範圍內，問題是邊界 case 卡 0.95 閾值。

---

## Slide 7. ★ V5 三層投票（修正流程）（2.0 min · 必講核心）

### 開場（10 秒）
「V5 對 getVote() 的具體改動：從 V3F 兩層升級為三層投票邏輯。」

### 主敘事（90 秒，必講核心，≥ 250 字）

「**Layer 1（germline first，紅圈）**：`if germlineHP1 > germlineHP2 → 1; elif germlineHP2 > germlineHP1 → 2; else 進入 Layer 1.5`。保持 germline first 核心原則，避免 somatic 自我決定，不破壞已驗證的 V3F 行為。這一層解 Bug 1（priority bug），把 somatic-first 翻轉成 germline-first。

**Layer 1.5（somatic fallback，黃圈，V5 新增）**：當 germline 平手時參考 somatic 投票。`if somaticTotal == 0 → hpResult = 0`；`elif somaticHP1 > 0.6 * somaticTotal → 對齊 HP1 方向`；`elif somaticHP2 > 0.6 * somaticTotal → 對齊 HP2 方向`；`else → hpResult = 33（真正 ambiguous）`。Confidence threshold 0.6 是經驗值（R4 caveat，未來 F2 用 vote log 直接驗證）。**這是 mid-low purity 的防禦層** — 在 0.6 純度 sample 上，Layer 1.5 會把弱 directional 訊號標 ambiguous，避免錯誤分配。

**Layer 2（encode hpResult，綠圈）**：germline=1 加 somatic_total>0 編 HP:i:11；germline=2 加 somatic_total>0 編 HP:i:21；germline=0 加 somatic_total>0 編 HP:i:33；germline only 編 HP:i:1 或 2；nothing 編 HP:i:0。**這層的 caller 端 enum 比較改為 integer literal 11/21/33**，解 Bug 2。

效果：AMB% 從 17.5% 降到 8.0%（−9.5 pp）；HP:i:33 reads 從 239,679 降到 110,197（−54%）；**介面契約零變動**（HaplotagProcess.h 第 66-68 行 method signature 一字未變，drop-in replace）。

副作用是預期的：部分 borderline somatic 從 33 變成 11/21，這正是 Layer 1.5 設計目的 — 把卡在 ambiguous 的 reads 解到具體 hap。」

### 過場（10 秒）
「下一張驗證這套三層投票通不通過硬性 sanity check。」

### 預期 Q&A
- Q：「Confidence threshold 0.6 怎麼決定的？」 → A：經驗值，未來 F2 用 vote log 直接驗證。對 mid-low purity 的保守防禦設計。
- Q：「介面契約是什麼意思？」 → A：HaplotagProcess.h 第 66-68 行三個函數簽章一字未變，下游 caller 不需改任何呼叫。drop-in replace。
- Q：「程式碼 diff 看不太懂可以給看嗎？」 → A：完整 diff 在 v2/notes/code_references.md，三處紅標分別對應 Bug 1（priority loop 拆解）、V5 新增（Layer 1.5）、Bug 2（enum → integer literal）。

---

## Slide 8. Sanity 15/15 PASS + 5 條獨立證據鏈（1.5 min）

### 開場（10 秒）
「驗證段第一張：硬性 sanity check 全通過 + 5 條獨立證據互相收斂。」

### 主敘事（70 秒，≥ 250 字）

「左欄 4 項硬性 sanity check（HCC1395 5kHz 15 sites，v5_audit_suite/06）：

(1) **守恆律 A · Δ-consistency**：ΔHP33 + (ΔHP11 + ΔHP21) = 0，tag 移轉量平衡，in = out，**15/15 PASS**。
(2) **守恆律 B · Germline 不變**：V3F 與 V5 的 HP1 / HP2 reads 逐 site 比對，**0 reads 漂移，15/15 PASS**。
(3) **Layer 1.5 期望 1 · 33→directional 精確守恆**：ΔHP11 等於 V3F=33→V5=11 的 read 數，ΔHP21 同理，**15/15 PASS**（V5max1 chr19:4639528 39 reads 精確守恆）。
(4) **Layer 1.5 期望 2 · 無 germline → HP33**：跨 15 sites pool transition table，germlineResult ≠ 0 時不出 HP33，**0 violation**。

右欄 5 條獨立證據鏈：

(1) **理論層** — germline het 隨機 → 預期 1:1（S3）。
(2) **全基因組層** — HP1:HP2 = 17.3:1 跨 23 chr 一致（S1）。
(3) **個別位點層** — SP1 chr19:17565944 baseline 113:0 全失衡（S10 climax 旁證）。
(4) **Sanity check** — 4 項 15/15 PASS（左欄）。
(5) **程式碼最小必要** — V5 累計 +68/-36 行集中於 3 函數，介面契約零變動。

5 條獨立路徑同步收斂於同一結論：**self-phasing 修補必要且充分**。

Caveat：15 sites 是 cherry-picked，需要 7 樣本擴展（F3）+ 100 隨機位點 cross-validate（F8），整體 sanity check 是必要不充分。」

### 過場（10 秒）
「下一張看修補後的量化指標 — AMB / HP33 / concordance 三項硬數字。」

### 預期 Q&A
- Q：「Sanity check 是誰寫的？」 → A：本工作新增（v5_audit_suite/06 Agent D 產出），未來建議加入 longphase-to-mod CI test suite。詳見 qa_11_questions.md Q5+Q7。
- Q：「為什麼只 15 sites？」 → A：cherry-picked（5 TP + 4 FP + 3 V5max + 3 SP-extreme），需 F3 7 樣本擴展 + F8 100 隨機位點補完。

---

## Slide 9. 量化指標：AMB ↓ / HP33 ↓ / Concordance ↑（1.0 min）

### 開場（10 秒）
「量化指標三聯橫幅，三個都同方向改善。」

### 主敘事（45 秒，≥ 200 字）

「(1) **AMB% 從 17.5% 降到 8.0%**（−9.5 pp）— V5 解了過半 ambiguous reads，回收成 HP1/HP2 directional reads。

(2) **HP:i:33 reads 數從 239,679 降到 110,197**（−54%）— 在 HCC1395 全基因組層級，超過一半的 ambiguous reads 被正確重分配。

(3) **全基因組 clean PS blocks paired GT concordance：V5 = 90.5% vs Baseline = 82.2%，+8.3 pp**。**口徑校準必須加 caveat「(clean PS blocks，全基因組)」**，不是全基因組 raw，不是 cherry-picked。Clean PS 的定義是 germline accuracy ≥ 70% 的 high-quality phase blocks。

對照 15-site cherry-picked clean PS：V5 88.2% vs BL 74.9%，+13.3 pp（小樣本更戲劇化）。Aggregate 全 15 sites pool：V5 78.85% vs BL 72.20%，+6.65 pp。

注意：**SEQC2 F1 不變是預期** — V5 不改 caller，改的是 BAM HP tag 編碼，與 caller 輸出 VCF 無關。F1 變動只可能來自下游 ISM SuggestFilter，但 ISM TO F1 增益本來只有 0.0124，V5 vs Baseline 最終 F1 差 −0.0003 在噪音範圍。**真實價值在 read-level concordance +8.3 pp，不是 F1**。」

### 過場（10 秒）
「下一張 climax — 39 reads 100% 重分配的視覺證據。」

### 預期 Q&A
- Q：「+8.3 pp 的 N 多少？信賴區間？」 → A：全基因組 clean PS N 在數萬到數十萬量級；推估 Wilson 95% CI ±1 pp 內。15 site clean PS N=11。詳見 Q8。
- Q：「跨樣本一致嗎？」 → A：concordance 跨樣本擴展 F3 待辦；但 self-phasing 方向 7/7 一致已驗證。

---

## Slide 10. 🎯 V5max1 climax — 39 reads 100% 重分配（1.5 min · 必講核心）

### 開場（10 秒）
「演講高潮。chr19:4639528 V5max1 是最戲劇化的視覺證據。」

### 主敘事（70 秒，必講核心，≥ 250 字）

「看這張 4-BAM 並列 IGV：在 V3-Fixed panel，**一大群紫色 reads 標為 HP:i:33（ambiguous），39 條**；V5 修補後這 **39 條全部正確重分配為 HP:i:11**（somatic on HP1）。守恆律驗證：V3F 39 reads HP33 → V5 39 reads HP11，in = out 完全平衡，sanity check 守恆律 1 跟 3 同時 PASS。

這是 Layer 1.5 somatic fallback 的具體效果 — 原本因為 V3F 過度保守卡在 ambiguous bucket 的 reads，V5 用 confidence ≥ 0.6 的 somatic 投票把它們解到具體 hap。從 baseline → V2b → V3F → V5 的 4-BAM 並列，視覺差異很明顯：baseline 紅藍分離但方向錯（self-phasing），V2b 翻轉到正確 paired 方向（解 phase scaffold），V3F 出現 39 條紫色 ambiguous（enum bug 修了但 V5max 過度保守），V5 39 條紫色全部回到正確 HP1。

**底部 caveat 防誤解（reviewer R2 採納，必須口說一句）**：「但 V5 不修 self-phasing 本身（V2b 已處理 phase scaffold）」。**口徑校準關鍵**：高潮圖容易讓教授誤以為 V5 全面解決 self-phasing，實際上 self-phasing 問題鏈是被分層處理的 — phase 層的 self-phasing scaffold 由 V2b 的 PON-only phasing 處理（改 phasing graph anchor 來源，從含 somatic 改為僅 PON-confirmed germline）；V3F 解 Bug 1+2 priority + enum mismatch；V5 解 V3F 過於保守，新增 Layer 1.5 directional reassignment 把卡在 33 的 reads 解到 11/21。

**所以 V5 的「39 reads → 100% reassigned」是 directional reassignment，不是 phasing graph 重建**。再強調一次：V5 不是「解所有」，V5 是「解 V3F 過於保守 + 整個分層的最後一塊」。

旁證 3 個 SP-extreme 位點：SP1 chr19:17565944 baseline 113:0 → V5/V2b 翻轉為 0:113，paired 確認 110:2 是 HP2 為真；SP2 SP3 同樣 3/3 一致翻轉至 paired 方向。」

### 過場（10 秒）
「下一張收束影響面：業界家族樹 + ISM 五大目標解鎖鏈。」

### 預期 Q&A
- Q：「V5 翻轉 SP1 113:0 → 0:113 是 V5 的功勞？」 → A：不是，SP1-3 的修補主要是 V2b（phase 層），V5 在這些位點 Δ ≈ 0；V5 主戰場是 V5max1 這種 V3F 過度保守的 33 → 11/21。
- Q：「怎麼確認 V5max1 的 39 reads 不是隨機翻轉？」 → A：守恆律 3 精確驗證 — ΔHP11 = n(V3F=33→V5=11) 逐 read 配對，不是統計上的「相近」。

---

## Slide 11. 業界家族樹 + 解鎖五大目標（1.5 min · 合併原 v2 S20+S22）

### 開場（10 秒）
「上半看本實作的業界定位，下半看解鎖鏈。」

### 主敘事（70 秒，≥ 250 字）

「**上半業界家族樹**：上游基底是 LongPhase（Lin 2022 *Bioinformatics*）處理 germline phasing。由此分出兩個分支：(1) **LongPhase-S（bioRxiv 2025）是同實驗室的 paired 版**，在 ClairS 上 SNV F1 +4.5%、indel +7.1%，提供「somatic 錨在 germline scaffold」的 anchoring 概念。(2) **longphase-to-mod V5 是本實作**，在 TO + PON 條件下的本地實作。WhatsHap 跟 HapCUT2 不支援 TO 是公開工具的 gap，本實作填補 TO + PON 場景，與 LongPhase-S 形成姐妹關係（paired vs TO）。

**口徑校準關鍵**：本實作不寫「業界共識」也不是「標準替代」，而是「**同實驗室相鄰工作**」。LongPhase-S 為 paired 場景參考，不重疊本工作 TO 場景。

**下半解鎖五大目標**：InterSubMod 五大研究目標 — 目標 1 per-CpG 多標籤關聯、目標 2 clone 結構、目標 3 二次打擊順序、目標 4 TO normal 補強、目標 5 整合 evidence panel 提升 F1。**目標 1/2/4 直接依賴 HP tag 正確**：per-CpG 關聯需 hap 分層、clone 結構需 4-bucket、TO normal 補強需 phasing 正確。目標 3 間接依賴目標 1/2；目標 5 部分依賴。

研究發展樹：self-phasing 修補完成（V5）→ 4-bucket 分群可信 → NGroups / HPSig / HP_Ratio 可信 → 五大目標解鎖。**Self-phasing 修補是 5 大目標解鎖前提**，特別是目標 1/2/4。本工作不是一次性修補，而是支撐整個 InterSubMod 的根基。

已有初步發現的方向（全部依賴 V5 BAM）：(1) Thread D LOH-constrained phasing — NG=2 cross-sample 6/6 POSITIVE，Wilcoxon p=0.0156，Inner same-hap 93-99%（obs18 6 樣本），論文主軸候選；(2) HPFineNGroups marker 機制重詮釋為 phasing signature；(3) Phase 2A normal methylation reference 解鎖。」

### 過場（10 秒）
「最後收尾：take-home + 3 個 P0 行動承諾。」

### 預期 Q&A
- Q：「跟 LongPhase-S 重疊多少？算 contribution 嗎？」 → A：LongPhase-S 是 paired，本實作是 tumor-only + PON；設計哲學一致但場景不重疊。詳見 qa_11_questions.md Q10。
- Q：「為什麼 LongPhase-S 不直接用在 TO？」 → A：LongPhase-S 假設有 normal BAM；當 normal 不可得，需 PON 替代，本實作即為 TO 變體。詳見 Q12（reviewer 補充）。

---

## Slide 12. Take-home + 3 P0 行動 + Q&A（1.0 min）

### 開場（10 秒）
「最後三句結語 + 三個短期行動承諾，下次 meeting 教授會看到具體進度。」

### 主敘事（45 秒，≥ 200 字）

「**Take-home 三點**：

(1) **TO 模式可用** — 透過 V5 修補，TO + PON 條件下 tag 品質接近 paired 水準，read-level concordance +8.3 pp（clean PS blocks，全基因組）。

(2) **Tag 層問題已解** — V3F 解 priority bug + enum mismatch，V5 補 Layer 1.5 directional reassignment，AMB% 17.5%→8.0%，HP:i:33 reads −54%，介面契約零變動。

(3) **五大目標解鎖** — InterSubMod 目標 1/2/4 直接依賴 HP tag 正確，V5 修補是研究可信度前提。

**口徑校準分層 caveat**：V5 不修 self-phasing 本身（V2b 已處理 phase scaffold）。

**Next step 3 條 P0 承諾**：

- **F1 Commit V5 working tree** — 切 2 獨立 commits（Layer 1.5 + countSNPHaplotype alt guard），完成後 tag v5.0，加 manifest haplotag_version。
- **F3 7 樣本 V5 BAM 全量重跑** — 已驗證 HCC1395_5kHz，待 DORADO / HCC1937 / HCC1954 / H1437 / H2009 / COLO829，預估 ~10 hr parallel。
- **F4 master × flag 重驗 HPFineNGroups** — 判定 phasing signature 機制本質，對應 Thread D 主軸。

以承諾收尾比口號收尾更有力。Q&A 歡迎隨時打斷，深度數據在 v2/source_materials 三份報告與 v2/notes 6 份 onboarding 文件。Thank you。」

### 過場（5 秒）
「以上是 V3 精簡版重點。詳細數據請看 v2 完整版的 source_materials；現在開放問答。」

### 預期 Q&A
- Q：「V5 什麼時候能 commit？」 → A：F1 待辦，切 2 獨立 commits + tag v5.0 是下次 meeting 前完成。
- Q：「7 樣本擴展什麼時候有結果？」 → A：F3 預估 ~10 hr parallel，下次 meeting 應有 5-7 樣本初步結果。

---

## 收尾（10 秒，講在 S12 之後進 Q&A 之前）

「以上是 V3 精簡版重點，詳細數據請看 v2 完整版的 source_materials；現在開放問答。」

---

## Q&A 預備清單（13 條速查）

| ID | 問題（簡述） | 預期出現於 slide | 簡答關鍵字 |
|----|-------------|----------------|-----------|
| **Q1** | 為什麼是 17.3？這數字怎麼算？ | S1 / S5 | HP1 family / HP2 family = 614K / 35.5K = 17.3 |
| **Q2** | 為什麼 paired 沒事？（Paired ~1:1 根本原因） | S4 | normal 提供 germline het ground truth；TO 無 normal → somatic 反客為主 |
| **Q3** | Jaccard=1.0 怎麼證明 phasing 不變？什麼是 LOH.bed？ | S3 / S4 | LOH.bed = LongPhase 從 VCF AD 偵測；Jaccard 1.0 = baseline vs V2b BED 完全相同 |
| **Q4** | 為什麼 paired 沒被發現 priority bug？ | S6 | paired germline votes 多掩蓋 + paired 有 normal 校正；PON-only 後立刻顯形 |
| **Q5** | 這 bug 不是該有 unit test？（enum vs integer literal） | S6 | LongPhase 原 codebase 無 HP tag round-trip test；C++ 隱式轉型不 warn；新增 4 項 sanity check 入 CI |
| **Q6** | 為什麼 PON 有資料還會 self-phasing？ | S6 | 區分 PON classification（VCF 層）vs phasing anchor（read 連結層）；baseline `--pon-only-phasing=false` |
| **Q7** | Sanity check 是誰寫的？覆蓋率多少？ | S8 | 本工作新增（v5_audit_suite/06 Agent D）；15 sites cherry-picked；F3 + F8 待擴展 |
| **Q8** | +8.3 pp concordance 的 N 多少？信賴區間？跨樣本一致？ | S9 | 全基因組 clean PS N 數萬-數十萬；推估 Wilson 95% CI ±1 pp；F3 跨樣本待辦 |
| **Q9** | F1 不變是預期 — 那為什麼還要做？對發論文有意義嗎？ | S9 / S11 | 下游 ISM 38% 受影響；HPFineNGroups 重詮釋；五大目標解鎖前提；同 LongPhase-S 形成姐妹工作 |
| **Q10** | 跟 LongPhase-S 2025 paper 重疊多少？算 contribution 嗎？ | S11 | LongPhase-S paired vs 本實作 TO；設計哲學一致；TO gap + 4 sanity check + 7/7 一致為新貢獻 |
| **Q11** | Thread D NG=2 6/6 POSITIVE 是什麼統計？p=0.0156 對 6 樣本可信嗎？ | S11 | Wilcoxon signed-rank W=21 exact p=0.0156；paired control gap 0.00003 p=0.578；evidence grade B |
| **Q12** | PON 雙路徑（Two-Pass）在不同 purity 下都有效嗎？ | S6 / S7 | 0.93 + 0.6 sample 都驗證有效；mid-low purity 是 conservative tagging 防禦 |
| **Q13** | Baseline somatic-first vs V5 germline-first — V5 真的更好還是只是不同？ | S6 / S7 | V5 更好；mid-low purity 加 confidence 0.6 threshold 是設計上的防禦層 |

**Q14（reviewer 補充）**：cnLOH 區為何仍未解？ → 雙親同源無 het variants，V5 fallback 也無法 anchor，需 cnLOH-aware filter（F6 future work）。

詳細答案見 `v2/notes/qa_11_questions.md`。

---

## 計時表（Reference）

| Slide | 主訊息 | 預期時長 | 累計 | 備註 |
|------|--------|---------|------|------|
| S1 | 17.3 : 1 → 1 : 1 | 1.0 min | 1.0 | **必講核心** · 開場衝擊 |
| S2 | 工作流程一覽 | 1.5 min | 2.5 | Pipeline 邊界 |
| S3 | HP tag 五值 + 三層證據 | 2.0 min | 4.5 | 合併原 v2 S3+S4 |
| S4 | Phasing OK，Tag 有問題 | 1.5 min | 6.0 | 合併原 v2 S5+S6+S7 |
| S5 | 29/14/42 features 受影響 | 1.5 min | 7.5 | ISM 影響量 |
| S6 | ★ 根因樹 | 2.0 min | 9.5 | **必講核心** · Purity 0.927 觸發 + 三 bug |
| S7 | ★ V5 三層投票 | 2.0 min | 11.5 | **必講核心** · germline-first + Layer 1.5 + encode |
| S8 | Sanity 15/15 + 5 證據鏈 | 1.5 min | 13.0 | 4 守恆律 + 多源驗證 |
| S9 | 量化指標 | 1.0 min | 14.0 | AMB / HP33 / concordance |
| S10 | 🎯 V5max1 climax | 1.5 min | 15.5 | **必講核心** · 39 reads 翻轉 + caveat |
| S11 | 業界家族樹 + 五大目標 | 1.5 min | 17.0 | LongPhase 世系 + 解鎖鏈 |
| S12 | Take-home + 3 P0 + Q&A | 1.0 min | **18.0** | 收束 + 行動承諾 |

**演講總時長**：18.0 分鐘（含過場與 take-home）；Q&A 預期 12-15 分鐘；總時長預估 30-33 分鐘。

---

## 緊急加速指引（時間不夠時）

| 場景 | 動作 |
|------|------|
| **時間餘 17 min** | 全照計時表講 |
| **時間餘 15 min** | S2 / S5 各壓 30 秒（背景說明簡化） |
| **時間餘 13 min** | 再壓 S8 / S9 各 30 秒（只報核心數字） |
| **時間餘 10 min** | 上述全壓；S11 業界家族樹只說一句「同實驗室相鄰工作」 |
| **時間餘 < 10 min** | **僅保留 5 張最低核心**：S1 / S6 / S7 / S10 / S12（其餘跳過或速講） |

**不可壓的 4 個必講核心**：S1（開場衝擊）、S6（根因樹）、S7（V5 三層投票）、S10（V5max1 climax）。

**口徑校準絕不能漏的 5 點**（任何時候都要講到）：
1. V5 不修 self-phasing 本身（V2b 已處理）— S10 caveat 必說
2. 「同實驗室相鄰工作」非「業界共識」— S11 必說
3. ISM 影響只列 count（29/14/42 共 85），來源報告 7% 有誤 — S5 注意說法
4. `Util.h:20-25`（不是 21-25）— S6 Bug 2 必說
5. +8.3 pp（clean PS blocks，全基因組）— S1 + S9 必加 caveat

---

## 演講者準備檢核清單

| 項目 | 確認 |
|------|:----:|
| 5 處口徑校準熟記 | ☐ |
| 4 個必講核心節奏熟練（S1/S6/S7/S10）| ☐ |
| S6 根因樹「Purity 0.927 ≤ 0.95 觸發 + 三 bug」一句話可說清 | ☐ |
| S7 V5 三層投票口語化（germline first / fallback 0.6 / encode 11/21/33） | ☐ |
| S10 climax 一句說清「V5 不修 self-phasing 本身」 | ☐ |
| 13 條 Q&A 預備答案各 30 秒節奏 | ☐ |
| 計時演練：18 分鐘內可講完 12 slides | ☐ |
| Backup：v2 PPTX 隨時可切換 | ☐ |

---

## 文件結束

如需深度資料請查：
- `v2/notes/speaker_script_v2.md` — 26 slides 完整講稿（30 分鐘）
- `v2/notes/key_metrics_table.md` — 所有重點數據與來源行號
- `v2/notes/qa_11_questions.md` — Q1-Q13 完整問答
- `v2/source_materials/` — 三份原始報告（IGV visual audit / V5 audit / longphase TO vs V5 技術報告）
