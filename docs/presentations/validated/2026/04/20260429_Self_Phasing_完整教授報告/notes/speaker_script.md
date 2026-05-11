<!--
build_date: 2026-04-29
audience: 演講者本人（講給教授聽）
status: validated
purpose: 逐 slide 演講稿，作為 PPT 講解時的具體 script
ppt_path: ../output/20260429_Self_Phasing_完整教授報告.pptx
storyboard: ../00_storyboard.md
-->

# Self-Phasing 完整教授報告 · Speaker Script

> **使用方式**：每張 slide 配 1-3 段口語化演講稿。重點關鍵字以 **粗體**標示。
> 預估總時長 35-45 分鐘（每張 slide ≈ 1 分鐘，視 Q&A 調整）。

---

## Section 0 · 開場 (S1)

### S1 · Self-Phasing：原因、修補、驗證

教授你好。今天的報告主題是 **self-phasing** — 在 ClairS-TO 模式下我們最近發現的一個機制問題，
以及我們是怎麼修、怎麼驗證的。

報告分七個段落：先用 5 張 slide 補完所有需要的基本定義，避免術語混用造成的誤解；
接下來是問題敘述、量化證據、修補方向、改動細節、驗證結果，最後是結論與後續工作。

兩個重點預告：
1. 這次的修補**不是動 InterSubMod**，是動另一個工具叫 **longphase-to-mod**，是獨立的 git repo。
2. 修完之後 **SEQC2 F1 變動只有 -0.0003**，這個我們等下會詳細解釋為什麼是預期行為，
   真實價值在 read-level 的 tag 品質提升 +8.3 個百分點。

歡迎隨時打斷提問。

---

## Section 1 · 基本定義 (S2-S6)

### S2 · 三條 pipeline 上下游關係

先看上下游關係。從左到右三個工具：
- **ClairS-TO** 是 caller，輸入是 tumor BAM，輸出 VCF（裡面有 FILTER、AF、GQ 這些 caller-level 訊息）
- **LongPhase-TO** 是 phasing 和 haplotag 工具，輸出 phased VCF、tagged BAM、LOH.bed
- **InterSubMod** 是下游消費者，讀 BAM 的 HP tag 算 region-level 特徵

這三個工具的職責切分**很重要**，過去常常被混為一談，例如把 ISM 算的 HP_Ratio LOH 當成 LOH.bed —
這是後面 S14 會詳談的常見誤解。

### S3 · Phasing 是什麼？

Phasing 的目的是把不同位點的 variants 串成同一條 haplotype。
做法是看 long reads 上是否同時帶有兩個 ALT — 如果是，就把它們判定為同一 haplotype。

圖中藍色 reads 跨多個 SNV 都帶 A，被歸到 HP1；橘色帶 T 歸到 HP2。
每個 phase block 內部 GT 用 **`|`** 分隔（例如 `1|0`），
**PS tag** 是 phase-set ID，同一個 PS 表示這些 variants 在同一個 phase block。

這個概念是後面要解釋 self-phasing 的基礎。

### S4 · 三層資料的角色

接下來這張很關鍵：phasing pipeline 有**三層資料**，三者**不可混用**。
- Caller layer 看 VCF：FILTER、AF、GQ
- Phasing layer 看 phased VCF GT、PS、LOH.bed
- Haplotag layer 看 BAM 的 HP:i tag

例如 ISM 算的 HP_Ratio LOH 是 read-level 統計，和 LOH.bed 的 region-level 標註是**不同層次**的東西 —
過去我們在某些報告中把兩者混為一談，是 LOH 結論誤判的根因。

### S5 · HP tag 的 5 個整數值

HP tag 是 BAM 檔內每條 read 的整數標記。在 longphase-to-mod 裡總共有 **6 種**可能值：
- HP:i:1 = germline HP1
- HP:i:2 = germline HP2
- HP:i:11 = somatic on HP1（germline HP1 之上的 somatic）
- HP:i:21 = somatic on HP2
- HP:i:33 = somatic ambiguous（兩 hap 都不確定）
- HP:i:0 = unphased（無 phase 資訊）

這個編碼是 longphase-to-mod 自訂，**不在 LongPhase 官方 spec**。
S22 會講 V5 final 的對應表，這裡先建立認知。

### S6 · InterSubMod 的下游角色

ISM 強依賴 HP tag — distance metric、HPSig、HPFineNGroups、HP_Ratio 都需要它。
**所以上游任何 bias 都會直接傳遞下來。**

但 ISM 自己**不修補** phasing/haplotag 問題（這是 longphase-to-mod 的職責）。
另一個重要事實：ISM 對 TO germline FP 的修復力是 0.0124（vs paired 0.0909），
所以 ISM 加再多特徵也**取代不了**上游修補。

---

## Section 2 · 問題敘述 (S7-S10)

### S7 · Self-phasing 是什麼？

進入問題本身。**Self-phasing 的核心機制**：
同一個 sub-clone 的 somatic variants 由相同的 read population 攜帶，
所以 long-read 跨多個 somatic 位點都同時帶 ALT。
phasing 演算法看到這種強連結，把這些 somatic 全部塞到同一個 phase block。

結果：**本應隨機 50:50 分到兩個 hap，變成全部偏向同一邊**（例如 94.6% 集中在 HP1）。
這個現象在 paired 模式不會出現，TO 模式才會爆發。

### S8 · AF = 0.3 走例

這張是教授必看的核心圖，用 AF=0.3 的具體例子走一遍：

- **Paired 模式**：有 matched normal 對照，phasing 用 normal 的 germline 當 anchor，
  HP_Ratio 維持 ~0.5（隨機分到兩 hap）
- **TO 模式**：沒有 normal，phasing 只能用 tumor 自己的 reads，
  同 sub-clone 的 somatic reads 連在一起，**HP_Ratio 跑到 0.94**（94% 都歸到一個 hap）

這就是 self-phasing 的具體展現。

### S9 · TO 模式為什麼特別嚴重？

TO 模式 self-phasing 嚴重度是**三個條件疊加**：
1. **無 matched normal** — phasing 失去外部 anchor
2. **PON 不完美** — PON 標已知 SNP，但 rare het 仍漏入 phasing
3. **預設 `--pon-only-phasing=false`** — somatic 直接進 phasing graph 當 anchor

V2b commit 把這個預設改成 true，是修補的第一步。

### S10 · 三層 bug 同時暴露

但更複雜的是，這**不是單一 bug**，是**三層獨立故障**在 PON-only 啟用後集中暴露：
- (A) Phase scaffold 用 somatic anchor — `convertNonGermlineToSomatic()` 不被觸發
- (B) `getVote()` priority bug — somatic vote 搶先 germline
- (C) enum vs HP integer literal mismatch — `if(hpResult != HAPLOTYPE1_1)` 永不匹配，HP:i:33 永不出現

三層獨立、需要三個獨立修補 — 這是後面 4-commit 漸進演進的背景。

---

## Section 3 · 量化證據 (S11-S14)

### S11 · 全基因組層證據 17.3 : 1

第一條證據：**全基因組層**。HCC1395 5kHz baseline：
- HP1 reads = 614,000
- HP2 reads = 35,500
- 比例 **17.3 : 1**（94.6% 集中在 HP1）

預期應該是 ~1:1 隨機。這個偏不是某個 chr 特例 — **跨 23 條染色體都一致**出現相同方向的偏。

### S12 · 個別位點層證據

第二條：**個別位點層**。chr19:17565944：
- baseline 113 reads 全部 HP1、HP2 = 0（**113:0 完全失衡**）
- paired 模式方向相反
- V5 修補後翻轉與 paired 一致

這是 IGV 4-BAM 並列圖（baseline / V2b / V3F / V5），S28 會仔細看。

### S13 · 跨樣本一致性 7 / 7

第三條：**跨樣本一致性**。7 個樣本 HP_Ratio 都在 0.91-0.95 區間，**7/7 同方向偏**，Cohen's d = -1.20。

另一個關鍵數字：**同位點跨模式 r = 0.001**（paired vs TO HP_Ratio 幾乎完全不相關），
也就是說 TO HP_Ratio 不反映真實 haplotype，是 self-phasing artifact。
TO-only 標的 LOH，paired 看起來 **86.5% 其實是平衡的**。

### S14 · LOH 的兩個層次（精確化）

回到 S4 提過的兩層次區分。
- **左 ISM HP_Ratio LOH**（read-level）：受 self-phasing 直接影響、**62% 是 artifact**
- **右 LOH.bed**（region-level）：longphase-to-mod 直接輸出、**Jaccard = 1.0 完全不變**

這是 PI 報告 4 的關鍵發現之一：要分清楚這兩個層次，self-phasing 只污染左邊，不污染右邊。

---

## Section 4 · 修補方向 (S15-S19)

### S15 · 修補位置在哪裡？

進入修補。**核心訊息：修補不是 InterSubMod，是 longphase-to-mod**。

路徑 `/big7_disk/liaoyoyo2001/longphase-to-mod/`，**獨立 git repo**，
跟 InterSubMod 是同層平行關係。今天討論的 4 個 commit 全部在 longphase-to-mod 內，
**InterSubMod 本 repo 完全沒有 C++ 改動**。

### S16 · 4-commit 漸進演進

4 個 commit：
- **V2b (8b8c1fd)** — 加 `--pon-only-phasing` flag
- **V3-Fixed (41ff147)** — 重寫 getVote() 兩層 + 修 enum mismatch
- **INDEL guard (380e8d2)** — 補 countINDELHaplotype() UB
- **V5 (working tree)** — Layer 1.5 somatic fallback + alt guard

每個 commit 解不同層次的 bug，**漸進式修補**，每一步都可回退、可驗證。

### S17 · 4 commit 對應 4 條 bug

完整對應表 — 這頁有 4 行表格，每個 commit 對應的 bug 與修補方式都列清楚。
重點：V3F 同時解兩個 bug（priority + enum），所以 4 commits 解了**至少 5 條 bug**。

### S18 · V5 三層投票邏輯

V5 的核心改動是 getVote() 從兩層升級為三層：
- Layer 1：germline first
- **Layer 1.5：somatic fallback（V5 新增）— `else if (somaticHP1>0 || somaticHP2>0)`**
- Layer 2：encode hpResult

V3F 因為缺 Layer 1.5，遇到純 somatic 直接 fallback 為 0；
V5 補上後 **AMB% 從 17.5% 降到 8.0%**、HP:i:33 reads 減 54%。

### S19 · 候選方案排除

動 V5 之前評估過兩個替代方案，都排除了：
1. **替換 LongPhase 為 WhatsHap/HapCUT2**：風險高、未驗證 TO 行為、外部依賴增加
2. **ISM 自己加 haplotag 邏輯**：違反單一職責、ISM F1=0.0124 對 TO germline FP 無修復力

最終選 **surgical fix** — 在 longphase-to-mod 內動 3 個函式。

---

## Section 5 · 改動細節 (S20-S24)

### S20 · 介面契約零變動

V5 改動的關鍵特性：**介面契約零變動**。
HaplotagProcess.h:66-68 三個 method signature **一字未變**。
總修改 **+68/-36 行**，集中於這 3 函式。
好處：上層 caller 不需重編譯、ABI 相容、回退僅需 git revert。

這是 surgical fix 的典範。

### S21 · getVote 兩層 → 三層

具體改動：getVote() 兩層 → 三層。
V5 邏輯：
```
if (germline > 0) { 依 germline }
else if (somaticHP1 > 0 || somaticHP2 > 0) { 依 somatic }   ← V5 Layer 1.5
else { 0 }
```
Layer 1.5 是純分支，**不破壞 germline first 的核心原則**。

### S22 · HP tag 對應表（V5 final）

V5 final 對應表：(germlineResult, somaticTotal>0) → hpResult
- (0, F) → 0
- (1, F) → 1
- (2, F) → 2
- (1, T) → 11
- (2, T) → 21
- (0, T) → 33

V3F 因為缺 Layer 1.5，最後一行 (0, T) 因為 enum mismatch 永不出現；V5 修完後對應表是穩定的。

### S23 · InterSubMod 的位置

**再次強調**職責切分（這頁要花時間講清楚，避免教授誤解）：
- 左 InterSubMod：本 repo，src/core/ 不含 HaplotagProcess，**0 行 C++ 改動**
- 右 longphase-to-mod：獨立 git repo，包含 HaplotagProcess.{h,cpp}，**4 commits × 3 函式**

ISM 只是讀新版 BAM 後**自動受惠**。

### S24 · V5 working tree 未 commit caveat

重要 caveat：V5 = 380e8d2 + **兩塊 working tree 修改（尚未 commit）**：
1. getVote() Layer 1.5 somatic fallback
2. countSNPHaplotype() alt guard 對稱化

**重現性風險**：誰 clone 都要手動 patch。建議切兩個獨立 commits（後續 F1）+ tag v5.0 +
manifest 記 haplotag_version。

---

## Section 6 · 驗證與差異 (S25-S29)

### S25 · 4 項硬性 sanity check

V5 修補完成後做了 4 項硬性檢查：
1. **Δ-consistency 守恆律** PASS
2. **Germline 不變** PASS
3. **Layer 1.5 期望 = 1** PASS
4. **無 germline → HP33** PASS

**15 / 15 PASS、0 violation**。獨立於下游 ISM / F1 結果的硬性內部驗證。

### S26 · Baseline / V3F / V5 量化對比

四個核心數字：
- AMB% **17.5 → 8.0%**
- HP:i:33 reads **-54%**
- 全基因組 clean PS concordance **V5 90.5% vs BL 82.2%（+8.3 pp）**
- 15 sites cherry-picked clean PS concordance **V5 88.2% vs BL 74.9%（+13.3 pp）**

從不同角度都驗證 V5 嚴格優於 baseline。

### S27 · IGV 4-BAM：V5max1 戲劇化

最戲劇化案例。chr19:4639528：
- V3F panel 一大群紫色 HP33 reads（39 條）
- V5 全部正確重分配為 HP11
- **守恆律**：39 條 HP33 → 39 條 HP11，in = out 平衡

這是 Layer 1.5 的具體效果。

### S28 · IGV 4-BAM：SP1 翻轉

chr19:17565944：
- baseline HP1 主導（113:0）
- V2b/V3F/V5 全部翻轉為 HP2 主導
- paired tumor 確認 HP2 是真實方向

3 個 SP-extreme（SP1/SP2/SP3）**全部一致翻轉** — 跨位點一致性證據。

### S29 · F1 為什麼幾乎不變？（重要釐清）

**最關鍵的概念修正**。為什麼 V5 修完 SEQC2 F1 幾乎不變？

- ClairS-TO **raw F1 = 0.7166** 對所有版本完全相同（V5 不改 caller）
- F1 變動只可能來自下游 ISM SuggestFilter
- ISM TO F1 增益本來只有 0.0124
- V5 vs Baseline 最終 F1 = -0.0003 **是統計噪音**

**結論：F1 不能衡量 V5 真實價值**。
真實價值在 read-level tag 品質：clean PS concordance +8.3 pp、AMB% 17.5→8.0%、HP:i:33 -54%。
這些是 read-level 指標，**無法用 F1 反映**。

教授可能會問「F1 沒改善代表沒用嗎？」這張就是回答。

---

## Section 7 · 結論與下一步 (S30-S34)

### S30 · ISM 影響 3-tier 分類

V5 對下游 ISM 特徵的影響分 3 tier：
- **嚴重 38%** — HP-依賴特徵（HP_Ratio、HPSig、HPFineNGroups、HP-stratified pairwise）必重跑
- **中度 7%** — QualityScore 等
- **不受影響 55%** — Pairwise distance（無 HP）、Allele、Caller、甲基化矩陣

重跑工作可聚焦在嚴重 + 中度的 45%。

### S31 · 修正了哪些舊說法

V5 修補也修正了 4 個過去說法：
- HPFineNGroups ≠ methylation bimodality（是 phasing bucket）
- HP:i:21 不必然是當前位點 ALT
- Self-phasing 不污染 LOH.bed
- V5 修補 ≠ caller F1 改善

這些都是過去文件中常出現的說法，需要明確修正。

### S32 · 已知限制 R1 - R9

9 條已知限制：cnLOH 區、AF>0.9 邊界、15 sites cherry-picked、Confidence threshold、
V5 working tree 未 commit、7 樣本擴展未做、L4 部分指標 V5 略遜、problem PS 不適用 read-level、
paired ground truth 自身用 LongPhase。

主要 open issue：**R1（cnLOH）、R5（commit）、R6（7 樣本）**。

### S33 · 後續工作 F1 - F8

8 條後續工作，高優先三項：
- **F1** Commit V5 working tree（高）
- **F3** 7 樣本擴展（高）
- **F8** 100 隨機位點 cross-validate（高）

中優先：F2 vote log unit test、F4 master×flag 重跑 HPFineNGroups、F6 cnLOH、F7 trio-phased。
低優先：F5 manifest 加 haplotag_version。

### S34 · TL;DR 一句話結論

**TL;DR**：
- Self-phasing 在 **LongPhase-TO 階段**以 4-commit 漸進修補
- **InterSubMod 是下游消費者**，本 repo 無 C++ 改動
- **V5 真實價值**：read-level tag 品質 +8.3 pp（clean PS concordance）
- **SEQC2 F1 不變是預期行為** — V5 不改 caller

歡迎提問。詳細數據在 source_materials 三份報告：
1. IGV visual audit
2. V5 audit
3. 技術報告（13-section 結構化版）

Thank you.

---

## 附錄 · 常見 Q&A 預備

| Q | A |
|---|---|
| F1 沒改善代表 V5 沒用嗎？ | 不對。V5 不改 caller，F1 結構上不會變；真實價值在 read-level tag 品質，看 concordance、AMB%、HP:i:33。詳 S29。 |
| 為什麼不直接換 WhatsHap？ | 風險高、外部依賴、TO 行為未驗證、self-phasing 是機制問題換工具不一定解。詳 S19。 |
| 7 個樣本都看過嗎？ | 全基因組層次 7/7 同方向（HP_Ratio Cohen's d=-1.20）；V5 完整 audit 目前只 HCC1395 5kHz，7 樣本擴展是 F3。 |
| 為什麼 InterSubMod 自己不修？ | 違反單一職責，ISM TO F1=0.0124 對 germline FP 無修復力，介面割裂。詳 S19、S23。 |
| 怎麼確認不是 cherry-pick？ | F8（100 隨機位點 cross-validate）就是補強這點；R3 是已知限制。 |
| LOH.bed 真的不變嗎？ | Jaccard = 1.0 完全不變；只有 ISM HP_Ratio LOH（read-level）受影響，62% 是 artifact。詳 S14。 |
| V5 何時可以正式發布？ | F1（commit working tree）+ F3（7 樣本擴展）完成後，目標 tag v5.0。 |
