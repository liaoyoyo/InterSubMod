# Speaker Script v3 — Self-Phasing 完整教授報告（24 slides v3，2026-04-30）

> 演講者用稿。每張平均 75 秒；總長 30 分鐘。整合 13 條教授可能 Q&A 預備答案（含 v3 新增 Q12 PON 雙路徑、Q13 V5 mid-low purity）。
> **v3 6 處口徑校準（嚴守）**：
> 1. V5 不修 self-phasing 本身（V2b 已處理 phase scaffold）。
> 2. 「同實驗室相鄰工作」非「業界共識」「標準替代」。
> 3. ISM 影響只列 count（29/14/42 共 85），來源報告誤寫 7%。
> 4. `Util.h:20-25`（不是 21-25）。
> 5. +8.3 pp（clean PS blocks，全基因組），不是全基因組 raw。
> 6. **★ v3 新增**：baseline 投票邏輯 = somatic-first 投票順序（不是「有 HP11/HP21 即 somatic」這個簡化敘述）；root cause 觸發條件 = baseline 估 purity = 0.927 ≤ 0.95 未觸發 Two-Pass → 走三條路（不是「baseline 誤判 purity」）。

---

## Section 1 · 高層次重點（S1-S2）

### S1 — 17.3 : 1 → 1 : 1（衝擊式 TL;DR + 動機小字）

開場 60 秒鎖定核心結論。本 PPT 是 v3 24-slide 版本，比 v2-26 更精簡。三行核心訊息：

1. **大字 17.3 : 1 → 1 : 1**：HP1 family : HP2 family read 比例。Baseline 17.3 倍偏在 HP1，修補後回到生物學期望 1:1。HP1 family = HP1 germline reads + HP1_1 somatic-on-HP1 reads，HP2 family 同理。HCC1395 5kHz 全基因組 baseline HP1 family = 614,000 reads，HP2 family = 35,500 reads，比值 17.3，跨 23 條染色體一致，94.6% 集中在 HP1。
2. **read-level concordance +8.3 pp**：全基因組 clean PS blocks 對 paired ground truth 的吻合度提升（V5 90.5% vs Baseline 82.2%，+8.3 pp）。**口徑校準必加 caveat「(clean PS blocks，全基因組)」**。
3. **重要分層**：self-phasing 問題鏈是被分層處理的。V2b commit `8b8c1fd` 解 phase scaffold（啟用 `--pon-only-phasing` flag），V3F commit `41ff147` + V5 working tree 解 tag 層（priority bug + enum mismatch + V3F 過度保守）。**V5 不修 self-phasing 本身**（V2b 已處理 phase scaffold）。

底部小字（v3 新增動機 strip）：「修補 self-phasing 是解鎖 InterSubMod 五大研究目標的前提」。InterSubMod 五大目標包括：目標 1 per-CpG 多標籤關聯、目標 2 clone 結構、目標 3 二次打擊、目標 4 TO normal 補強、目標 5 F1 panel；目標 1/2/4 直接依賴 HP tag 正確。

本 PPT 24 slides 分 6 段：高層次重點 (S1-S2) → 觀察問題 (S3-S7) → 為何重要 (S8-S9) → 解釋與原因 (S10-S13) → 改動驗證 (S14-S20) → 未來規劃結語 (S21-S24)。本工作焦點在 longphase-to-mod fork（獨立 git repo @ `/big7_disk/liaoyoyo2001/longphase-to-mod/`），InterSubMod 本 repo 無 C++ 改動。

### S2 — Pipeline 與 4 階段工作

本張展示三個系統的 pipeline 切分與本工作 4 階段的關係。

**Pipeline**：tumor BAM → ClairS-TO (caller) → longphase-to-mod (phasing + haplotag, V5 修補位置) → tumor_tagged.bam (HP:i tag) → InterSubMod (read-level ISM)。

**強調**：InterSubMod 不是修補位置，本 repo 無 C++ 改動。HaplotagProcess.cpp 在 longphase-to-mod fork（獨立 git repo），不在 ISM 內。ISM 只是讀新版 BAM 後下游分析自動受惠。

**4 階段工作**：(1) 機制定位 — 從理論預期 1:1 對照觀察到的 17.3:1，推導 self-phasing 機制；(2) 4-commit 修補 — V2b → V3F → INDEL guard → V5，每個 commit 解一層 bug；(3) 驗證 — 4 項硬性 sanity check 15/15 PASS + 5 條獨立證據鏈；(4) 影響評估 — ISM 3-tier (29/14/42 共 85) 與五大研究目標銜接。

後續 PPT 段落依此 4 階段順序展開。v2-26 原 S2 業界家族樹已吸收進 v3 S20。

---

## Section 2 · 觀察到問題（S3-S7）

### S3 — HP tag 定義（5 整數值 + PS / LOH 兩層）

建立 HP tag 認知。BAM 內每條 read 一個整數 HP tag，6 種可能值：HP:i:0 (unphased)、HP:i:1 (germline HP1)、HP:i:2 (germline HP2)、HP:i:11 (somatic on HP1)、HP:i:21 (somatic on HP2)、HP:i:33 (somatic ambiguous)。

這個編碼是 longphase-to-mod 自訂，不在 LongPhase 官方 spec，與 C++ enum (`Util.h:20-25`) HAPLOTYPE_UNDEFINED=-1, HAPLOTYPE1=1, HAPLOTYPE2=2, HAPLOTYPE1_1=3, HAPLOTYPE2_1=4, HAPLOTYPE3=5 不同 — 這個型別語意失配是後面 S11 enum bug (Bug 2) 的伏筆。

PS = phase-set ID，同一個 PS 表示這些 reads/variants 屬於同一個 phase block。

LOH 兩層必須區分（後面 S6+S7 詳述）：
- **ISM HP_Ratio LOH**：read-level，從 BAM HP tag 統計（受 self-phasing 影響，62% artifact）
- **LOH.bed**：region-level，從 VCF AD 偵測（不受 self-phasing 影響，V5 前後 Jaccard=1.0）

### S4 — 證據三層（理論 / 全基因組 / 個別位點）

證據三層並列（v3 合併原 v2-26 S5 觀察 17.3:1 + S6 理論 1:1 為單一 slide）。

**第一層 理論**：為什麼 17.3:1 是 artifact 而非真實生物？Germline het 是個體 zygosity 層的差異，在染色體上隨機分散，沒有方向性偏好；Long-read 平均 5-30 kb 跨多個 phase block，理論上 HP1 / HP2 應各分到一半 reads。**正確預期 HP1 : HP2 = 1 : 1**。

**第二層 全基因組**：HCC1395 5kHz baseline 硬數據：HP1 family = 614,000，HP2 family = 35,500，比例 17.3:1；94.6% 集中於 HP1；跨 23 染色體一致。fig01d Panel D 的 6 圖證明這是全基因組層的系統性偏差，不是個別 hotspot 也不是統計噪音。

**第三層 個別位點 (SP1 chr19:17565944)**：baseline HP2:HP1 = 113:0 完全失衡（113 reads 全 HP1，HP2 = 0）；paired tumor BAM 確認真實方向是 HP2（110:2）— baseline 完全反向。

三層獨立證據彙整於同一結論：self-phasing artifact 真實存在，不可能是隨機誤差。

### S5 — Paired 對照確實 ~1 : 1

核心對照證據。同樣本 paired 與 TO 的 HP_Ratio 對照：Paired (有 normal) 全 7 樣本 HP_Ratio 在 0.49-0.51 之間，接近預期 0.5；TO baseline (無 normal) 全 7 樣本 HP_Ratio 在 0.91-0.95，平均 0.94，**Cohen's d = -1.20** 巨大效應量。

關鍵統計：同位點跨模式 (paired vs TO) HP_Ratio 相關係數 r = 0.001，n = 288K pairs — 幾乎完全無相關，TO HP_Ratio 完全不反映真實 haplotype，是 self-phasing artifact。TO-only 標記為 LOH 的位點，paired 模式下 86.5% 是平衡的，TO 的 LOH 標記大部分是 self-phasing artifact 而非真實 LOH。

**Q2 預備**「為什麼不直接用 paired？」：normal BAM 不一定可得（臨床、archive、cell-free DNA），TO 是必要研究方向。本實作目標就是讓 TO 模式達到接近 paired 的 tag 品質。**機制**：paired normal 提供 germline het 的 ground truth scaffold；TO 模式無 normal，phasing 必須從 tumor 自己推 germline，但 tumor 含 somatic → somatic 反客為主進 phasing graph → self-phasing。

### S6 — 拆鎖：phasing 沒問題，是 tag 的問題

拆鎖關鍵：self-phasing 影響的是哪一層？答案是 tag 層，不是 phasing 層（v3 從原 v2-26 S8 拆出）。

**左欄 Phasing 層 (LOH.bed)**：路徑 = VCF allele depth (AD) → LongPhase region 偵測，輸出 BED region-level LOH。V5 前後 Jaccard = 1.0 完全不變，V2b PON-only 已解 phase scaffold self-phasing。

**右欄 Tag 層 (BAM HP_Ratio)**：路徑 = BAM HP:i tag → ISM per-region HP_Ratio 統計。Baseline 17.3:1 嚴重偏差，62% ISM HP_Ratio LOH 是 artifact，V5 修補對象就是這層。

兩套 LOH 系統使用不同數據源 (VCF AD vs BAM HP tag)，kappa = 0.670 由此解釋。**Q3 預備**「Jaccard=1.0 怎麼證明？」：baseline LOH.bed 與 V2b PON-only LOH.bed 的 region-level Jaccard 相似度 = 1.0，兩 BED 完全相同。本工作影響的是 BAM HP tag → ISM HP_Ratio LOH，不是 LOH.bed。

### S7 — LOH 兩層精確化（kappa = 0.670 解釋）

LOH 兩層精確細節（v3 從 S8 拆出，避免常見誤解）。

**ISM HP_Ratio LOH (左欄)**：數據源 BAM HP:i tag；計算 HP_Ratio < 0.1 or > 0.9 標 LOH。受 self-phasing 影響嚴重：62% 是 artifact，AF 0.1-0.8 範圍近 100% artifact，TO TP 中 86.5% paired 平衡。本工作修補後改善。

**LOH.bed (右欄)**：數據源 VCF allele depth (AD)；計算 LongPhase region 偵測。不受 self-phasing 影響：V5 前後 Jaccard = 1.0。本工作不影響此層。

**kappa = 0.670 解釋**：兩套系統用不同數據源，本就期待中度一致性 (kappa 0.4-0.7)；0.670 不是 bug，而是兩個獨立衡量的自然差異。本工作 self-phasing 修補只影響 ISM HP_Ratio LOH（左欄），不改變 LOH.bed（右欄）。

教授若混淆兩層，會誤以為「修補後 LOH.bed 也變」，這是錯的。

---

## Section 3 · 為何重要（S8-S9）

### S8 — TO 根基 + ISM 強依賴 HP tag

為何重要的學術意義。**左欄 TO 模式為何是研究根基**：(1) 臨床、archive、cell-free DNA、舊樣本等情境下 normal BAM 不可得，TO 是必要研究方向；(2) TO 純數據是 TP/FP 觀察的根基，不能依賴 paired 比對外推；(3) 本實作目標就是讓 TO + PON 條件下 tag 品質接近 paired 水準。

**右欄 ISM 4-bucket 分群核心**：ISM 把 reads 分成 HP1 (germline HP1)、HP1-1 (somatic on HP1)、HP2 (germline HP2)、HP2-1 (somatic on HP2) 四個 bucket，計算 region-level 特徵。依賴 4-bucket 的關鍵特徵：NGroups (subclone count)、HPSig (HP signature)、HP_Ratio_norm、HPMergedDelta、HPFineNGroups (subclone marker)。

如果 HP tag 偏差 (17.3:1)，4-bucket 分群錯誤，所有依賴特徵的結論都不可信。整個 ISM 研究都建立在這個基礎上，因此 self-phasing 修補是研究可信度的前提。

### S9 — 影響範圍 29 / 14 / 42（共 85）

V5 修補對下游 ISM 特徵的影響量化為 3 tier（來源 source 03 §4.3）：

- **嚴重影響 (HP-依賴) 29 個**：HP_Ratio (62% 假 LOH)、Potential_LOH、HPMergedDelta/Sig (方向反轉)、HPFineNGroups (含 NG=2 LOH-constrained phasing)。**必須在 V5 BAM 上重跑才能信任**，舊結論需加註版本。
- **中度影響 (間接污染) 14 個**：QualityScore (AUC 0.497)、GlobalP (取 HP/Allele 最小值)、CramersV (取 HP/Allele 最大值)、VerificationClass (label_sig 含 HP 成分)。重評多數影響微弱。
- **不受影響 (無 HP 依賴) 42 個**：PairwiseMean/MedianDist、AlleleDelta/AlleleP、Caller 特徵、甲基化矩陣、CpG 座標、region_methyl_mean。**結論穩固不需重測**。

**口徑校準關鍵**：slide 上只列 count 不列 %。來源報告寫 14/85 = 7% 但 14/85 實際 16.5%，來源百分比有誤。為避免傳遞錯誤數字，全部以 count 表達。重跑工作可聚焦在嚴重 + 中度的 43 個，不需全 85 重跑。

---

## Section 4 · 解釋問題與發生原因（S10-S13）

### S10 — Prerequisite（4 函數 + PON + Purity 觸發鎖）

v3 合併原 v2-26 S11 (4 函數表) + S12 (PON 概念) 為單一 prerequisite slide，並補上 Purity 0.95 閾值概念。

**左半 longphase-to-mod 4 函數**：(1) `getVote()` — 統計 reads 對 HP1/HP2/somatic 的投票，V3F + V5 修補主戰場（HaplotagProcess.cpp:512-563）；(2) `judgeHaplotype()` — vote → hpResult integer，V3F enum 修在這（line 697）；(3) `countSNPHaplotype()` — SNP 位點各 hap read 計數，V5 補 alt path UNDEFINED guard（line 489-494）；(4) `countINDELHaplotype()` — INDEL 版本，380e8d2 commit 補 UB guard（line 497-510）。

**右半 PON 概念**：PON (Panel of Normals) 是群體 germline 變異資料庫（1000 Genomes、CoLoRSdb、dbSNP、gnomAD）；TO 模式無 normal BAM 時 PON 替代 normal 作 phasing anchor。V2b 啟用 `--pon-only-phasing` 解 phase scaffold self-phasing。

**★ Purity 0.95 觸發鎖（v3 補充，S11 主者預告）**：`PhasingProcess.cpp:197` 硬編碼 `purity > 0.95` 閾值，決定是否啟動 baseline 內建的 Two-Pass 路徑。HCC1395 5kHz 真實 purity > 0.95，但 baseline 估計 0.927（四捨五入 0.93，**非誤判，是邊緣 case**），0.927 ≤ 0.95 → 未觸發 Two-Pass → 走「三條路」標準流程。

**Q6 預備**「PON 有資料還會 self-phasing？」：PON classification (VCF 層) 與 phasing anchor (read 連結層) 是獨立步驟。Baseline bug：PON 雖標出 somatic 但 phasing 階段仍把 somatic 當 anchor。V2b 啟用後 phasing 才只用 PON-confirmed germline 做 anchor。「修一個 bug 暴露另兩個 bug」：PON-only 啟用後 germline votes 急減，原本被 paired germline 訊號掩蓋的 tag 層 priority bug + enum bug 立刻顯形。

### S11 — ★ Baseline 根因樹（觸發條件主者 + 三 bug）

**v3 ★ KEY slide** — 加入 Purity 0.95 觸發條件主者作為 root-cause tree 的根。Self-phasing 不是單一 bug，而是「特定觸發條件下三層獨立故障的集中暴露」。

**根（觸發條件主者，紅色框）**：Baseline 估 purity = 0.927 ≤ 0.95（純樣本標準）→ 未觸發 Two-Pass → 走 baseline 標準三條路流程。

**口徑校準關鍵**：不是「baseline 誤判 purity」（0.927 是正確估計，是邊緣 case），而是 baseline 設計上的硬編碼閾值 0.95 對近似純樣本不夠保守。HCC1395 5kHz 真實 purity > 0.95（純 tumor，無 normal contamination），但 baseline 的 polynomial calculator 估出 0.927，這個 0.927 數值本身正確，只是落在 0.95 閾值的邊緣下方就觸發了 baseline 標準流程而非 Two-Pass 路徑。

**三葉（bug）**：

- **Bug 1（紅色，getVote priority）**：variantKeys 順序 {HP1_1, HP2_1} 排第一，任一 somatic count > 0 即 break，germline 完全跳過。為何 paired 沒事？paired germline votes 多但 priority bug 仍存在；TO PON-only 後 germline votes 急減立刻顯形（HCC1395 99.9% reads 拿 HP21）。**V3F commit 41ff147 修**。
- **Bug 2（橘色，enum vs int literal mismatch）**：C++ enum 定義在 `Util.h:20-25`（注意是 20-25 不是 21-25，HAPLOTYPE_UNDEFINED=-1 起於 line 20），HAPLOTYPE1_1 = 3。但 BAM HP tag integer = 11。caller 端 `if(hpResult != HAPLOTYPE1_1)` 用 enum=3 比較 hpResult=11/21/33，fallback 死分支永不執行，HP:i:33 永不出現（baseline 0/15 sites）。C++ 把 enum 隱式轉 int 不會 warn。**V3F 同 commit 改為 integer literal 比較**。
- **Bug 3（灰色，scaffold，已由 V2b 解）**：`PhasingProcess.cpp:154-157 convertNonGermlineToSomatic()` 不被觸發，somatic 當 phasing anchor。**已由 V2b commit 8b8c1fd 解決**（透過 `--pon-only-phasing` flag），本 PPT 焦點不在此，視覺只佔 20%。

**因樣本實為純無 normal contamination 而流程假設有 → 暴露 tag 層 somatic-first 投票 bias → 結果 self-phasing 17.3:1 artifact**。驗證來源：`source_materials/04_purity_calculator_failure_root_cause.md`（v5_audit_suite/18 複製）。

### S12 — ★ V5 三層投票流程（細化定位）

**v3 ★ 強化版** — speaker note 整合 Q2 答（V5 在 mid-low purity 比 baseline 更好）。V5 對 `getVote()` 的具體改動：從 V3F 兩層升級為三層投票邏輯。

**Layer 1 (germline first)**：`if (g1 > g2) → 1; elif (g2 > g1) → 2; else → 進入 Layer 1.5`。**這裡是 baseline 的 somatic-first 投票順序的關鍵翻轉**：baseline `getVote()` 用 variantKeys loop 把 {HP1_1, HP2_1} 排第一，somatic 投票優先；V5 直接看 germline 投票數，翻轉投票優先序。

**Layer 1.5 (somatic fallback, V5 新增, mid-low purity 防禦層)**：當 germline 平手時參考 somatic 投票。`if (sTotal == 0) → 0; elif (s1 > 0.6*tot) → align HP1; elif (s2 > 0.6*tot) → align HP2; else → 33 (true ambiguous)`。Confidence threshold 0.6 是經驗值，**設計上是 mid-low purity 的防禦層** — 在 0.6-0.93 中等純度範圍，弱 directional 訊號標 33 ambiguous（保守 +10pp HP33 比例）避免錯誤分配。

**Layer 2 (encode hpResult)**：`germline=1 + somatic_total>0 → HP:i:11; germline=2 + somatic_total>0 → HP:i:21; germline=0 + somatic_total>0 → HP:i:33; germline only → HP:i:1 or 2; nothing → HP:i:0`。

**★ Q2 答（V5 vs baseline 在 mid-low purity）**：V5 更好。

- 純樣本（HCC1395 5kHz）：baseline somatic-first 在 PON-only 啟用後 germline votes 急減 → somatic 主導 → 17.3:1 artifact；V5 germline-first 不受此影響。
- **Mid-low purity 場景（0.6 sample, t30_n20）**：baseline 仍 somatic-first 但 self-phasing 自然弱化（normal contamination 平衡 hap）；V5 加入 confidence 0.6 threshold → 對「弱 directional 訊號」標 HP33 ambiguous → 09 報告：0.6 sample HP33 比例 12.4% vs Baseline 2%（+10pp 保守）→ 避免錯誤分配。

效果：AMB% 17.5% → 8.0% (-9.5 pp)；HP:i:33 reads 239,679 → 110,197 (-54%)；介面契約零變動（HaplotagProcess.h:66-68 一字未變，drop-in replace）。

Caveat：V5 在 SP-extreme self-phasing 位點仍不修（V5 設計就不修 phase 層問題，那是 V2b PON-only 的責任）。

### S13 — 程式碼 diff（Baseline vs V5）

**v3 修正**：baseline 欄上方加紅標 caveat「投票邏輯 = somatic-first 優先序（非「有 HP11/HP21 即 somatic」這個簡化敘述）」。

**左欄 Baseline (有 bug)**：
```cpp
for (auto& vk : variantKeys) {
  // {HP1_1, HP2_1} listed first
  if (count[vk.first] > 0
      || count[vk.second] > 0)
    break;  // somatic preempts
}
if (hpResult != HAPLOTYPE1_1)
  hpResult = 0;  // enum=3
  // never matches HP tag=11
```

**★ V3 紅標 caveat**：baseline 投票邏輯不是「有 HP11/HP21 即 somatic」這個簡化敘述，而是「somatic-first 投票優先序設計缺陷」 — variantKeys 排序把 somatic key 排在前，任一 somatic 投票觸發 break，結果 germline 投票被忽略。

這個區分對教授很重要：簡化敘述會讓人誤以為「baseline 一旦看到 somatic tag 就分類為 somatic」，但實際是「投票順序決定優先權」，V5 修法是「翻轉投票優先序為 germline-first」，不是「忽略 somatic」。

**右欄 V5 (修補後)**：
```cpp
// Layer 1: germline first
if (g1 > g2) germlineRes = 1;
else if (g2 > g1) germlineRes = 2;
// Layer 1.5: somatic fallback
else if (s1 + s2 > 0) {
  if (s1 > 0.6*tot) align HP1;
  else if (s2 > 0.6*tot) align HP2;
}
if (hpResult != 11)  // int literal
```

**三處紅標對應前面 S11 root-cause tree**：(1) priority loop 拆解（somatic-first → germline-first 翻轉，Bug 1）；(2) Layer 1.5 純分支新增（V5 新工作，mid-low purity 防禦層）；(3) enum HAPLOTYPE1_1 改為 integer literal 11（Bug 2）。

**Q13 預備**「程式碼 diff 看不太懂」：請看這三處紅標，每行對應前面 S11 一條 bug。完整程式碼參考見 `notes/code_references.md`，V5 commit hash 待 working tree commit 後補（R5）。

---

## Section 5 · 改動驗證與結論（S14-S20）

### S14 — ★ Sanity 4 + 5 證據鏈（PON 雙路徑有效性整合）

**v3 ★ 強化版** — speaker note 整合 Q1 答（PON-only Two-Pass 在 0.93 與 0.6 sample 都有效）。

**左欄 4 項硬性 sanity check**（HCC1395 5kHz 15 sites，v5_audit_suite/06 報告）：

1. 守恆律 A · Δ-consistency：ΔHP33 + (ΔHP11 + ΔHP21) = 0，tag 移轉量平衡，**15/15 PASS**。
2. 守恆律 B · Germline 不變：V3F 與 V5 的 HP1/HP2 reads 逐 site 比對，0 reads 漂移，**15/15 PASS**。
3. Layer 1.5 期望 1 · 33→directional 精確守恆：ΔHP11 == n(V3F=33→V5=11)，**15/15 PASS**（V5max1 chr19:4639528 39 reads 精確守恆）。
4. Layer 1.5 期望 2 · 無 germline → HP33：跨 15 sites pool transition table，**0 violation**。

**右欄 5 條獨立證據鏈**：

1. 理論層 — germline het 隨機 → 預期 1:1。
2. 全基因組層 — HP1:HP2 = 17.3:1 跨 23 chr 一致。
3. 個別位點層 — SP1 chr19:17565944 baseline 113:0 全失衡。
4. Sanity check 4 項 15/15 PASS（左欄）。
5. 程式碼最小必要 — +68/-36 行集中於 3 函式（getVote / countSNPHaplotype / countINDELHaplotype）。

5 條獨立路徑同步收斂於同一結論：self-phasing 修補必要且充分。

**★ Q1 答（PON-only Two-Pass 在不同 purity 下都有效）**：

- **0.93 純樣本（HCC1395 5kHz）**：PON-only Two-Pass 可執行，sanity 15/15 PASS，Aggregate paired GT concordance V5 78.85% vs Baseline 72.20% (+6.65pp)，Clean PS paired GT V5 88.2% vs Baseline 74.9% (+13.3pp)，全基因組 clean PS V5 90.5% vs BL 82.2% (+8.3pp)。
- **0.6 sample (t30_n20)**：PON-only Two-Pass 可執行，09 報告 V5 HP33 比例 12.4% vs Baseline 2% (+10pp 保守 tagging) → 在 mid-low purity 也有改善，設計上的 conservative tagging 防禦層。

**Caveat**：0.6 sample 沒有 paired reference（synthetic mix 無 ground truth）無法直接做 concordance，但 HP33 比例分析顯示 V5 確實在 mid-low purity 提供保守 tagging 防禦。

**Q7 預備**「sanity check 是誰寫的？覆蓋率多少？」：本工作新增（v5_audit_suite/06 Agent D），覆蓋限制是 15 sites cherry-picked (R3)，需 7 樣本擴展 (F3) + 100 隨機位點 cross-validate (F8)；整體 sanity check 必要不充分。

### S15 — 量化指標（AMB / HP33 / Concordance）

v3 從原 v2-26 S17 拆出，專注量化指標。下一張 S16 才講 4-commit timeline。

量化指標彙整：(1) AMB% (HP:i:33 比例) 從 17.5% 降到 8.0% (-9.5 pp)，V5 解了過半 ambiguous reads；(2) HP:i:33 reads 數從 239,679 降到 110,197，減 54%；(3) 全基因組 clean PS blocks paired GT concordance：V5 = 90.5% vs Baseline = 82.2%，**+8.3 pp**；(4) 15 sites cherry-picked clean PS concordance：V5 = 88.2% vs BL = 74.9%，+13.3 pp。

**口徑校準關鍵**：+8.3 pp 必須加 caveat「(clean PS blocks，全基因組)」，不是全基因組 raw，不是 cherry-picked。

**Q8 預備**「N 多少？信賴區間？跨樣本一致？」：15-site cherry-picked N=11 clean PS；全基因組 clean PS N 在數萬-數十萬 reads 量級（具體 N 在 PI 報告 4 §3.7）；跨樣本擴展未做（F3 待辦）；信賴區間 source 03 未直接給，推估 Wilson 95% CI ±1 pp 內。Caveat：clean PS blocks 篩選 germline accuracy ≥ 70%；problem PS blocks 上 V5 與 Baseline 接近隨機。

### S16 — 4-commit 漸進演進 timeline

v3 從原 v2-26 S17 拆出，專注 4-commit timeline 視覺化。

4-commit 漸進演進（路徑：longphase-to-mod fork @ `/big7_disk/liaoyoyo2001/longphase-to-mod/`，獨立 git repo）：

1. **V2b commit `8b8c1fd`**（PON-only）— 解 phase 層 self-phasing scaffold，Phasing.cpp +9/-2、PhasingGraph.cpp +34/-0、PhasingProcess.cpp +25/-3；HaplotagProcess.cpp 未動。啟用 `--pon-only-phasing` flag。
2. **V3F commit `41ff147`**（priority + enum）— 解 tag 層 Bug 1 + Bug 2，HaplotagProcess.cpp:506-541 (getVote 重寫) + :697 (caller 端)；+36/-25。
3. **INDEL guard commit `380e8d2`** — 補 UB 漏洞（countINDELHaplotype），HaplotagProcess.cpp:497-510；+8/-4。
4. **V5 working tree (未 commit)** — 解 V3F 過於保守，HaplotagProcess.cpp:489-494 + :512-563；+24/-7。

累計：+68/-36 行集中於 3 函式（getVote / countSNPHaplotype / countINDELHaplotype）；介面契約 HaplotagProcess.h:66-68 一字未變（drop-in replace）。

**V5 working tree 未 commit (R5 caveat)**，後續 P0 行動 F1 切 2 獨立 commits 完成 + tag v5.0。

### S17 — 🎯 CLIMAX：V5max1（39 reads, 100% reassigned）

演講 CLIMAX。chr19:4639528 V5max1 是最戲劇化的視覺證據。

V3-Fixed panel 看到一大群紫色 reads 標為 HP:i:33 (ambiguous)，39 條；V5 修補後這 39 條全部正確重分配為 HP:i:11 (somatic on HP1)。守恆律驗證：V3F 39 reads HP33 → V5 39 reads HP11，in = out 完全平衡，sanity check PASS。這是 Layer 1.5 somatic fallback 的具體效果 — 原本因為 V3F 過於保守卡在 ambiguous bucket 的 reads，現在透過 confidence 0.6 threshold 被正確歸到具體 hap。

**底部 caveat 防誤解（reviewer R2 採納）**：「但 V5 不修 self-phasing 本身（V2b 已處理 phase scaffold）」。

**口徑校準關鍵**：高潮圖容易讓教授誤以為 V5 全面解決 self-phasing，但 self-phasing 問題鏈是被分層處理的：phase 層的 self-phasing scaffold 由 V2b PON-only phasing 處理；V3F/V5 解的是 tag 層。V5 的「39 reads → 100% reassigned」是 directional reassignment 不是 phasing graph 重建。speaker 必須在 60-90 sec 內一句說清這個分層。

### S18 — HP 翻轉證據（Bar chart + IGV 縮圖）

S17 climax 用單一 IGV 圖視覺單調；S18 用 bar chart 加上 IGV 縮圖打破單調（climax 延伸）。

Bar chart：3 個 SP-extreme 位點，每位點 3 組 bar（Baseline / V5 / Paired）。SP1 chr19:17565944 baseline 113:0 → V5/V2b 0:113 → paired 2:110 確認 HP2。SP2 chr19:12452332 baseline 109:1 → V5 1:109 → paired 3:105。SP3 chr19:12467180 baseline 108:0 → V5 0:108 → paired 1:100。3/3 SP-extreme 一致翻轉至 paired 方向。

**重要**：為什麼 baseline orientation 跟 paired 完全相反？因為 self-phasing 把 reads 歸到錯邊（sub-clone somatic 集合反客為主），修補後 phasing 正確以 germline 為 anchor，reads 回到正確 hap。**SP1-3 修補主要是 V2b 主導（phase 層），V5 在這些位點不再變動（Δ ≈ 0）；V5 主戰場是 V3F 過度保守的 33 → 11/21（S17 V5max1）**。

### S19 — 矛盾或盲點（F1 釐清 + cnLOH 邊界）

教授必問（**Q9**）：F1 不變是預期 — 那為什麼還要做這修補？對發論文有意義嗎？

**為什麼 F1 不變**：ClairS-TO raw F1 = 0.7166 對所有版本（Baseline / V3F / V5）完全相同，因為 V5 不改 caller，改的是 BAM HP tag 編碼，與 caller 輸出 VCF 無關。F1 變動只可能來自下游 ISM SuggestFilter（RegionProcessor.cpp:1120, 1269），但 ISM TO F1 增益本來只有 0.0124（對 caller germline FP 修復力先天不足），V5 vs Baseline 最終 F1 差 -0.0003 噪音。**結論**：F1 不能衡量本實作品質，真實價值在 read-level concordance +8.3 pp。

**為什麼仍要做**：(1) 下游 ISM 影響 — 29 個 HP-依賴特徵必須在 V5 BAM 上才可信；(2) 生物學詮釋正確性 — HPFineNGroups marker 過去解讀為 methylation bimodality，V5 修補後重詮釋為 phasing bucket signature（Thread D 主軸候選）；(3) 跨樣本研究可信度 — 7/7 樣本同方向影響全部研究方向（五大目標 1/2/4）。

**對發論文意義**：本工作不是新 caller / 新 phaser，是填補 LongPhase-TO 在 PON-only 啟用後的 tag 層 bug；與 LongPhase-S 2025 paired 形成姐妹工作；真實 contribution 是 read-level tag 品質 +8.3 pp，對下游 epigenetic analysis 是 enabling step；學術定位是 InterSubMod 五大目標的解鎖前提。

**Caveat 4 條**：cnLOH 區未獨立評估 (R1)；AF > 0.9 邊界 (R2)；Confidence threshold 0.6 未直接驗證 (R4)；V5 working tree 未 commit (R5)。

### S20 — 業界家族樹 + 2x2 比較表（合併）

v3 合併 — 上半業界家族樹（原 v2-26 S2）+ 下半 2x2 比較表（原 v2-26 S21）為單一 slide。

**上半業界家族樹**：上游基底是 LongPhase（Lin 2022 Bioinformatics）。兩個分支：(a) LongPhase-S（bioRxiv 2025.11.20.689492v1）是同實驗室的 paired 版，ClairS SNV F1 +4.5%；(b) longphase-to-mod V5（本實作）是 TO + PON 條件下的本地實作。

**下半 2x2 matrix (Source × Mode)**：公開工具 + Germline = LongPhase / WhatsHap / HapCUT2；公開工具 + Tumor = LongPhase-S（Paired only）；本實作 + Germline = 不做 germline；**本實作 + Tumor = longphase-to-mod V5（Tumor-Only + PON，+8.3 pp clean PS concordance，TO + PON 條件下的本地實作，填補 TO 場景 gap）**。

**口徑校準關鍵**：本實作填補 tumor-only 在公開工具中的 gap（WhatsHap 不支援 TO），**不是「業界共識」「標準替代」，而是「同實驗室相鄰工作」**。LongPhase-S +4.5% 為 paired 場景參考（注意非直接可比）。

**Q10 預備**「跟 LongPhase-S 重疊多少？算 contribution 嗎？」：LongPhase-S 是 paired，本實作是 tumor-only；設計哲學一致，本工作填補 TO gap。Contribution 區隔：4-commit 漸進修補在 longphase-to-mod fork（獨立 git repo），4 項 sanity check 是新貢獻，ISM 下游影響 3-tier 分類是新貢獻，ISM 跨樣本一致性 7/7 + Cohen's d=-1.20 是新貢獻。

---

## Section 6 · 未來目標與規劃 + 結語（S21-S24）

### S21 — 後續可動 + 已初步現方向（合併）

v3 合併 — 上半「後續可動」(原 S22 4 row) + 下半「已初步現」(原 S24 3 cards) 為單一 slide。

**上半 4 row 可一起做的事**：
1. **Phase 2A — Normal methylation reference**：依賴 V5 BAM；HP tag 正確後 normal-only methylation reference 可用，對應五大目標 4 解鎖。
2. **HPFineNGroups marker — master × flag 重驗**：subclone marker ⭐4→⭐3 已降級因 HCC1395 TO ClairS-TO raw split 無法重現 master 89.1% (Fisher odds 0.913 反向 p=3.5e-3) + flag=on NG≥3=0；機制重詮釋為 phasing signature。
3. **Archive haplotag_version 標記**：manifest.yaml 加 haplotag_version 欄位（baseline / V5），對應 P0 行動 F5。
4. **Thread D — LOH-constrained phasing（論文主軸候選）**：NG=2 cross-sample 6/6 POSITIVE (Wilcoxon p=0.0156, median gap 0.365)，Inner same-hap 93-99% (obs18 6 樣本)；需獨立 phasing-vs-methylation 驗證；V5 BAM 是前提。

**下半 3 卡片 已初步現**：
- (a) **Thread D NG=2 6/6 POSITIVE**：B1 報告 2026-04-23 Wilcoxon signed-rank W=21（6 樣本配對最強顯著性），exact p=0.0156；配 paired control (B3) p=0.578 不顯著 → 確認 LOH-constrained phasing 是 TO-specific；HCC1954 outlier (B2) 解析為 caller FP 背景 84% 非 phasing failure。
- (b) **HPFineNGroups 重驗中**：評級 ⭐4→⭐3；Phase 2B 計畫 master × flag 機制本質判定。
- (c) **LOH-constrained phasing 論文主軸候選**：Inner same-hap 93-99%，需獨立驗證。

**Q11 預備**「Thread D 6 樣本可信嗎？」：Wilcoxon 6 樣本配對最強顯著性是 6/6 同方向 W=21 p=0.0156，本實驗達此上限；evidence grade B（待 Phase 2B 升 A）；多重比較校正未做（只 1 假設不需 Bonferroni）。

### S22 — 五大目標銜接 + 研究發展樹

InterSubMod 五大研究目標（來源 `InterSubMod/docs/architecture/20260327_InterSubMod研究願景定錨_01.md`）：

- **目標 1** — per-CpG 甲基位點多標籤關聯性評分。
- **目標 2** — clone 結構分析（sub-clone + 共演化）。
- **目標 3** — 二次打擊與事件順序推論。
- **目標 4** — TO normal 資訊補強。
- **目標 5** — 整合 evidence panel 提升 F1。

**研究發展樹**：self-phasing 修補完成 (V5) → 4-bucket 分群可信 → NGroups/HPSig/HP_Ratio 可信 → 五大目標解鎖。目標 1/2/4 直接依賴 HP tag 正確；目標 3/5 部分依賴。

**Self-phasing 修補是 5 大目標解鎖前提**，特別是目標 1/2/4。本張展示研究藍圖，讓教授看到本工作不只是一次性修補，而是支撐整個 InterSubMod 的根基。與 S1 底部小字呼應。

### S23 — 短期 P0 行動清單

短期 P0 行動只列 3 條高優先，其餘 F5-F8 入 speaker note。

- **F1 (高) Commit V5 working tree**：切 2 獨立 commits — Commit A getVote() Layer 1.5 somatic fallback（行 512-563）；Commit B countSNPHaplotype() alt guard（行 489-494）。完成後 tag v5.0，在 InterSubMod manifest.yaml 加 haplotag_version。
- **F3 (高) 7 樣本 V5 BAM 全量重跑**：HCC1395_5kHz (已驗證)、HCC1395_DORADO、HCC1937、HCC1954、H1437、H2009、COLO829；預估 ~10 hr parallel。
- **F4 (高) master × flag 重驗 HPFineNGroups**：對應 S21 機制本質判定，確認是 phasing signature 還是 methylation bimodality。

**F5-F8 (中低優先入 speaker note)**：F5 (中) Manifest 加 haplotag_version 欄位；F6 (中) cnLOH 區獨立評估方案 (R1 應對)；F7 (低) Trio-phased 第二 ground truth (R9 應對)；F8 (高) 100 隨機位點 cross-validate (R3 cherry-pick 風險應對)。

Priority 透明展示：高優先三項 F1 / F3 / F4 是這次 meeting 後優先要做的。

### S24 — Take-home + Next step + Q&A 感謝

結語 v3 強化版 — 倒過來收尾：take-home + next step 承諾 + Q&A 感謝，不重複 S1 TL;DR。

**Take-home**：TO 模式可用、Tag 層問題已解、五大目標解鎖。

三個重點：
1. TO 模式可用 — 透過 V5 修補，TO + PON 條件下 tag 品質接近 paired 水準，read-level concordance +8.3 pp（clean PS blocks，全基因組）。
2. Tag 層問題已解 — V3F 解 priority bug + enum mismatch (Bug 1+2)，V5 補 Layer 1.5 directional reassignment，AMB% 17.5%→8.0%，HP:i:33 reads -54%，介面契約零變動。
3. 五大目標解鎖 — InterSubMod 目標 1/2/4 直接依賴 HP tag 正確，V5 修補是研究可信度前提。

**口徑校準分層 caveat**：V5 不修 self-phasing 本身（V2b 已處理 phase scaffold）；self-phasing 問題鏈是被分層處理的（phase 層由 V2b PON-only，tag 層由 V3F/V5）。

**Next step 3 條承諾**（下次 meeting 教授會看到的進度）：
- F1 Commit V5 working tree（切 2 獨立 commits + tag v5.0）。
- F3 7 樣本擴展 V5 全量重跑。
- F4 HPFineNGroups master × flag 重驗（判定 phasing signature 機制本質）。

以承諾收尾比口號收尾更有力，教授知道下次會看到具體進度。

**v3 結尾強化**：底部 Q&A 歡迎條 + Thank you，雙語並列（CJK + Latin per-char fallback）。

Q&A 隨時歡迎打斷，詳細數據在：
- `v2/source_materials/` 三份報告（01 IGV visual audit、02 V5 audit、03 longphase TO vs V5 技術報告）+ 04 Purity calculator failure root cause（v3 新增）
- `v2/notes/` 6 份 onboarding 文件（00 background、qa_11_questions（含 v3 Q12 Q13）、key_metrics_table、code_references、industry_references、glossary）
- `v2/00_storyboard_v2.md` v3 24-slide ground truth

Thank you。

---

## 附錄 · v3 24-slide 結構速查

| # | 標題 | speaker note 字數 |
|:-:|------|:-:|
| S1 | 17.3:1 → 1:1 + 動機小字 | 996 |
| S2 | Pipeline + 4 階段 | 884 |
| S3 | HP tag 5 整數值 + PS/LOH 兩層 | 964 |
| S4 | 證據三層 (理論+全基因組+SP1) | 820 |
| S5 | Paired 對照確實 ~1:1 | 849 |
| S6 | 拆鎖 phasing vs tag | 795 |
| S7 | LOH 兩層精確化 | 732 |
| S8 | TO 根基 + 4-bucket | 699 |
| S9 | 影響範圍 29/14/42 | 839 |
| S10 | Prerequisite (4 函數+PON+Purity) | 1361 |
| **S11** | **★ 根因樹 (觸發條件主者+三 bug)** | **1565** |
| **S12** | **★ V5 三層 (mid-low purity 防禦層)** | **1775** |
| S13 | 程式碼 diff (somatic-first 紅標) | 1389 |
| **S14** | **★ Sanity 4+5 證據 (Q1 PON 雙路徑)** | **1702** |
| S15 | 量化指標 | 942 |
| S16 | 4-commit timeline | 1187 |
| S17 | 🎯 V5max1 climax + caveat | 1037 |
| S18 | HP 翻轉 bar + IGV | 944 |
| S19 | F1 釐清 + cnLOH 邊界 | 1131 |
| S20 | 業界家族樹 + 2x2 合併 | 1244 |
| S21 | 後續可動 + 已初步現 合併 | 1368 |
| S22 | 五大目標 + 樹 | 721 |
| S23 | P0 行動清單 | 917 |
| S24 | Take-home + Q&A 感謝 | 1242 |

**S11/S12/S14 三張 ★ KEY slide 都 ≥ 500 字**（要求達標）。
**全 24 張 speaker note 都 ≥ 350 字**（要求達標）。

來源依據：v2/source_materials/03 longphase_TO_vs_V5_技術報告.md、04 purity_calculator_failure_root_cause.md、v5_audit_suite/06+07+08+09+18 報告。
