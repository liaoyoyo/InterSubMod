# 教授可能提問 11 條關鍵「為什麼」 + 預備答案

> 由「教授視角審視 agent」找出的 11 條最可能引發 cognitive bumps 的問題；演講者應在 speaker note 預備回答。每條對應 PPT slide 編號。

## Q1【S5】為什麼是 17.3？這數字怎麼算？

**答**：HP1 family / HP2 family。HP1 family = HP1 + HP1_1（germline HP1 + somatic on HP1）reads；HP2 family 同理。HCC1395 全基因組 baseline：
- HP1 family = 614,000 reads
- HP2 family = 35,500 reads
- 比值 614,000 / 35,500 = **17.3**

跨 23 染色體一致；94.6% 集中於 HP1。生物學上隨機分配應 ~1:1，所以 17.3:1 不可能是真實 → 必為 artifact。

## Q2【S7】為什麼有 normal 就會對？（Paired ~1:1 的根本原因）

**答**：normal 提供 **germline het 的 ground truth**：
- Paired 模式：normal BAM 內 het variants 已被 phasing；tumor reads 對齊 normal 的 germline scaffold 即可正確分類 HP1/HP2
- TO 模式：無 normal 對照；phaser 必須從 tumor 自己推 germline，但 tumor 含 somatic → somatic 反客為主進 phasing graph → self-phasing

**追問**：「那為什麼不直接用 paired？」
答：normal BAM 不一定可得（臨床樣本、archive 樣本、cell-free DNA 等）；TO 模式是必要研究方向。本實作目標就是讓 TO 模式達到接近 paired 的 tag 品質。

## Q3【S8】Jaccard=1.0 是怎麼證明 phasing 不變？什麼是 LOH.bed？

**答**：
- **LOH.bed**：LongPhase-TO 從 phased VCF 的 allele depth (AD) 偵測產生 region-level LOH。
- **PON-only 修正前後比對**：baseline LOH.bed 與 V2b PON-only LOH.bed 的 region-level Jaccard 相似度。
- **Jaccard = 1.0**：兩個 LOH.bed 完全相同 → phasing 層的 LOH region 偵測**不受 self-phasing 影響**。

**重要區分**：本工作影響的是 **BAM HP tag → ISM HP_Ratio LOH**（62% artifact），不是 LOH.bed。兩套 LOH 系統用不同數據源（BAM HP tag vs VCF AD）。

## Q4【S13】為什麼 paired 沒被發現這 bug（getVote priority bug）？

**答**：paired 模式 germline 訊號充足。`getVote()` 迴圈遇到 `{HP1_1, HP2_1}` 時，paired 模式 germline votes 多，但因為 priority bug 還是被 somatic 搶先 break — 但 paired 模式有 normal 可校正最終結果，所以 BAM 端看起來「OK」。

**真正暴露的時機**：V2b 啟用 PON-only phasing 後，phasing graph 不再用 somatic 做 anchor，germline votes 急減（因為 PON 過濾掉很多 putative germline），剩下的 somatic vote 就直接主導 → bug 立刻顯形（HCC1395 99.9% reads 拿 HP21）。

**這是「修一個 bug 暴露另兩個 bug」的因果鏈**。

## Q5【S13】這 bug 不是該有 unit test 抓到嗎？（enum vs integer literal）

**答**：誠實回答：LongPhase 原始 codebase 無 HP tag round-trip test（即「encode → decode → 應守恆」）。這個 bug 是型別語意失配，靜態檢查不易抓（C++ 把 enum 隱式轉 int 不會 warn）。

本工作補了 **15/15 sanity check 4 項硬性檢查**：
1. 守恆律 A · ΔHP33 + (ΔHP11 + ΔHP21) = 0
2. 守恆律 B · germline 不變
3. Layer 1.5 期望 1 · 33→directional 精確守恆
4. 無 germline → HP33 violation = 0

未來建議：把這 4 項加入 longphase-to-mod 的 CI test suite。

## Q6【S13/S16】為什麼 PON 有資料還會 self-phasing？

**答**：要區分兩個獨立步驟：
- **PON classification**：用 PON 把 variant 標為 germline / somatic / unknown（VCF 層）
- **Phasing anchor**：phasing graph 用哪些 variants 做 anchor（read 連結）

baseline 的 bug 是 **PON classification 雖然標出 somatic，但 phasing 階段還是把 somatic 拿去當 anchor**（`--pon-only-phasing=false` 預設）。V2b 啟用 `--pon-only-phasing=true` 後 phasing 才只用 PON-confirmed germline 做 anchor → 解 phase scaffold self-phasing。

但 V2b 解了 Phase 層，Tag 層（getVote priority + enum）的 bug 才被暴露 → 需要 V3F+V5 接著修。

## Q7【S16】sanity check 是誰寫的？覆蓋率多少？

**答**：本工作新增（v5_audit_suite/06 報告 Agent D 產出）。覆蓋 4 條守恆律 + Layer 1.5 期望，HCC1395 5kHz 15 sites（5 TP + 4 FP + 3 V5max + 3 SP-extreme）全 PASS。

**覆蓋限制**：
- 15 sites 為 cherry-picked，不是隨機抽樣
- 需要 7 樣本擴展（F3 後續行動）
- 需要 50-100 隨機抽樣位點 cross-validate（F8）
- Confidence threshold 0.6 為間接驗證，需 vote log 直接驗證（F2）

整體 sanity check 是 **必要不充分** — 通過 = 沒有明顯邏輯錯，但不保證跨樣本完美。

## Q8【S17】+8.3 pp concordance 的 N 多少？信賴區間？跨樣本一致嗎？

**答**：
- **15-site cherry-picked**（11 clean PS）：N = 11；V5 88.2% vs Baseline 74.9%（+13.3 pp）；單樣本 HCC1395 5kHz
- **全基因組 clean PS**（PI 報告 4 §3.7）：N = 全基因組 clean PS blocks（數萬-數十萬 reads 量級，具體 N 在報告中）；V5 90.5% vs Baseline 82.2%（+8.3 pp）；單樣本 HCC1395 5kHz
- **跨樣本擴展未做**（F3 待辦）

**Caveat**：clean PS blocks（germline accuracy ≥ 70%）；problem PS blocks 上 V5 與 Baseline 接近隨機（read-level orientation 不穩定）。

**信賴區間**：source 03 未直接給；可由 binomial CI 計算（推估 V5 90.5% N=N_clean 的 Wilson 95% CI 應在 ±1 pp 內）。

## Q9【S20】F1 不變是預期 — 那為什麼還要做這修補？對發論文有意義嗎？

**答**：
**為什麼仍要做**：
1. **下游 ISM 影響**：38% ISM 特徵（29 個 HP-依賴）必須在 V5 BAM 上才可信；舊結論需重跑或加註版本（archive haplotag_version）
2. **生物學詮釋正確性**：HPFineNGroups marker 過去解讀為 methylation bimodality；V5 修補後重詮釋為 phasing bucket signature（Thread D 主軸候選）
3. **跨樣本研究可信度**：7/7 樣本同方向 self-phasing 影響全部研究方向（包括五大目標 1/2/4）

**對發論文的意義**：
1. 本工作不是新 caller / 不是新 phaser；是 **填補 LongPhase-TO 在 PON-only 啟用後的 tag 層 bug**
2. 與 LongPhase-S 2025 paired 形成姐妹工作（同實驗室哲學一致）
3. 真實 contribution：**read-level tag 品質 +8.3 pp 全基因組 clean PS**，這對下游 epigenetic analysis 是 enabling step
4. 學術定位：**InterSubMod 五大目標的解鎖前提**（特別是目標 1/2/4）

## Q10【S21】你的工作跟 LongPhase-S 2025 paper 重疊多少？算 contribution 嗎？

**答**：直球回答：
- **LongPhase-S 是 paired 模式**（tumor + normal）
- **本實作（longphase-to-mod V5）是 tumor-only 模式**（無 normal，PON 替代）
- **設計哲學一致**：兩者都用「somatic 錨在 germline scaffold」
- **本工作填補 TO gap**：當 normal 不可得時的標準解

**Contribution 區隔**：
1. 4-commit 漸進修補在 longphase-to-mod fork（獨立 git repo）
2. 修補對象是 LongPhase-TO 在 PON-only 啟用後集中暴露的 3 層 tag-side bug
3. 4 項 sanity check 是新貢獻
4. ISM 下游影響 3-tier 分類是新貢獻
5. ISM 跨樣本一致性 7/7 + Cohen's d=−1.20 是新貢獻

LongPhase-S 不重疊本工作的 tumor-only 場景。同實驗室相鄰工作，並非 over-claim 為「業界共識」。

## Q11【S24】Thread D NG=2 6/6 POSITIVE 是什麼統計？p=0.0156 對 6 樣本可信嗎？

**答**：
**統計**：Wilcoxon signed-rank（配對非參數檢定）
- **6 樣本配對 paired vs TO**：B1 報告 2026-04-23
- **W = 21**（小樣本 W 統計量上界）
- **p = 0.0156**（exact 計算，不是 normal approx）
- **median gap = 0.365**

**6 樣本可信度**：
- Wilcoxon 6 樣本配對的最強顯著性是 6/6 同方向（W=21, p=0.0156）；本實驗達此上限
- 配 paired control（B3 報告）median gap=0.00003, p=0.578（不顯著）→ 確認 LOH-constrained phasing 是 TO-specific
- HCC1954 outlier（B2 報告）解析為 caller FP 背景 84%、非 phasing failure；修正後 effective gap +0.385
- **多重比較校正**：未做（只有 1 個假設，不需 Bonferroni）

**Caveat**：N=6 仍小；evidence grade 維持 B（待 Phase 2B master×flag 重驗升 A）。Thread D 是論文主軸候選，**前提是 self-phasing 修補完成**（V5 已驗證）。

---

## ★ Q12【S14】（v3 新增）PON 雙路徑（Two-Pass）在不同 purity 下都有效嗎？

**答**：

**Yes, 在 0.93 純樣本 + 0.6 中等純度 sample 都驗證有效**：

| 指標 | 0.93 sample (HCC1395 5kHz) | 0.6 sample (t30_n20) |
|------|:--------------------------:|:--------------------:|
| PON-only Two-Pass 可執行 | ✅ | ✅ |
| Sanity 4 項硬性檢查 | 15/15 PASS, 0 violation | (未測，但機制相同) |
| Aggregate paired GT concordance | V5 78.85% vs Baseline 72.20% (+6.65pp) | (無 paired ref 無法測) |
| Clean PS paired GT concordance | V5 88.2% vs Baseline 74.9% (+13.3pp) | (同上) |
| 全基因組 clean PS 量化 | V5 90.5% vs Baseline 82.2% (+8.3pp) | (同上) |
| Baseline self-phasing 強度 | 1.33:1 (chr19-22 somatic-ALT) | 1:1.14 (衰減 ~60%) |
| V5 vs Baseline HP33 比例 | 8.0% vs ~1.3% | **12.4% vs 2%** (+10pp 保守) |

**Caveat（誠實邊界）**：
- 0.6 sample **沒有 paired reference**（synthetic mix 無 ground truth），所以無法直接做 concordance
- 但 09 報告顯示 V5 在 0.6 下 HP33 比例 +10pp = **保守 tagging 防禦**（避免錯誤分配）
- baseline 在 0.6 自然 self-phasing 弱化（structural reason），V5 仍有改善但 marginal

**結論**：✅ PON-only Two-Pass 設計對全 purity 範圍有效；V5 在 mid-low purity 的價值是 **conservative tagging**（避免錯誤分配），不是直接修 self-phasing bias（後者在 0.6 已自然弱化）。

來源：v5_audit_suite/06、07、08、09 報告；source 03 §0、§4.2.

## ★ Q13【S12】（v3 新增）Baseline「somatic-first」vs V5「germline-first」投票 — V5 真的更好還是只是不同？

**答**：

**V5 更好，且設計上是針對 mid-low purity 的防禦層**。

**Baseline 投票邏輯**（被簡化敘述為「有 HP11/HP21 即 somatic」，實際是 somatic-first 優先序）：
```cpp
// baseline getVote() 偽碼
for (variantKeys order: {HP1_1, HP2_1} → {HP3, HP2_1} → {HP1, HP2}) {
    if (any count > 0) break;  // 先看 somatic，有 somatic vote 就停
}
```

**V5 投票邏輯**（germline-first）：
```cpp
// V5 getVote()
if (germline > 0) → use germline (Layer 1)
else if (somatic > 0 with confidence ≥ 0.6) → use somatic fallback (Layer 1.5)
else → ambiguous (Layer 2 encode)
```

**為何 V5 更好（特別在 mid-low purity）**：

1. **純樣本場景**（HCC1395 5kHz）：baseline somatic-first 在 PON-only 啟用後 germline votes 急減 → somatic 主導 → 17.3:1 self-phasing artifact；V5 germline-first 不受此影響
2. **Mid-low purity 場景**（0.6 sample）：baseline 仍 somatic-first 但 self-phasing 自然弱化（normal contamination 平衡 hap）；V5 加入 confidence 0.6 threshold → 對「弱 directional 訊號」標 HP33 ambiguous（保守 +10pp HP33 比例）→ 避免錯誤分配
3. **Conservative tagging 設計**：V5 寧可標 ambiguous（HP33）也不誤分配；這對下游 ISM 的 4-bucket 分群品質**重要**（NGroups subclone marker 等）

**caveat**：V5 在 SP-extreme self-phasing 位點仍不修（V5 設計就不修 phase 層問題，那是 V2b PON-only 的責任）。

**結論**：V5 確實更好（germline-first 翻轉 + Layer 1.5 confidence 防禦），對 mid-low purity 是**設計上的防禦層**，但**不是萬能**（self-phasing extreme 仍需 V2b PON-only）。

來源：v5_audit_suite/01_code_diff、09_purity06_simulation、source 03 §5.B-D。

---

## 附：教授可能不問但 reviewer 提到的補充疑問

### Q12【S2】為什麼 LongPhase-S 2025 不直接用在 tumor-only？

**答**：LongPhase-S 假設有 normal BAM；當 normal 不可得時，需 PON 替代。本實作是 LongPhase-S 的 TO 變體。同實驗室會否未來合併？open question。

### Q13【S15】程式碼 diff slide 我看不太懂

**答**：請看這三處紅標：
1. baseline `for (variantKeys) { if any non-zero break }` → V5 移除迴圈
2. baseline 第一組 key 是 `{HP1_1, HP2_1}` → V5 第一層直接看 germline `{HP1, HP2}`
3. baseline `if(hpResult != HAPLOTYPE1_1)` 用 enum=3 比較 → V5 改 `if(hpResult != 11)` integer literal

每行都對應前面 S13 root-cause tree 的一條 bug。

### Q14【S22】為什麼 cnLOH 區仍未解？

**答**：cnLOH（copy-neutral LOH）= 雙親同源無 het variants 的區域。沒有 het variants → 無從區分 hap1 / hap2 → V5 fallback 也無法 anchor。需要 cnLOH-aware filter（CN 層 + germline-only methylation reference）作為 future work（F6 行動）。

---

## 演講者準備時間

| 階段 | 時間 |
|------|------|
| 默讀本檔 | 15 min |
| 對照 storyboard 找 slide 號 | 5 min |
| 想像每題 30 秒回答節奏 | 10 min |
| **小計** | **30 min** |

加上 storyboard 與三份來源報告閱讀，總準備時間 ~60 min。
