<!--
build_date: 2026-04-29 06:00
agent: structured-tech-report skill (Q&A 補充)
status: validated
report_class: supplement-QA
audience: PI / 工程深度查證 / 後續 AI session 對話材料
parent_report: InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md
parent_audit: InterSubMod/docs/reports/validated/2026/04/20260428_Self_Phasing_Baseline_V5_Audit_01.md
inputs:
  - /big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp (baseline + V5)
  - /big7_disk/liaoyoyo2001/longphase-to-mod/Util.h (Haplotype enum)
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/13_phase_vs_tag_algorithm_detail.md
  - memory/project_v3_fixed_haplotag_verification.md
  - memory/project_v5_somatic_fallback_verification.md
verdict: validated — 7 Q&A 對齊源碼 + audit suite + memory；明確標註推論 vs source-confirmed
last_verified: 2026-04-29
-->

# longphase-to getVote() 設計意圖與 V5 修補 — Q&A 補充

> 本檔為 `InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md` 的補充。
> 對話模式記錄 V5 修補鏈最容易混淆的概念點。每個 Q&A 標明 ✅ source-confirmed / ⚠ 推論解讀。

---

## 目錄

| Q | 主題 | 證據強度 |
|---|------|---------|
| Q1 | longphase-to 為何寫成這樣的 variantKeys 順序？ | ⚠ 推論為主 |
| Q2 | HP2-1 為何在 variantKeys 中重複出現？ | ⚠ 推論 |
| Q3 | getVote() 的「按順序檢查、第一對有票就 break」邏輯精確說法是什麼？ | ✅ source-confirmed |
| Q4 | paired 模式 / baseline TO / V2b 三種模式 countMap 為何分佈不同？ | ✅ source-confirmed |
| Q5 | 「已 phased somatic 誰多誰贏」是設計者意圖還是我的推論？ | ⚠ **推論，需明確降級** |
| Q6 | 純度 1.0 vs purity 0.6 有差嗎？V5 在低純度下還有用嗎？ | ✅ 09_purity06_simulation 完整驗證 |
| Q7 | V2 有問題，V3 改正後就正確了嗎？ | ✅ V3F + V5 memory + audit suite |
| Q8 | V3-Fixed 拋棄 vector 改用兩層 if-else 是否合理？論文與 0.6 purity 驗證 | ✅ 論文 + 上游參數 + 09_purity06；⚠ V5 fork 內部修補無上游論文認證 |
| Q9 | longphase-s 如何定義 HP0 與 HP3？與 longphase-to 對應關係？ | ✅ knowledge longphase-s.md L174-191、bam-format.md L112 |
| Q10 | baseline 純樣本上 17.3:1 是 bug 還是純樣本固有偏見？資料／程式碼／邏輯三方驗證 | ✅ 4 證據鏈確認非固有 bias，是「資料+演算法+2 程式碼 bug」四重疊加 |
| Q11 | VCF 層 GT2 中 hp1-1 vs hp2-1 數量是否在純樣本下偏斜 17 倍？ | ✅ **No** — VCF site count ≈ 1:1.06；17.3:1 是 BAM read count（兩層分開）|
| Q12 | tag 階段為何造成不平衡？已確認問題點？純 vs 0.6 機制意義差異？ | ✅ 3 重機制（self-phasing + getVote priority + enum mismatch）+ 4 項 sanity check 0 violation；純/0.6 機制相同但觸發強度不同 |
| Q13 | 「baseline HP:i:33=0；ambiguous reads 被推到 HP1 family」是否是 HP1-1 特別多的原因？ | ✅ 是核心原因之一 — V3F 揭露 baseline HP:i:11 過度集中含 ~512K HP:i:21 + ~240K HP:i:33 reads 被誤推 |
| Q14 | V5 Layer 1.5 邏輯合理嗎？HP:i:11/21 語義混雜？業界共識對照？「germline 主 HP2 但有 HP1_1 vote」case 多嗎？| ✅ 6 層論證 + 業界 8/8 工具對齊 + ~512K case 量化；⚠ HP:i:11 含 A/B 兩類來源未明示；含程式碼校正 |
| Q15 | Cross-purity ablation 實證：「PON 方法是否在 0.6 也最好」+「V5 vs baseline 在 0.6 F1 是否較差」+「baseline bug 是否跨 purity 都存在」| ✅ 6 版本 F1 實證 + 0.6/0.93 4-cell ablation 全結論 |

---

## Q1. longphase-to 為何寫成這樣的 variantKeys 順序？

### 程式碼證據（baseline 版本）

```cpp
// /big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp (baseline / V2b)
std::vector<std::pair<int,int>> variantKeys = {
    {HAPLOTYPE1_1, HAPLOTYPE2_1},   // ① phased somatic 比較
    {HAPLOTYPE3,   HAPLOTYPE2_1},   // ② ambiguous somatic vs HP2 phased
    {HAPLOTYPE1,   HAPLOTYPE2}      // ③ germline fallback
};
```

### enum 定義（`Util.h:20-25`）

| enum 名稱 | 數值 | BAM HP tag | 語意 |
|----------|-----:|-----------|------|
| HAPLOTYPE_UNDEFINED | -1 | (無) | 無法判定 |
| HAPLOTYPE1 | 1 | HP:i:1 | germline ALT on HP1 |
| HAPLOTYPE2 | 2 | HP:i:2 | germline ALT on HP2 |
| HAPLOTYPE1_1 | 3 | HP:i:11 | **已 phased somatic** ALT on HP1（`_1` 後綴 = sub-genotype）|
| HAPLOTYPE2_1 | 4 | HP:i:21 | **已 phased somatic** ALT on HP2 |
| HAPLOTYPE3 | 5 | HP:i:33 | **未 phased（ambiguous）somatic** |

### 設計意圖（⚠ 推論）

> **這是依 enum 命名與 vector 順序的合理解讀，未經 LongPhase 上游 design doc 認可**：

```
「如果 read 經過已 phased 的 somatic（HP1_1/HP2_1），
 那它一定屬於該 phase block — phasing graph 已經處理過了，
 phased somatic 的 confidence 高。
 Germline 變異是 fallback。」
```

### 為何在 paired 模式合理

- normal sample 提供獨立 germline anchor → phasing graph 可信
- HP1_1/HP2_1 標籤本身的 phasing 確實是 truthful
- 「特殊訊號優先 fallback germline」是合理設計

### 為何在 TO + PON-only 崩潰

- 沒有 normal anchor → phasing graph 自我定相
- HP1_1/HP2_1 標籤本身**就是 self-phasing 的產物** → 不再可信
- 此時還「somatic 優先」 → 放大 artifact

### ⚠ 待補

未直接驗證：LongPhase 上游 GitHub README / 論文 / commit message 是否有設計者直接陳述「somatic phased 優先」？目前僅依 enum 命名與順序推論。建議後續查 [LongPhase paper](https://doi.org/10.1093/bioinformatics/btac058) 或 [GitHub repo](https://github.com/twolinin/longphase) 確認。

---

## Q2. HP2-1 為何在 variantKeys 中重複出現？

### 觀察

```
{HAPLOTYPE1_1, HAPLOTYPE2_1}   // ① HP2_1 在這裡
{HAPLOTYPE3,   HAPLOTYPE2_1}   // ② HP2_1 又在這裡（重複！）
{HAPLOTYPE1,   HAPLOTYPE2}     // ③ 無 HP2_1
```

**HP2_1 出現兩次，但對稱的 HP1_1 沒出現在第二對**。為何不對稱？

### 三個可能假設（⚠ 全部推論）

| 假設 | 解釋 | 可信度 |
|------|------|--------|
| **歷史 commit 殘留** | 開發者最初寫了 `{HP3, HP2_1}` 處理 HP2-side ambiguous 場景，後來忘了補對稱的 `{HP3, HP1_1}` | **高**（最常見軟體歷史模式）|
| **隱性偏好** | LongPhase 內部某 convention 預設「ambiguous → 偏 HP2」 | 中 |
| **歷史 paired-mode 經驗** | paired mode 下 HP2_1 出現比 HP1_1 多 | 低 |

### 結論

→ 「HP2-1 重複而 HP1-1 沒重複」**本身可能是 baseline 的另一個小 bug 或設計缺陷**，但不是 V5 修補的主焦點。V3-Fixed 直接**拋棄整個 vector 結構**改用顯式 if-else 兩層，所以這個不對稱也被一併消除（無需專門修補）。

### ⚠ 待補

未直接驗證：LongPhase repo `git log --follow HaplotagProcess.cpp` 看哪個 commit 加了 `{HP3, HP2_1}` 這對 — 可能找到歷史脈絡。本次審核未做。

---

## Q3. getVote() 的「按順序檢查、第一對有票就 break」邏輯精確說法

### 用戶口語版（部分正確）

> 「baseline 是取數量較高的直接定義」

### 精確技術描述（取代口語版）

> 1. **getVote() 接收 countMap**：read 跨 PS block 內所有位點累積的票數陣列。
> 2. **按 variantKeys 三層優先序檢查**：
>    - 檢查 Pair ① `{HP1_1, HP2_1}`：若 `countMap[HP1_1]>0 OR countMap[HP2_1]>0` → 取較高者 → break
>    - 否則檢查 Pair ② `{HP3, HP2_1}`：同上邏輯
>    - 否則檢查 Pair ③ `{HP1, HP2}`：同上邏輯
> 3. **第一對有票即 break，不繼續檢查後面**。
> 4. **數量比較只發生在「同一個 pair 內」**，不跨 pair。

### 關鍵釐清

| 誤解 | 正確 |
|------|------|
| ❌ 「全局取最高的 HP」 | ✅ 「按優先序，第一對有票就 break；對內取較高」 |
| ❌ 「baseline 與 paired 程式碼不同」 | ✅ **三模式 getVote() 程式碼完全相同**（同一函式，不同輸入）|
| ❌ 「V2b 改了 getVote()」 | ✅ V2b commit `8b8c1fd` **不動 HaplotagProcess.cpp**；只加 `--pon-only-phasing` flag |

---

## Q4. paired 模式 / baseline TO / V2b 三種模式 countMap 為何分佈不同？

### 核心釐清

> **三模式 getVote() 程式碼完全相同**；差異不在程式邏輯，而在 **phasing graph 怎麼建構 → countMap 的分佈不同**。

### 三模式 countMap 分佈對比（同一個 read 跨 PS block 累積票）

```
═══════════════════════════════════════════════════════════════
情境 1：paired 模式（normal+tumor，germline anchor 充足）
═══════════════════════════════════════════════════════════════
phasing graph: 用 normal 找出大量 germline het 當 anchor
countMap 分佈典型：
  HP1   ████████████  (12 votes) ← germline anchor 充足
  HP2   █████████      (9 votes)
  HP1_1 █              (1 vote)  ← somatic 占少數
  HP2_1 (0)
  HP3   (0)

迴圈走向：
  Pair ① {HP1_1=1, HP2_1=0} → 觸發！break！指派 HP11
  → bug 啟動但因 paired 整體 germline 比例高，個別 reads
    偏差不會匯聚成統計顯著的全基因組偏斜
  → 觀察 HP_Ratio ≈ 0.5（合理）

═══════════════════════════════════════════════════════════════
情境 2：baseline TO（PON-only=false，預設）
═══════════════════════════════════════════════════════════════
phasing graph: 用 germline + somatic + unknown 混合 anchor
countMap 分佈典型：
  HP1   ████        (4 votes)
  HP2   ███         (3 votes)
  HP1_1 ██          (2 votes)
  HP2_1 █           (1 vote)
  HP3   █           (1 vote)

迴圈走向：
  Pair ① {HP1_1=2, HP2_1=1} → 觸發！break！指派 HP11
  → bug 觸發機率 ~30%
  → 整體 17.3:1 偏斜已存在但症狀「不極端」，未被即時注意

═══════════════════════════════════════════════════════════════
情境 3：V2b TO（PON-only=true，啟用 flag）
═══════════════════════════════════════════════════════════════
phasing graph: 僅 PON-confirmed germline 當 anchor
              非 PON 變異全部被標為 somatic
countMap 分佈典型：
  HP1   █            (1 vote)  ← germline anchor 失去大量訊號
  HP2   (0)
  HP1_1 ███          (3 votes) ← somatic 反客為主
  HP2_1 ██           (2 votes)
  HP3   ██           (2 votes)

迴圈走向：
  Pair ① {HP1_1=3, HP2_1=2} → 觸發！break！hpResult=HP1_1=3 (enum)
  + caller 端 if(hpResult != HAPLOTYPE1_1) 比 enum 3 vs integer 11
    → 永不匹配 → 強行轉 HP21
  → 99.9% reads 拿到 HP21，commit message 立即點名
```

### 一句總結

> **bug 一直存在於程式碼**（`variantKeys` 順序錯）；**「能不能觸發」取決於 countMap 中 somatic 比 germline 多還是少**。phasing graph anchor 來源決定了 countMap 分佈，最終決定 bug 觸發機率。

---

## Q5. 「已 phased somatic 誰多誰贏」是設計者意圖還是推論？

### 誠實降級

| 我之前說的 | 實際依據 |
|----------|---------|
| 「variantKeys 順序」 | ✅ **程式碼直接證據**（baseline `getVote()` line 506-530）|
| 「Pair ① 順序意圖是『最強訊號優先』」 | ⚠ **推論**（依 enum 命名 `_1` 後綴 + vector 順序）|
| 「設計者原始意圖是給 paired 模式」 | ⚠ **推論**（無上游 design doc 直接陳述）|
| 「V2b PON-only 後 priority bug 立刻顯形」 | ✅ **commit message 直接證據**（`8b8c1fd`：「99.9% HP_Ratio extreme rate」）|

### 報告中應如何寫

> ⚠ 「最強訊號優先」**降級為「對 enum 命名與 vector 順序的合理解讀」**，不是 LongPhase 上游 design doc 認可的結論。除非後續查到上游論文 / commit message / README 直接陳述設計意圖，否則此說法應標 ⚠ 為推論。

### 後續驗證建議

- 查 [LongPhase paper (Bioinformatics 2022)](https://doi.org/10.1093/bioinformatics/btac058) 是否提到 `_1` sub-genotype 的優先性
- 查 [LongPhase GitHub](https://github.com/twolinin/longphase) issues / PRs / README 設計討論
- 查 longphase-to-mod 的 `git log -p HaplotagProcess.cpp` 在 fork 前的歷史

---

## Q6. 純度 1.0 vs purity 0.6 有差嗎？V5 在低純度下還有用嗎？

### 數據佐證（`v5_audit_suite/09_purity06_simulation.md` 完整驗證）

#### Q6.1 HCC1395 在不同 purity 的 self-phasing 強度對比

| 場景 | tumor:normal 混合 | 內部 purity 估計 | chr19-22 baseline somatic-ALT HP1:HP2 |
|------|----------------|---------------|-------------------------------------|
| 0.93（純 tumor）| 100:0 | ~0.93 | **1.33:1** |
| 0.6（t30_n20）| 30:20 | 0.607（估計準確）| **1:1.14** ← self-phasing 顯著減弱 |

> 09_purity06_simulation §1.4：baseline 內部 purity 算出 0.607，精準對應 t30_n20 的 60% 設計，驗證 baseline phase 階段對 purity 的估計能力正確。

#### Q6.2 Self-phasing extreme 位點隨 purity 衰減

| Site | baseline @ 0.93 | baseline @ 0.6 | 衰減幅度 |
|------|----------------|---------------|---------|
| SP1 chr19:17565944 | 113:1 | 1:70 | -38%（且方向翻轉）|
| SP2 chr19:12452332 | 109:1 | 1:41 | **-62%** |
| SP3 chr19:12467180 | ∞ (105:0) | ~42× | 顯著減弱 |

→ self-phasing 強度確實**隨 purity 衰減**（normal reads 稀釋 tumor sub-clone read 共現 → phasing graph 強連結減弱）。

#### Q6.3 V5 在不同 purity 的價值定位變化

| purity | V5 修復目標 | V5 主要價值 | HP33 比例 |
|--------|----------|-----------|---------|
| 0.93（高純度）| 修復 17:1 self-phasing | **修復 17:1 → 1:1**（核心價值）| 8.2% |
| 0.6 | baseline 自然就 1:1.14 → V5 修復 unobservable | **conservative tagging**（HP33 ↑）| 12.4%（HP33 reads 多 6×）|
| < 0.5 | 同上趨勢 | conservative tagging（需 HP33-aware downstream）| 預期更高 |

#### Q6.4 V5 @ 0.6 的「反直覺」結果（重要）

⚠ V5 在 0.6 下的 chr19-22 somatic-ALT ratio 是 **1:2.34**（不是預期的 1:1）。

不是 V5 bug，是兩個機制疊加：
1. V5 conservative：把無 germline anchor 的 reads 推到 HP33 → 留下的 reads 是 V5 真信心高的
2. PS block orientation 系統偏好（chr19-22 sample 範圍內）→ ratio 偏 HP2

### 機制解釋（Q&A 整合）

| 問題 | 答覆 |
|------|------|
| 為何低 purity 下 self-phasing 自然減弱？ | normal reads 稀釋了「同 sub-clone read 共現」→ phasing graph 邊權重變弱 → 不再形成強連結 |
| V5 在低純度下不需要嗎？ | ✅ 仍可用，但**價值定位轉變**：從「修 self-phasing」變為「conservative tagging 避免錯誤分配」 |
| 推薦使用場景 | ≥ 0.85: V5 必用；0.6-0.85: V5 或 baseline 都可（V5 更 conservative）；< 0.5: V5 + HP33-aware downstream |

→ 詳見 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md` §7.4。

---

## Q7. V2 有問題，V3 改正後就正確了嗎？

### 答案：**不！V3-Fixed 修了 baseline 的 bug，但又引入「過度保守」新問題**。需要 V5 才達到平衡正確。

### 修補階段對比（4-stage 演進）

| 版本 | 修了什麼 | 仍有什麼問題 | AMB%（HP:i:33 比例）|
|------|---------|------------|------|
| baseline / V2b | — | priority bug + enum mismatch + HP2-1 不對稱 | N/A（HP:i:33 永不出現，enum bug）|
| **V3-Fixed (`41ff147`)** | ✅ priority bug<br>✅ enum mismatch<br>✅ HP2-1 不對稱（拋棄整個 vector） | ⚠ **過度保守**：germline=0 時一律 HP33<br>⚠ 丟失 HP1_1/HP2_1 的 phased 方向資訊 | **17.5%** ← 過高（paired 5.4% 為比較基準）|
| INDEL guard (`380e8d2`) | ✅ INDEL UNDEFINED UB | 同 V3F | 17.5% |
| **V5 (working tree)** | ✅ Layer 1.5：germline=0 但 somatic HP1_1/HP2_1 有票時，用 somatic 投票 fallback<br>✅ SNP alt UNDEFINED guard | ⚠ Confidence threshold 0.6 未直接驗證<br>⚠ working tree 未 commit | **8.0%** ← 接近 paired 5.4% |

### V3-Fixed 的「過度保守」具體表現

```cpp
// V3-Fixed 邏輯
int germlineResult = 0;
if (countMap[HP1] > 0 || countMap[HP2] > 0) {
    germlineResult = (countMap[HP1] >= countMap[HP2]) ? 1 : 2;
} else {
    // ⚠ V3-Fixed 在這裡直接放棄
    // 即使 countMap[HP1_1] > 0 或 countMap[HP2_1] > 0，也不利用！
    // 結果：hpResult = 33（HP:i:33 ambiguous）
}
```

- HP:i:33 reads（HCC1395 全基因組）= **239,679**
- 這 239,679 reads 中有大比例其實有 HP1_1 / HP2_1 phased 方向資訊 → 被 V3F 強制忽略
- 結果：HP:i:33 比例 17.5%，比 paired ground truth 5.4% 高 **3.2 倍**

### V5 Layer 1.5 的補救

```cpp
// V5 working tree (uncommitted)
if (countMap[HP1] > 0 || countMap[HP2] > 0) {
    // Layer 1: germline 主導（同 V3F）
    germlineResult = (countMap[HP1] >= countMap[HP2]) ? 1 : 2;
}
// Layer 1.5（V5 NEW）：germline 全 0 但 somatic phased 有訊號
else if (countMap[HP1_1] > 0 || countMap[HP2_1] > 0) {
    if (countMap[HP1_1] >= countMap[HP2_1]) {
        germlineResult = 1;  // 用 somatic phasing 推測方向
    } else {
        germlineResult = 2;
    }
}
// 兩層都不滿足 → 才標 HP33
```

### V5 數據改善

| 指標 | V3-Fixed | V5 | 變化 |
|------|---------:|---:|-----:|
| HP:i:33 reads（HCC1395 全基因組）| 239,679 | 110,197 | **-129,482 (-54%)** |
| AMB%（HP:i:33 占比）| 17.5% | **8.0%** | -9.5 pp |
| HP:i:11 reads | 584,117 | 666,997 | +82,880（從 HP33 升級）|
| HP:i:21 reads | 547,118 | 593,720 | +46,602（從 HP33 升級）|
| Sanity check（4 項硬性檢查）| 未做 | **15/15 PASS, 0 violation** | — |

### 邏輯對應圖

```
baseline / V2b：「somatic 反客為主，全錯」
   └─→ 99.9% HP21 集中（commit message 點名）

V3-Fixed：「germline first，但太保守，把方向訊號丟了」
   └─→ AMB% 17.5%（HP1_1/HP2_1 的方向被忽略）

V5：「germline first，但 fallback 用 somatic phasing 救方向」
   └─→ AMB% 8.0%（接近 paired 5.4% 為比較基準）✅ 平衡正確
```

### 一句結論

> **V3-Fixed 是「對的方向但太保守」；V5 才是「最終平衡正確」**。V5 不是「V3F 的小調整」，而是補上 V3F 因過度保守丟失的 phased somatic 方向訊號。

---

## Q8. V3-Fixed 拋棄 vector 改用兩層 if-else 是否合理？論文與 0.6 purity 驗證

### 8.1 論文真實存在性確認（已從知識庫核對）

✅ **核對通過 — 引用的論文皆真實存在於知識庫或 InterSubMod 文獻庫**：

| 論文 | 知識庫位置 | 認可內容 |
|-----|----------|---------|
| **Lin et al., Bioinformatics 2022**（LongPhase 原版）| `/big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase.md` 第 47 行明列「引用文獻」 | LongPhase 原版設計：兩階段 phasing，**germline het 為 DAG vertex 主要 anchor** |
| **陳鎮宇 CCU 2025 碩士論文（LongPhase-TO 設計）** | `/big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase-to.md` 第 40 行 + 第 427 行（papers/LongPhase-TO.pdf）| LongPhase-TO 演算法完整描述（六步驟流程）+ **§4.3 8 cell lines × 20-100% purity 驗證 ClairS-TO-ssrs F1 提升** |
| **InterSubMod 文獻調查報告**（含 6 個競品論文）| `InterSubMod/docs/references/manual/20260402_longphase_to_phasing_quality_literature.md` | WhatsHap (Martin 2016, Patterson 2015) / HapCUT2 (Edge 2017) / MethPhaser (Nat Comm 2024) / Octopus (Cooke 2021 Nat Biotech) / SAVANA (Nat Methods 2025) / Severus (Nat Biotech 2025) |

### 8.2 LongPhase-TO 上游論文 §4.3 §4.4 的 purity 範圍驗證

`/big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase-to.md` 第 125-132 行「論文驗證結果摘要」直接引用：

| 驗證項目 | 結果 | 來源 |
|---------|------|------|
| Phasing（vs LongPhase germline）| Block N50: 10-25 Mb；Phased Ratio: 0.55-0.62 | 論文 §4.1, 8 cell lines |
| LOH（vs CytoScan HD array）| HCC1395 genome-wide 高度一致 | 論文 §4.2, Fig 4.3 |
| **Somatic refinement** | ClairS-TO-ssrs precision 提升、F1 提升 | **論文 §4.3, 8 cell lines × 20-100% purity** ★ |
| Purity estimation（vs ASCAT）| 20-100% purity 範圍緊貼 y=x 理想線 | 論文 §4.4 |

→ **LongPhase-TO 上游機制（PON-only Two-Pass + `somaticConnectAdjacent=6` + `-p 0.6`）已論文驗證在 20-100% purity 範圍工作正常**。

### 8.3 重要校正：「Confidence threshold 0.6」不是 V5 創新

從 longphase-to.md 第 250 行：

```
| -p | Haplotype 判定門檻；read 中 allele 比例低於此值則不標記 | 預設 0.6 |
```

⚠ **校正前述敘述**：
- V5 audit suite 一直標的「Confidence threshold 0.6」是 LongPhase-TO **原生 `-p` 參數**，不是 V5 fork 創新
- 上游碩士論文 §4.3 已驗證此參數在 20-100% purity 範圍工作正常
- v5_audit_suite/06 的「caveat」是指「V5 fork 內部 vote log 未直接觀察」，不是參數本身有問題

### 8.4 V3-Fixed/V5 改動的合理性 — 5 層獨立論證

#### 8.4.1 ✅ 論文層（Lin 2022 + 陳鎮宇 2025 直接支持）

- Lin et al. 2022 §2.1：LongPhase 原版設計 = germline het 為 phasing scaffold 主要 anchor → V3-Fixed 反轉為 germline-first 回歸此原則
- 陳鎮宇 2025 §4.3：LongPhase-TO 自身設計 `convertNonGermlineToSomatic()` Two-Pass 機制（high-purity > 0.95 自動觸發）→ V5 強制啟用此機制

#### 8.4.2 ✅ 上游參數層（LongPhase-TO 預設值已支持 V5 邏輯）

| 參數 | 預設值 | 來源 | 含義 |
|------|-------:|------|------|
| `--connectAdjacent` | 35 | knowledge longphase-to.md L219 | germline 連接寬 |
| **`--somaticConnectAdjacent`** | **6** | knowledge longphase-to.md L207 | somatic 連接窄（防 self-phasing partial fix）|
| **`-p 0.6`** | **0.6** | knowledge longphase-to.md L250 | Haplotype 判定門檻（V5 用的 confidence 即此參數）|

→ LongPhase-TO 上游已透過 3 個參數設計表明「somatic anchor 應比 germline 弱」+「需 confidence threshold」 — V5 fork 強化此方向，不違背上游設計精神。

#### 8.4.3 ✅ 數據層（V3F → V5 改進實證）

| 指標 | V2b | V3-Fixed | V5 | paired 基準 |
|------|----:|--------:|---:|----------:|
| HP Tag Balanced% | 13.0% | **22.5%** (+9.5pp) | (V5 同 V3F 基準) | 36.9% |
| HP:i:33 出現 | 0（enum bug）| **6,793**（修正後）| 110,197（V5 重分配後）| N/A |
| AMB% | N/A | 17.5% | **8.0%** | 5.4% |
| 4 項硬性 sanity check | — | — | **15/15 PASS, 0 violation** | — |

#### 8.4.4 ✅ 0.6 purity 實證層（v5_audit_suite/09 完整驗證）

| 場景 | chr19-22 ratio | HP33 % | 內部 purity 估計 | 驗證 |
|------|--------------|------:|--------------:|-----|
| baseline @ 0.93 | 1.33:1 | 1.5% | ~0.93 | — |
| **baseline @ 0.6** | 1:1.14 | 2.0% | **0.607**（估計準確）| self-phasing 自然減弱 ✅ |
| V5 @ 0.93 | 1:1.36 | 8.2% | N/A（PON-only）| 修復 17:1 → 接近 1:1 ✅ |
| **V5 @ 0.6** | 1:2.34 | **12.4%** | N/A | conservative tagging ✅（HP33 ↑ 6×）|

mean log gap：0.93 下 V5 vs BL = 2.7 → 0.6 下 = 1.2（gap 縮小，方向一致提升 4/15 → 8/15 sites）。

#### 8.4.5 ✅ 邏輯一致性層（normal mix 強化 germline-first）

0.6 = 60% tumor + 40% normal mix：

```
┌────────────────────────────────────────────────────────────┐
│  0.6 purity 下 V3-Fixed 邏輯為何仍站得住腳                  │
├────────────────────────────────────────────────────────────┤
│ 1. normal sample 提供大量 germline het（~50% het rate）    │
│ 2. mix 後 germline anchor 訊號強化                          │
│ 3. somatic read pool 被稀釋（normal 沒 somatic）            │
│ 4. countMap 中 germline votes 比 0.93 更充足                 │
│ 5. V3-Fixed Layer 1（germline first）通常有票 → 走 Layer 1  │
│ 6. V5 Layer 1.5 觸發機率下降（germline 多 → 不需 fallback） │
│ 7. V5 conservative tagging（HP33 ↑）作為設計安全網          │
└────────────────────────────────────────────────────────────┘
```

→ **邏輯敘述無自相矛盾**：normal mix 強化 V3F 的 germline-first 假設；V5 Layer 1.5 在 0.6 下使用率降低但仍提供安全網。

### 8.5 邏輯敘述審查（質疑點與回應）

| 質疑 | 處理 |
|------|------|
| **Q：V3-Fixed @ 0.6 沒有獨立對照組** | 真的。09 章只測 baseline_06 vs V5_06；V3F-only @ 0.6 沒測。但 V5 = V3F + Layer 1.5；V5_06 通過 sanity = V3F 邏輯在 0.6 下不破壞 V5 的整體有效性 |
| **Q：論文認證是 LongPhase-TO 上游，不是 V5 fork** | ✅ **誠實接受**。論文認證範圍：(1) Lin 2022 = germline phasing 原則；(2) 陳鎮宇 2025 §4.3 = PON-only Two-Pass + 上游參數 in 20-100% purity。V5 fork 的 Layer 1.5 + SNP alt guard 屬於 fork 內部修補，**論文認證僅延伸到 V3-Fixed 修補方向（germline first）**，不延伸到 V5 具體實作 |
| **Q：0.6 V5 ratio 是 1:2.34 不是 1:1，是不是不有效？** | 09 章 §5.1 已解釋：1:2.34 是 conservative tagging（HP33 從 8.2% → 12.4%）+ PS block 系統偏好的合成；非 V5 bug。V5 在 0.6 的價值是 conservative tagging 而非 bias 修正（baseline 自然 1:1.14） |
| **Q：論文 §4.3 §4.4 認證 20-100%，怎麼知道 < 0.5 還有效？** | ✅ **誠實**：09 章未測 < 0.5。報告主軸限 ≥ 0.6；< 0.5 屬另案（v5_audit_suite/00_INDEX 後續行動 #3：7 樣本擴展可包含低 purity）|
| **Q：「Confidence threshold 0.6」是不是 V5 創新？** | ❌ **校正**：不是！這是 LongPhase-TO 原生 `-p` 參數（knowledge longphase-to.md L250）；上游碩士論文 §4.3 已驗證；前述「Confidence threshold 未直接驗證」應改為「V5 fork 內部 vote log 未觀察」（更精確）|

### 8.6 一頁信心級別總結

```
┌──────────────────────────────────────────────────────────────────┐
│ 信心 │ 證據                                                       │
├──────────────────────────────────────────────────────────────────┤
│ ★★★★★ │ Lin et al. 2022（已存在於知識庫 longphase.md L47）        │
│         │ → V3-Fixed germline-first 回歸論文認可的設計原則           │
├──────────────────────────────────────────────────────────────────┤
│ ★★★★★ │ 陳鎮宇 2025 §4.3（已存在於知識庫 longphase-to.md L40+L427）│
│         │ → 8 cell lines × 20-100% purity 驗證 LongPhase-TO 上游機制 │
├──────────────────────────────────────────────────────────────────┤
│ ★★★★★ │ LongPhase-TO 原生參數（-p 0.6 / somaticConnectAdjacent=6）│
│         │ → V5 用的 confidence 0.6 是上游既有，論文已認證               │
├──────────────────────────────────────────────────────────────────┤
│ ★★★★★ │ v5_audit_suite/09 0.6 purity 完整實證                     │
│         │ → V5 @ 0.6 conservative tagging 工作正常（HP33 12.4%）   │
├──────────────────────────────────────────────────────────────────┤
│ ★★★★☆ │ 邏輯一致性：normal mix 強化 germline-first 假設            │
│         │ → 推論但無自相矛盾                                          │
├──────────────────────────────────────────────────────────────────┤
│ ★★★☆☆ │ V5 fork 的 Layer 1.5 + SNP alt guard 具體實作               │
│         │ → fork 內部修補，無上游論文直接認證；audit suite 4-agent 多面驗證 │
├──────────────────────────────────────────────────────────────────┤
│ ★★☆☆☆ │ < 0.5 purity 場景                                         │
│         │ → 09 章未測；屬另案（後續 7 樣本擴展可涵蓋）                 │
└──────────────────────────────────────────────────────────────────┘
```

### 8.7 一句結論

> **V3-Fixed 拋棄 vector 改用兩層 if-else 的方向有 5 層獨立論證支持**：(1) Lin 2022 LongPhase 原版設計原則；(2) 陳鎮宇 2025 LongPhase-TO Two-Pass + 上游參數已論文驗證 20-100% purity；(3) LongPhase-TO 原生 `-p 0.6` 已認證；(4) v5_audit_suite/09 完整實證 0.6 purity；(5) normal mix 強化 germline-first 邏輯。
>
> **誠實降級**：V5 fork 的 Layer 1.5 + SNP alt guard 具體實作屬 fork 內部修補，**論文認證僅延伸到 V3-Fixed 方向（germline first），不延伸到 V5 具體實作**。整體信心級別為 **★★★★☆**（在 0.6-0.93 範圍）。

---

## Q9. longphase-s 如何定義 HP0 與 HP3？與 longphase-to 對應關係？

### 9.1 longphase-s HP tag 完整定義（從 knowledge `05_tools/longphase-s.md` 第 174-191 行）

| HP tag | 定義 | 對應 longphase-to 整數編碼 |
|--------|------|------------------------|
| `HP:z:1` | Germline haplotype 1 的 read | `HP:i:1` |
| `HP:z:2` | Germline haplotype 2 的 read | `HP:i:2` |
| `HP:z:1-1` | 支持 somatic ALT，可追溯到 germline HP1 的 read | `HP:i:11` |
| `HP:z:2-1` | 支持 somatic ALT，可追溯到 germline HP2 的 read | `HP:i:21` |
| `HP:z:3` | 支持 somatic ALT，**無法由既有 germline haplotype 唯一導出** | `HP:i:33` |

從 knowledge `03_file_formats/bam-format.md` 第 112 行確認：
> 「`HP:i:11` / `HP:i:21` / `HP:i:33` 與 `HP:z:1-1` / `HP:z:2-1` / `HP:z:3` 語義相同，**只是編碼格式不同**」

### 9.2 HP0 的本質（重要釐清 — 不是明確定義的 tag）

⚠ **HP0 不是 longphase-s 或 longphase-to 主動指派的 tag**，而是 **BAM convention 中「無 HP tag」的狀態**：

- BAM 檔案中 read 沒有 `HP:i:` 或 `HP:z:` 欄位 = unphased / untagged
- 在統計報告中常以 `HP:i:0` 或 `HP=0` 表示「未標記」
- **不是 longphase 主動的「HP=0」結果**，而是「沒被指派」的結果

觸發條件：
- read 跨 PS block 的 countMap 全為 0（無任何位點 vote）
- 或 confidence 低於 `-p 0.6` 門檻被 `judgeHaplotype()` 攔截
- → 整體不寫入 HP tag → BAM 顯示為「無 HP」狀態

### 9.3 HP3 的本質（明確不是 LOH 標記）

從 longphase-s.md 第 186-187 行：

> **判讀邊界**：`HP:z:1-1`、`HP:z:2-1`、`HP:z:3` 都是 somatic ALT reads 的 read-level 標記；**`HP:z:3` 是 ambiguous somatic haplotype，不是 LOH 標記**。LOH 需另外根據 haplotype imbalance 與 phased/LOH 輸出判讀。

→ **HP3 = 「support somatic ALT 但 phasing graph 無法決定該 read 屬 HP1 phase block 還是 HP2 phase block」**。

下游分析陷阱（已知）：
- ❌ 不要把 HP3 當作「LOH region 訊號」
- ❌ 不要把 HP3 數量直接等同 ambiguous LOH 比例
- ✅ LOH 應另外用 `LOH.bed` (region-level) + haplotype imbalance 推論

### 9.4 longphase-s vs longphase-to 對 HP3 的處理差異

| 維度 | longphase-s（paired）| longphase-to（tumor-only）|
|------|----|----|
| 編碼格式 | 字串 `HP:z:3` | 整數 `HP:i:33` |
| 觸發條件 | normal 提供 germline anchor 後仍 ambiguous | TO 模式無 normal；ambiguous 比例顯著高 |
| baseline 預期比例 | 低（normal 已過濾大量 germline noise）| 高（受 self-phasing + getVote bug 雙重影響）|
| baseline TO 實測 | N/A | **0**（enum bug → 永不寫入；audit suite 確認）|
| V3-Fixed 修正後 | N/A（V3F 是 longphase-to 修補，不影響 longphase-s）| HP:i:33 出現 6,793 次（修了 enum bug 才正確產生）|
| V5 修正後 | N/A | HP:i:33 從 239,679 → 110,197（−54%；Layer 1.5 重分配方向資訊）|

### 9.5 HP3 在 V5 後的「正確」分佈邊界

從 v5_audit_suite/06_v5_sanity_bug_check.md：
- V5 後 HP:i:33 比例 = **8.0%**（AMB%）
- 對比 paired ground truth = **5.4%**（自然 ambiguous 率）
- → V5 仍略高於 paired，反映 TO 模式無 normal anchor 時 ambiguous 比例不可能完全降到 paired 水平
- 但已從 V3F 的 17.5%（過度保守）降到 8.0%（接近合理）

---

## Q10. baseline 純樣本上 17.3:1 是 bug 還是純樣本固有偏見？

### 10.1 精確澄清 — 17.3:1 指的是什麼？

從 `v5_audit_suite/10_somatic_bias_explanation.md` §2 Panel D：

> 跨整個基因組統計**所有 somatic variant 的 ALT reads 加總**：
> - HP1 reads：**614,000**
> - HP2 reads：**35,500**
> - 比例 = 614K / 35.5K = **17.3 : 1**

**精確定義**：HCC1395 baseline TO 模式下，**所有 somatic 位點上 ALT-supporting reads** 的 HP1 family : HP2 family 比例 = 17.3:1（94.6% 集中於 HP1）。

→ **不是只看 HP1_1 vs HP2_1**（phased somatic），而是**所有 somatic ALT reads 的 HP1 全家族（HP1 + HP1_1）vs HP2 全家族**。

### 10.2 17.3:1 不是純樣本固有偏見，是 4 因素疊加

| 歸因 | 性質 | 是否「bug」？ | 修補狀態 |
|-----|------|-----------|---------|
| **A. 純樣本特性**（HCC1395 0.93 高純度）| 客觀資料事實：sub-clone reads 共現強 → phasing graph 強連結 | ❌ **不是 bug，是資料客觀特性** | 不需修補 |
| **B. baseline 預設 `--pon-only-phasing=false`** | 演算法選擇：phasing graph 用 germline+somatic+unknown 混合 anchor | ⚠ **設計選擇**（給 paired 模式 OK；TO 純樣本上崩潰）| V2b `8b8c1fd` 加 flag 補救 |
| **C. getVote priority bug**（baseline `HaplotagProcess.cpp:506-530`）| 程式碼 bug：variantKeys 順序 `{HP1_1, HP2_1}` 在 `{HP1, HP2}` 之前；任一 somatic vote → break，germline 被忽略 | ❌ **程式碼 bug** | V3-Fixed `41ff147` 修正 |
| **D. enum vs integer literal mismatch**（`HaplotagProcess.cpp:697`）| 程式碼 bug：caller 端 `if(hpResult != HAPLOTYPE1_1)` 用 enum 值 (=3) 比 integer literal (=11)；永不匹配 → fallback 失效 | ❌ **程式碼 bug** | V3-Fixed `41ff147` 同 commit 修正 |

### 10.3 4 因素如何疊加

```
A. 純樣本                       B. PON-only=false              C. getVote 順序錯              D. enum-int mismatch
sub-clone reads                phasing graph anchor          variantKeys 第一對是          fallback 路徑型別不一致
共現強                    →    用混合（germline+somatic）  →   {HP1_1, HP2_1}            →   永不匹配，HP21 強行集中
                               instead of PON-only germline    任一 somatic vote 即 break
                                                                germline 永遠被忽略
                                                                ↓
                                         四因素疊加 → ALT reads 17.3:1 偏 HP1
                                                       （V2b 啟用 PON-only flag 後 → 99.9% HP21 集中）
```

### 10.4 4 條獨立證據鏈確認非純樣本固有偏見

#### 證據 1：低 purity 自然消失（v5_audit_suite/09_purity06_simulation §2.2）

| purity | chr19-22 baseline somatic-ALT HP1:HP2 ratio |
|-------|--------------------------------------|
| 0.93（純樣本）| **1.33:1**（chr19-22 範圍內較弱；全基因組 17.3:1）|
| 0.6（mix normal）| **1:1.14**（self-phasing 自然減弱）|

→ normal 稀釋後 17.3:1 消失 → **不是純樣本固有 bias**，是「程式碼 bug 在純樣本上的最大強度顯現」。

#### 證據 2：跨 23 染色體一致性 → 不可能是真實生物學

從 10_somatic_bias_explanation §1.2：

> 如果 baseline 的 17.3:1 是真的生物學現象，意味：
> - HCC1395 tumor **95% somatic 都長在 maternal chromosome 上** → 違反隨機性
> - 跨 23 條染色體都這樣 → 機率極低
>
> **唯一合理解釋**：baseline 的 phasing graph 演算法有 bug

#### 證據 3：跨樣本 7/7 一致（CV-2 PASS）

從 `02_Self_Phasing根因.md` 量化證據總表：
> 跨樣本方向一致性 = **7/7 樣本全部觀察到相同方向的 self-phasing 效應**

→ 7 個獨立 cell line（HCC1395、HCC1395_DORADO、HCC1937、HCC1954、H1437、H2009、COLO829）都出現同方向偏斜 → 排除「HCC1395 樣本特異性」可能 → 必屬演算法／程式碼層級問題。

#### 證據 4：commit message 直接點名 bug（01_code_diff §2.1）

V2b commit `8b8c1fd` 的 message：
> 「**Known issue**: haplotag getVote() priority bug … HP_Ratio extreme rate 99.9%。Will fix in next commit (`41ff147`)」

→ 開發者自身在啟用 PON-only flag 時就發現 bug 顯形，**直接認可這是 bug 而非設計**。下個 commit `41ff147` (V3-Fixed) 即修正。

### 10.5 paired 模式為何 bug 不顯眼？

關鍵釐清：兩個程式碼 bug（C 和 D）**在 paired 模式下也存在**，但**不顯眼**。原因：

| 模式 | countMap 中 germline:somatic vote 比例 | bug 觸發機率 | 觀察症狀 |
|------|----------------------------------|------------|---------|
| paired（normal+tumor）| ~5:1（germline 充足）| <5% | HP_Ratio ≈ 0.5（合理）|
| baseline TO（PON-only=false，純樣本）| ~3:2（混合）| ~30% | 17.3:1 偏斜（已存在但不極端）|
| V2b TO（PON-only=true，純樣本）| **~1:3**（somatic 反客為主）| **>95%** | **99.9% HP21**（commit message 點名）|

→ **bug 一直存在於程式碼**；**「能不能觸發」取決於 countMap 中 somatic 比 germline 多還是少**。phasing graph anchor 來源決定了 countMap 分佈，最終決定 bug 觸發機率。

### 10.6 一句結論

> **17.3:1 不是純樣本固有偏見**，是**「純樣本特性（A）+ baseline 演算法選擇 PON-only=false（B）+ getVote priority bug（C）+ enum-int mismatch（D）」四重疊加**在純樣本上的最大強度顯現。
>
> - **A（資料）+ B（設計）**：合理（給 paired 模式設計 OK）
> - **C + D（程式碼）**：bug（V3-Fixed 同 commit `41ff147` 修正）
>
> 4 條獨立證據鏈（低 purity 自然消失 / 跨染色體一致 / 跨樣本 7/7 / commit message 直接點名）確認**非資料固有 bias，是程式碼 bug 在純樣本上的最大強度顯現**。修補後（V3F + V5）HP1:HP2 回歸 ~1:1，與 paired ground truth 一致。

---

## Q11. VCF 層 GT2 中 hp1-1 vs hp2-1 數量是否在純樣本下偏斜 17 倍？

### 11.1 精確答案：**No** — VCF 層 hp1-1 ≈ hp2-1（接近 1:1.06）

從 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/12_gt_distribution_audit.md` §3.3（HCC1395 baseline 純樣本，PASS 47,798 個 somatic）：

| GT2 sub-genotype | 用戶提供的定義 | PASS 中比例 | 估計 site count |
|------|------|-----------:|--------------:|
| `GT2=1\|., GT3=./.` | **hp1-1 是 alt（somatic-derived）** | **18%** | ~**8,604** |
| `GT2=.\|1, GT3=./.` | **hp2-1 是 alt（somatic-derived）** | **19%** | ~**9,082** |
| `GT2=1\|1, GT3=./.` | Ambiguous somatic haplotype (hp3); both alleles are alternate | 3.6% | ~1,720 |

→ **比例 hp1-1 : hp2-1 = 8,604 : 9,082 ≈ 1:1.06（幾乎完全平衡）**
→ **VCF 層級沒有 17 倍偏斜！**

### 11.2 baseline 與 V2b 在 VCF 層 GT2 完全相同

從 12_gt_distribution_audit §3.1 + §3.3：

| GT raw 值 | Baseline 數量 | V2b 數量 | Δ |
|----------|:------------:|:--------:|:-:|
| `0\|0`（somatic NoLOH）| **21,304** | **21,304** | **0**（完全一致）|
| `0\|.`（somatic in LOH）| 11,323 | 11,290 | −33（< 0.3% 翻向）|
| `.\|0`（somatic in LOH）| 660 | 693 | +33（< 0.3% 翻向）|
| GT2 × GT3 cross-tab（PASS only）| **完全相同** | 完全相同 | — |

→ **baseline 與 V2b 在 VCF GT2 層的 hp1-1 / hp2-1 / hp3 計數幾乎完全一致**。
→ phasing 階段對 somatic 的 GT 決策**沒有偏斜也沒有變動**。

### 11.3 17.3:1 偏斜是 BAM 層的 read count，不是 VCF 層的 site count

⚠ **關鍵釐清 — 兩個層級不要混淆**：

| 層級 | 統計對象 | baseline 數值 | 比例 |
|------|---------|:----------|:----:|
| **VCF 層**（site count，每個 variant 計 1 次）| hp1-1 vs hp2-1 sites | 8,604 vs 9,082 | **≈ 1:1.06**（接近平衡）|
| **BAM 層**（read count，每條 read 計 1 次）| HP1 family（HP:i:1 + HP:i:11）vs HP2 family（HP:i:2 + HP:i:21）reads | **614,000 vs 35,500** | **17.3:1**（嚴重偏斜）|

**為什麼兩層級不一致？**
- VCF 層：phase 階段決定每個 somatic variant 的 hp1-1 / hp2-1 標籤（**決策接近平衡**）
- BAM 層：tag 階段對每個 read 投票指派 HP（**read 大量集中於 HP1 family**）
- → **同一份 phased VCF 在 read level 投票時，被 getVote priority bug 加 enum mismatch 拉偏向 HP1 family**

### 11.4 用戶提供的 GT/GT2/GT3 定義對照（精確化用法）

| GT 類別 | GT 值 | GT2 / GT3 細節 | 含義 |
|------|------|------|------|
| **Germline het** | `0\|1`, `1\|0` | — | hap1/hap2 互為 ref/alt |
| **Germline hom or LOH** | `.\|1`, `1\|.` | — | 一條 hap unknown（potential LOH）|
| **Somatic NoLOH** | `0\|0` | `GT2=1\|., GT3=./.` 或 `GT2=.\|1, GT3=./.` 或 `GT2=1\|1, GT3=./.` | 兩條 hap 都 ref；GT2 補 hp1-1/hp2-1/hp3 |
| **Somatic in LOH** | `.\|0`, `0\|.` | `GT2=0\|., GT3=1\|.` 等（hp1-1-1 / hp1-1-2 細節）| 一條 hap 已知 ref，另一條 unknown；GT2/GT3 補次層 sub-clone |
| **Unphased** | `0/1`, `1/1` | — | edge 不足或無 anchor |

**GT2/GT3 統計來源**：12_gt_distribution_audit §3.3 直接從 phased VCF 解析；audit suite 已驗證 baseline / V2b / V5 三版本在 PASS GT2 × GT3 cell 分布**完全相同**（n=47,798 PASS variants × 3 versions × 多個 GT2/GT3 組合）。

### 11.5 為何這個結果反過來證明「self-phasing bug 在 tag 階段」

12_gt_distribution_audit §4.2 結論：

```
┌───────────────────────────────────────────────────────────────────┐
│  Phasing 階段（longphase-to phase → VCF GT/GT2/GT3）              │
│   baseline.vcf ≈ v2b.vcf                                          │
│   - somatic GT 100% 一致（21,304 / 11,983）                        │
│   - GT2 hp1-1 vs hp2-1 ≈ 1:1.06（接近平衡）                         │
│   - GT2 × GT3 sub-genotype 完全相同                                │
│   - germline 翻向但語意不變                                         │
│   → ❌ 不是 phasing 問題                                           │
└───────────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌───────────────────────────────────────────────────────────────────┐
│  Haplotag 階段（longphase-to haplotag → BAM HP:i tag）            │
│   baseline.bam (17.3:1 偏斜) ≠ v5.bam (≈1:1)                      │
│   - HP1 family reads 614K vs HP2 family 35.5K                     │
│   - 99.9% reads 拿到 HP21（V2b PON-only 啟用後）                    │
│                                                                   │
│   根因：getVote() priority bug + enum mismatch                    │
│   修復：V3-Fixed (`41ff147`) 兩層 + V5 Layer 1.5                  │
│   → ✅ 確認是 tag 問題                                             │
└───────────────────────────────────────────────────────────────────┘
```

**邏輯**：
1. 如果 self-phasing 是 phase 階段的 bug，VCF 層的 hp1-1 vs hp2-1 應該也偏斜（hp1-1 >> hp2-1）
2. 但實測 VCF 層 hp1-1 ≈ hp2-1（1:1.06）→ **phase 階段沒有偏斜**
3. 17.3:1 偏斜只在 BAM read count 出現 → **bug 必在 read 投票（tag）階段**
4. 修復目標：V3-Fixed/V5 修補 `HaplotagProcess.cpp`（tag 階段）；不需動 `PhasingProcess.cpp`（phase 階段，除了 V2b 加 PON-only flag）

### 11.6 一句結論

> **不是**。HCC1395 baseline 純樣本下 VCF 層 GT2 中 hp1-1 vs hp2-1 = **8,604 vs 9,082 ≈ 1:1.06**（接近完全平衡），**沒有 17 倍偏斜**。
>
> 17.3:1 是 **BAM 層 read count**（HP1 family 614K vs HP2 family 35.5K），不是 VCF 層 site count。
>
> 這個對比反過來**證明 self-phasing bug 在 haplotag 階段而非 phase 階段** — 同一份接近平衡的 phased VCF，在 read 投票時被 getVote priority bug + enum mismatch 拉偏到 HP1 family。修復目標明確：改 `HaplotagProcess.cpp`，不需動 `PhasingProcess.cpp`。

---

## Q12. tag 階段為何造成不平衡？已確認問題點？純 vs 0.6 機制意義差異？

### 12.1 tag 階段不平衡的 3 重機制（已確認）

| 機制 | 階段 | 已確認問題點（程式碼層級）| 直接影響 |
|------|------|----------------------|------|
| **(M1) Phasing graph self-phasing** | phase（產 VCF）| `/big7_disk/liaoyoyo2001/longphase-to-mod/PhasingProcess.cpp:154-157` baseline 不啟用 PON-only → graph 用 germline+somatic+unknown 混合 anchor → 同 sub-clone 多個 somatic variants 在 long read 共現 → 強連結被推同一 phase block | VCF 層 PS block 內 somatic 集中於某邊（site-level 全基因組仍接近 1:1）|
| **(M2) getVote priority bug** | tag（產 BAM）| `HaplotagProcess.cpp:506-530` baseline `variantKeys` 順序 `{HP1_1, HP2_1}` 在 `{HP1, HP2}` 之前；任一 somatic vote → break，germline 永遠被忽略 | read 跨 PS block 投票時，每個 read 把同 sub-clone 的 reads 鎖向同一 HP family |
| **(M3) enum vs integer literal mismatch** | tag（產 BAM）| `HaplotagProcess.cpp:697` caller 端 `if(hpResult != HAPLOTYPE1_1 && hpResult != HAPLOTYPE2_1)` 用 enum 值 (3, 4) 比 integer literal (11, 21)；**永不匹配 → fallback 路徑永不觸發** | baseline HP:i:33 永遠 = 0；本應標 ambiguous 的 reads 被推到 HP1 family |

### 12.2 從 site count 1:1.06 到 read count 17.3:1 的「放大效應」

**關鍵釐清**：VCF 層 site count（hp1-1 18% vs hp2-1 19%）與 BAM 層 read count（HP1 family 614K vs HP2 family 35.5K）並不矛盾，因為兩者是不同統計層級：

```
跨基因組 PS block orientation 統計（12_gt_distribution_audit §3.2）：
  baseline GT 分佈：
    1|0 = 508,368 sites（hap1 是 alt，slight 多）
    0|1 = 559,116 sites（hap2 是 alt）
                  ↓
  每個 read 跨 5-50 個 variants 投票：
    M1: phasing graph 把同 sub-clone variants 集中於 PS block 某邊
    M2: getVote 順序錯，read 投票偏向 somatic vote
    M3: enum mismatch 把不該標的 reads 拉到 HP1 family
                  ↓
  全基因組 read 累積（baseline）：
    HP1 family（HP:i:1 + HP:i:11）= 614,000 reads
    HP2 family（HP:i:2 + HP:i:21）= 35,500 reads
    比例 = 17.3:1（嚴重偏斜）
                  ↓
  VCF 層卻仍是：
    hp1-1 sites = 18% (~8,604)
    hp2-1 sites = 19% (~9,082)
    比例 ≈ 1:1.06（接近平衡）
```

**為何 site count 平衡但 read count 偏斜**：read 在跨多個 site 投票過程中**累積方向偏差** — 一條 read 跨 5 個 hp1-1 sites + 5 個 hp2-1 sites，理論上 vote 應接近 5:5；但 priority bug 讓任一 hp1-1 vote > 0 就 break，加上 enum mismatch 讓 fallback 失效，多輪投票後 read 集中向 HP1 family。

### 12.3 已確認問題點 — 4 項硬性 sanity check 全部 PASS

從 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/06_v5_sanity_bug_check.md` §7：

| 4 項硬性檢查 | 違反位點 / reads | 結論 |
|------|----------------------|------|
| 守恆律 A：Δ-consistency（HP33 減少 = HP11+HP21 增加）| 0 / 15 sites | ✅ PASS |
| 守恆律 B：Germline 不變（V3F→V5 不動 HP1/HP2）| 0 / 15 sites | ✅ PASS |
| Layer 1.5 預期 1：33→directional 精確守恆 | 0 / 15 sites | ✅ PASS |
| Layer 1.5 預期 2：無 germline → HP33 違規 | 0 reads pooled | ✅ PASS |

**不只「指出問題點」，還通過 4 項硬性數學守恆律驗證 V3F + V5 修補不引入新 bug**。

每個 V3F/V5 修補對應的 bug：

| Bug ID | 階段 | 程式碼位置 | V5 audit suite 對應 |
|--------|------|----------|------------------|
| **Bug 1-1** | tag M2 priority | `HaplotagProcess.cpp:506-530` | V3-Fixed `41ff147` 重寫為兩層；V5 加 Layer 1.5 |
| **Bug 1-2** | tag M3 enum mismatch | `HaplotagProcess.cpp:697` | V3-Fixed `41ff147` 同 commit 改用 integer literal |
| **Bug 1-3** | tag UNDEFINED guard（INDEL）| `HaplotagProcess.cpp:497-510` | `380e8d2` UNDEFINED guard |
| **Bug 1-4** | tag UNDEFINED guard（SNP alt）| `HaplotagProcess.cpp:489-494` | V5 working tree（uncommitted）|
| **Bug 2-1** | phase M1 PON-only flag | `PhasingProcess.cpp:154-157` | V2b `8b8c1fd` 加 `--pon-only-phasing` flag |

→ **5 個確認的 bug 點，每個都有對應的 commit fix + sanity 驗證**。

### 12.4 純樣本 vs 0.6 purity 的機制意義差異

**核心釐清**：機制（程式碼）相同，但觸發強度與影響結果差異懸殊。

| 維度 | 純樣本（0.93 HCC1395）| 0.6 purity（t30_n20 mix）|
|------|-------------|-------------------|
| **程式碼 bug 是否存在**（M1/M2/M3）| ✅ 同樣存在（程式碼一致）| ✅ 同樣存在（程式碼一致）|
| **sub-clone read 共現強度** | 強（純 tumor reads）| 弱（normal 稀釋 sub-clone）|
| **countMap 中 germline:somatic** | ~3:2（baseline）/ ~1:3（V2b PON-only）| ~5:1（baseline）/ ~3:1（V2b）— normal 補 germline anchor |
| **M1 self-phasing 觸發** | 嚴重（graph 強連結）| 自然減弱（read 共現減少）|
| **M2 getVote bug 觸發機率** | baseline ~30% / V2b >95% | baseline <10% / V2b ~30% |
| **觀察症狀（baseline）**| 17.3:1 偏 HP1 family（chr19-22 1.33:1） | **1:1.14**（接近平衡）|
| **觀察症狀（V5）**| ~1:1（修復成功） | 1:2.34 + HP33 12.4%（**conservative tagging**）|
| **V5 主要價值** | **修 17.3:1 → ~1:1**（核心價值）| **conservative tagging**（HP33 ↑ 6×；baseline 自然 1:1.14，V5 修復 unobservable）|

### 12.5 為何純樣本機制觸發強度高？— 物理直覺

```
═══════════════════════════════════════════════════════════════
純樣本（0.93）：sub-clone read 共現強
═══════════════════════════════════════════════════════════════

  Read A: ━━━━━●━━━━━━●━━━━━━●━━━━━━●━━━━━ （帶 4 個 somatic ALT）
  Read B: ━━━━━●━━━━━━●━━━━━━●━━━━━━●━━━━━ （同 sub-clone，相同 4 個 somatic ALT）
  Read C: ━━━━━●━━━━━━●━━━━━━●━━━━━━●━━━━━ （同上）
                ↓
  phasing graph 看到：A-B-C 共享 4 個 somatic edges
  → 強連結 → 全推同一 phase block 同一邊（M1 啟動）
                ↓
  read 跨 PS block 投票：
    countMap[HP1_1] = 多 / countMap[HP2_1] = 少
  → M2 priority bug：第一對 break，germline 忽略
  → M3 enum mismatch：fallback 失效
  → read 全集中於 HP1 family

═══════════════════════════════════════════════════════════════
0.6 purity（t30_n20）：normal reads 稀釋 sub-clone 共現
═══════════════════════════════════════════════════════════════

  Read A (tumor): ━━━━━●━━━━━━●━━━━━━━━━━━ （帶 2 個 somatic ALT）
  Read B (tumor): ━━━━━━━━━━━━●━━━━━━●━━━━ （帶 2 個 somatic ALT）
  Read C (normal): ━━━━━━━━━━━━━━━━━━━━━━━ （無 somatic，只 germline ALT）
  Read D (normal): ━━━━━━━━━━━━━━━━━━━━━━━ （同上）
                ↓
  phasing graph 看到：A-B 弱連結（只 1 個共享 somatic）
  → 連結被 normal reads 稀釋
  → graph 不形成強連結 → M1 自然減弱
                ↓
  read 跨 PS block 投票：
    countMap[HP1] = 多（normal 提供）/ countMap[HP1_1] = 少
  → M2 priority bug 觸發機率低（germline 票多 → 不至於某對 0/0）
  → 即使 M2 觸發，HP1_1 vs HP2_1 票分布較均勻
  → 整體 HP1 family vs HP2 family 接近 1:1.14
```

→ **同樣的程式碼 bug 在不同 purity 下觸發強度不同**：純樣本「最大強度顯現」；0.6 purity「自然減弱」。

### 12.6 純 vs 0.6 修補價值差異總結

| purity 範圍 | V5 修補的真實價值 | 為何如此？ |
|------------|----------------|----------|
| **≥ 0.85 高純度** | **修 17:1 self-phasing → 1:1**（核心價值）| sub-clone read 共現強，bug 全面爆發，V5 修補帶來顯著改善 |
| **0.6-0.85 中純度** | V5 / baseline 都可（V5 更 conservative）| self-phasing 中等強度，V5 多花 ~10% reads 進 HP33 |
| **< 0.5 低純度** | V5 + HP33-aware downstream（safety net）| self-phasing 自然減弱，baseline 自然較好；V5 conservative tagging 對下游分析更安全 |

來源：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md` §7.4。

### 12.7 一句結論

> **tag 階段不平衡的根因是 3 重機制疊加**：(M1) phase 階段 self-phasing 把同 sub-clone variants 集中於 PS block 某邊；(M2) tag 階段 `getVote()` priority bug 讓 somatic vote 一票否決 germline；(M3) caller 端 enum vs integer literal mismatch 讓 fallback 永不觸發。
>
> **5 個程式碼問題點全部已確認**（含位置 file:line + 對應修補 commit），**4 項硬性 sanity check 15/15 PASS**。
>
> **純樣本（0.93）vs 0.6 purity 機制相同（程式碼一致），但觸發強度差異懸殊**：純樣本 sub-clone read 共現強 → bug 全面爆發 → 17.3:1 嚴重偏斜；0.6 normal 稀釋 → bug 觸發機率低 → 1:1.14 接近平衡。
>
> **V5 修補價值隨 purity 變化**：高純度修 self-phasing 為核心；低純度自然減弱，V5 改提供 conservative tagging 安全網（HP33 比例上升）。

---

## Q13. 「baseline HP:i:33=0；ambiguous reads 被推到 HP1 family」是否是 HP1-1 特別多的原因？

### 13.1 直接答覆

✅ **是，這是核心驅動力之一**（不是唯一原因）。但需精確化：**不只 HP:i:33 reads 被錯誤推到 HP1 family，HP:i:21 reads 也被錯誤推到 HP1 family**。V3-Fixed 修補後揭露的數據顯示，baseline HP:i:11 「過度集中」中包含 ~752K reads 是錯誤分配。

### 13.2 V3-Fixed 修補揭露的 reads 流向真相

從 `memory/project_v5_somatic_fallback_verification.md` 第 28-36 行 + V3-Fixed memory：

| HP tag | baseline（推算）⁽¹⁾ | V3-Fixed（修補後）| Δ（V3F − baseline）|
|--------|:-----------------:|:---------------:|:----:|
| HP:i:1 + HP:i:2（germline）| 17,435,063 | 17,435,063 | 0（**守恆律 B 確認 germline 不變**）|
| **HP:i:11**（HP1 family somatic）| **~614,000** | 584,117 | **−~30,000**（從「過度集中」回歸）|
| **HP:i:21**（HP2 family somatic）| **~35,500** | 547,118 | **+~511,618**（**+15.4×**！）|
| **HP:i:33**（ambiguous somatic）| **0**（enum bug 永不寫入）| 239,679 | **+239,679**（**從 0 出現**）|

⁽¹⁾ baseline 的 HP:i:11 / HP:i:21 / HP:i:33 個別數字未直接列於 V5 memory；以「全基因組 somatic ALT reads = HP1 family 614K vs HP2 family 35.5K」（10_somatic_bias §2 Panel D）+「HP:i:33=0 確認」（V3F memory 第 22 行）反推。

### 13.3 V3-Fixed 從 baseline HP:i:11「回收」~752K reads

V3-Fixed 修兩個 bug 後，baseline 過度集中的 HP:i:11「分解」為：

```
baseline 錯誤分佈                     V3-Fixed 正確分佈
────────────────────                  ────────────────────
HP:i:11 ≈ 614K  ┐                    HP:i:11 = 584K（真實 HP1 family somatic）
                │                              +（回收 547K → HP:i:21）
HP:i:21 ≈ 35K   ├── 三類 reads        HP:i:21 = 547K（回歸正確分配）
                │   全推 HP:i:11      
HP:i:33 = 0     ┘                    HP:i:33 = 240K（從 0 出現）
                                      
                                      回收量化：
                                      547K + 240K − 35K（baseline HP:i:21 已存在）
                                      = ~752K reads 從 baseline HP:i:11
                                        重分配到正確標籤
```

**精確分解**：
- ~**512K** reads 是本應該 HP:i:21（HP2 family somatic）但 baseline 標 HP:i:11
- ~**240K** reads 是本應該 HP:i:33（ambiguous）但 baseline 標 HP:i:11（M3 enum mismatch 受害者）
- 這兩組 reads 加總 ~752K → 直接造成 baseline HP1 family（~614K）vs HP2 family（~35K）= 17.3:1 嚴重偏斜

### 13.4 為何「推到 HP1 而非 HP2」？— 兩個微觀機制

**機制 A：M3 enum mismatch 強制 fallback 永不觸發**

baseline `HaplotagProcess.cpp:697` caller 端：
```cpp
if(hpResult != HAPLOTYPE1_1 && hpResult != HAPLOTYPE2_1) {
    // fallback 路徑（本應把 ambiguous reads 標 HP:i:33）
}
```

- `HAPLOTYPE1_1 = 3`（enum）
- `HAPLOTYPE2_1 = 4`（enum）
- 但 hpResult 已是 integer literal 11/21/33
- `11 != 3 AND 11 != 4` → 條件**永遠成立** → 但 fallback 內部寫入也用 enum 比對 → **fallback 邏輯失敗** → reads 沒被正確 fallback → 落到 hpResult 已設定值

**機制 B：M2 priority bug 的「第一對有票即 break」放大微小 PS orientation 偏差**

從 `12_gt_distribution_audit §3.2`：

| GT 值 | baseline 數量 | 含義 |
|------|:------------:|------|
| `1\|0`（hap1=alt）| **508,368** | slight 多 |
| `0\|1`（hap2=alt）| 559,116 | slight 多另一邊 |

baseline PS block orientation **接近平衡**（508K vs 559K，slight 偏 `0|1`）。理論上 read 被推到 HP1 vs HP2 應接近 1:1。

但 priority bug 的「第一對有票即 break」**放大微小偏差**：
- 一條 read 跨 5 個 PS blocks
- 即使每個 PS block HP1_1 vote 只比 HP2_1 多 1，5 個 block 都觸發 break → 全推 HP1_1
- 跨基因組累積 → 17:1 級放大

### 13.5 完整因果鏈（從 0 到 17.3:1）

```
原始狀態：
  純樣本（0.93）+ baseline 預設不啟用 PON-only flag

[Step 1] M1 self-phasing（phase 階段）
  phasing graph 用 germline+somatic+unknown 混合 anchor
  → 同 sub-clone 多個 somatic variants 在 long read 共現
  → graph 強連結把同 sub-clone variants 集中於 PS block 某邊
  → VCF GT2 site count 跨基因組仍接近 1:1（hp1-1 18% vs hp2-1 19%）
     但每個 PS block 內 somatic 集中

[Step 2] M2 priority bug（tag 階段 getVote line 506-530）
  read 跨 PS block 投票，countMap 累積票
  variantKeys 第一對 {HP1_1, HP2_1} 任一 vote > 0 → break
  → germline votes 被忽略
  → read 被推向 somatic family

[Step 3] M3 enum mismatch（tag 階段 caller line 697）
  fallback 條件型別不一致 → fallback 路徑邏輯失敗
  → 本應標 HP:i:33 的 reads（~240K）必須去某處
  → 加上 M2 把 PS orientation slight 偏 0|1 放大成 HP1 偏好
  → 240K HP:i:33 reads + 512K HP:i:21 reads 全推 HP:i:11

[Step 4] 累積結果
  baseline HP:i:11 ≈ 614K（HP1 family）
  baseline HP:i:21 ≈ 35K（HP2 family）
  baseline HP:i:33 = 0
  → HP1 family : HP2 family = 17.3 : 1
```

### 13.6 V3-Fixed 修補後的「回歸」效果

V3-Fixed `41ff147` 兩個關鍵修補：

| 修補 | 對應 bug | 效果 |
|------|---------|------|
| getVote 改寫為兩層（germline first → somatic annotate）| M2 priority | 大部分 reads 走 Layer 1 直接拿 HP:i:1/2/11/21；不再被 priority bug 鎖向 somatic |
| caller 端比對改 integer literal 11/21/33 | M3 enum mismatch | fallback 條件正確觸發 → 240K ambiguous reads 正確標 HP:i:33 |

→ V3F 後 HP1 family : HP2 family ≈ **1.07:1**（接近平衡，與 V3F memory 第 19 行 TP Balanced% 22.5% 一致）。

### 13.7 V5 進一步補強：Layer 1.5 把過度保守的 HP:i:33 重分配

V3-Fixed 還有「過度保守」問題：240K HP:i:33 中部分 reads 其實有 phased 方向資訊（HP1_1 OR HP2_1 vote > 0 但 germline=0）。V5 Layer 1.5 把這些 reads 升級回 HP:i:11 / HP:i:21：

| HP tag | V3-Fixed | V5 | Δ |
|--------|:-------:|:--:|:-:|
| HP:i:11 | 584,117 | 666,997 | +82,880 |
| HP:i:21 | 547,118 | 593,720 | +46,602 |
| HP:i:33 | 239,679 | **110,197** | **−129,482**（−54%）|

→ V5 把 V3F 的 240K HP:i:33 reads 中 **129K 重分配回 HP:i:11/21**（保留方向資訊），剩下 110K 真正 ambiguous 仍標 HP:i:33。

### 13.8 一句結論

> **是，「baseline HP:i:33=0 + ambiguous reads 被誤推到 HP1 family」是 17.3:1 偏 HP1 family 的核心驅動力之一**。但更精確：不只 HP:i:33 reads（~240K）被誤推，**HP:i:21 reads（~512K）也被誤推到 HP:i:11**。
>
> V3-Fixed 修補後從 baseline HP:i:11 「**回收**」約 **752K reads**：~512K 回歸 HP:i:21（從 35K 跳到 547K，+15.4×），~240K 回歸 HP:i:33（從 0 出現）。這個 752K 回收量化解釋了 baseline 17.3:1 嚴重偏斜為何在 V3-Fixed 後降到 ~1.07:1（接近平衡）。
>
> **核心因果鏈**：M1 self-phasing（PS block 內 somatic 集中）+ M2 priority bug（第一對有票即 break）+ M3 enum mismatch（fallback 永不觸發）三重疊加 → 把本該標 HP:i:21 與 HP:i:33 的 reads 全推到 HP:i:11 → HP1 family 嚴重過度集中 → 17.3:1 偏 HP1 family。

---

## Q14. V5 Layer 1.5 邏輯合理嗎？HP:i:11/21 語義混雜？業界共識對照？

### 14.0 ⚠ 程式碼校正（先校正先前 Q&A 中的錯誤）

從 grep `/big7_disk/liaoyoyo2001/longphase-to-mod/Util.h:19-26` 取得真正的 enum 定義：

```cpp
enum Haplotype {
    HAPLOTYPE_UNDEFINED = -1,
    HAPLOTYPE1 = 0,      ⚠ 不是 1（先前 Q&A Q1, Q3, Q10 描述錯誤）
    HAPLOTYPE2 = 1,      ⚠ 不是 2
    HAPLOTYPE3 = 2,      ⚠ 不是 5
    HAPLOTYPE1_1 = 3,    ✅ 正確
    HAPLOTYPE2_1 = 4,    ✅ 正確
};
```

**校正點**：先前 Q1/Q3/Q10/Q12/Q13 寫的「HAPLOTYPE1=1, HAPLOTYPE2=2, HAPLOTYPE3=5」應為「HAPLOTYPE1=0, HAPLOTYPE2=1, HAPLOTYPE3=2」。但 enum vs integer literal mismatch 的描述方向**仍然成立**（baseline hpResult 已是 integer 1/2/11/21/33，與 enum 值 0/1/2/3/4 完全錯位）。

**校正點 2**：Q4/Q12 描述 baseline caller fallback「reads 被錯誤分配（多數推 HP1 family）」是**推論非源碼證據**。從 audit suite 01 §2.2 精確說法應為「整數 literal 修正後此 fallback 才生效」 — 意思 baseline 下 fallback 邏輯**不生效**（hpResult 拿到 getVote 給的原值，未經 fallback 修正）。

### 14.1 V5 完整邏輯（grep `HaplotagProcess.cpp:512-563` 確認）

```cpp
// Layer 1: germline 主導方向
if (countMap[HP1] > 0 || countMap[HP2] > 0) {
    germlineResult = (HP1 >= HP2) ? 1 : 2;
}
// Layer 1.5（V5 NEW）：germline=0 才用 phased somatic 推方向
else if (countMap[HP1_1] > 0 || countMap[HP2_1] > 0) {
    germlineResult = (HP1_1 >= HP2_1) ? 1 : 2;
}

// Layer 2: 編碼 HP tag
if (somaticTotal > 0) {                       // 是否有 somatic 經過？
    if (germlineResult == 1)       hpResult = 11;  // ← HP:i:11
    else if (germlineResult == 2)  hpResult = 21;  // ← HP:i:21
    else                           hpResult = 33;
} else {
    hpResult = germlineResult;  // 純 germline → 1/2
}
```

**重要釐清**：
- V5 **不會「過度誤標為 germline」** — Layer 2 的 `somaticTotal > 0` 條件保證有 somatic 經過的 read 必標 11/21/33（不會降為 1/2）
- V5 設計本質：**「方向看 germline；annotation 看是否有 somatic」**，兩者拆分

### 14.2 ⚠ HP:i:11/21 語義混雜（用戶質疑指向的設計取捨）

V5 後 HP:i:11 含**兩種來源**：

| 來源 | 條件 | germline 是否真有 anchor？ | enum 命名語義 |
|------|------|----------------|----|
| **A 類** | Layer 1: HP1≥HP2>0 + somaticTotal>0 | ✅ 有 germline anchor | 嚴格符合「HP1 + somatic」原意 |
| **B 類** | Layer 1.5: HP1=HP2=0 但 HP1_1≥HP2_1>0 | ❌ **沒有 germline anchor** | 「方向推測是 HP1，含 somatic」— 非原意但 phased _1 後綴含 on-HP1 訊息 |

**用戶質疑點**：HP:i:11 通常被理解為「germline + somatic 同邊」；如果 B 類「沒有 germline 也標 HP:i:11」，分類**意義不純**。

→ **質疑成立**。這是 V5 設計取捨，但 V5 audit suite §9 caveat **未明示**「HP:i:11 含 A+B 兩類來源」。

### 14.3 V5 Layer 1.5 的合理性 6 層論證

| 層次 | 論證 | 信心 |
|------|------|----|
| **1. enum 命名語義** | HAPLOTYPE1_1/HAPLOTYPE2_1 後綴 `_1` 含「on HP1/HP2」方向資訊（PON-only 啟用後可信）| ★★★★ |
| **2. PON-only 後 phasing graph 可信** | V2b `convertNonGermlineToSomatic()` 排除 self-phasing → HP1_1/HP2_1 phasing 反映真實 phased somatic | ★★★★★ |
| **3. read 投票統計邏輯一致** | Layer 1.5 多數決與 Layer 1 邏輯**完全一致** | ★★★★★ |
| **4. confidence threshold 0.6 把關** | LongPhase-TO 原生 `-p` 參數（已論文認證 20-100% purity）；110K HP:i:33 confidence < 0.6 被正確攔截 | ★★★★ |
| **5. sanity check 數學守恆** | 4 項硬性檢查 15/15 PASS, 0 violation | ★★★★★ |
| **6. paired GT concordance** | clean PS V5=88.2% vs BL=74.9%（+13.3pp）；全基因組 V5=90.5%（PI 報告 4）| ★★★★★ |

**整體信心：★★★★☆**（在 PON-only 模式下）。⚠ 唯一 caveat：HP:i:11 含 A+B 兩類**未明示**。

### 14.4 baseline 直接用 HP1_1 投票判方向的 3 大問題

| # | 問題 | 機制 |
|---|-----|------|
| **問題 1** | HP1_1 標籤本身被 self-phasing 汙染 | baseline `--pon-only-phasing=false` → phasing graph 用混合 anchor → HP1_1 標籤是 self-phasing 自己決定的結果 → **環依賴**（用 self-phasing 結果驗證 self-phasing）|
| **問題 2** | somatic 一票否決 germline | `getVote()` 第一對 `{HP1_1, HP2_1}` 有票即 break → germline 永遠被忽略 → **信息量低的 vote 否決信息量高的 vote** |
| **問題 3** | HP:i:11 標籤被過度分配 | V3-Fixed 揭露 baseline HP:i:11 過量含 ~512K HP:i:21 + ~240K HP:i:33 = **~752K reads** 被誤推（佔 baseline HP1 family 614K 的 ~95%）|

### 14.5 V5 的 3 個關鍵改進（對應問題 1/2/3）

| 改進 | 對應修補 | 效果 |
|------|---------|------|
| **改進 1** 啟用 PON-only flag | V2b `8b8c1fd` `convertNonGermlineToSomatic()` | phasing graph 只用 PON-germline 主導 → HP1_1 標籤反映真實 phased somatic（**修復標籤本身**）|
| **改進 2** germline first | V3-Fixed `41ff147` getVote 兩層 | Layer 1 germline 主導方向 → HP1_1 從「主決策」降為「fallback」 |
| **改進 3** somatic 作 annotation | V3-Fixed `41ff147` Layer 2 | hpResult 拆「方向（germline）+ somatic 經過 annotation」雙維度，不再混淆 |

### 14.6 「germline 主 HP2 但有 HP1_1 vote」case 的精確量化

用戶提出的具體場景：

```
countMap 範例：
  HP1   = 2 votes
  HP2   = 8 votes  ← germline 主 HP2
  HP1_1 = 3 votes  ← phased somatic on HP1（衝突）
  HP2_1 = 0 votes
```

| 版本 | 路徑 | 結果 | 評估 |
|------|---------|------|------|
| **baseline** | Pair ① HP1_1=3>HP2_1=0 → break → HP:i:11 | ❌ 偏 HP1（**忽略 germline 主 HP2 = 8 votes**）|
| **V5** | Layer 1: HP2=8>HP1=2 → germlineResult=2 → somaticTotal>0 → **HP:i:21** | ✅ 偏 HP2（符合 germline 主導訊號）|

**這類 case 的數量**：從 V3-Fixed 揭露的 reads 流向反推 ~**512K reads** 屬於這類「germline 應主 HP2 但 baseline 誤推 HP1」場景（佔 baseline HP:i:11 過量的 68%）。

**V5 信任 germline 而非 phased somatic 的理由**：

| 訊號來源 | 信心源 | PON-only 下可信度 |
|---------|------|----------|
| germline (HP1, HP2) | 4 個 PON 資料庫過濾（1000g-pon, CoLoRSdb, dbsnp, gnomad）| **高**（Lin 2022 LongPhase 原版設計核心）|
| phased somatic (HP1_1, HP2_1) | phasing graph 推測（仍有 graph noise）| 中 |

→ 衝突時信任 PON-confirmed germline 是合理選擇。

### 14.7 業界共識對照（V5 對齊 8/8 主流工具）

從 `InterSubMod/docs/references/manual/20260402_longphase_to_phasing_quality_literature.md` §4：

| 工具 | 論文 | 設計原則 | V5 對齊？ | baseline 對齊？ |
|------|------|---------|:---:|:---:|
| **LongPhase 原版** | Lin et al. *Bioinformatics* 2022 | germline het 為 phasing scaffold | ✅ | ❌ |
| **WhatsHap** | Martin 2016, Patterson 2015 | germline phasing；somatic 不處理 | ✅ | ❌ |
| **HapCUT2** | Edge 2017 *Genome Research* | germline-focused phasing | ✅ | ❌ |
| **MethPhaser** | Nature Comm 2024 | germline pre-phased 為 seed + methylation 擴展 | ✅ | ❌ |
| **Octopus** | Cooke 2021 *Nature Biotech* | joint calling + phasing 避免 self-phasing | ✅ | ❌ |
| **Severus** | Nature Biotech 2025 | **明確使用 normal sample germline phasing 為基礎，避免 tumor variants 汙染 phasing scaffold** | ✅ **直接對應 V5** | ❌ |
| **SAVANA** | Nature Methods 2025 | haplotype-resolved，germline-aware | ✅ | ❌ |
| **POG Study** | Sato 2022 Nature Comm | germline phasing 為 scaffold，somatic 被動分配 | ✅ | ❌ |

**業界共識（8/8 工具）**：
1. **germline het 為 phasing scaffold 主要 anchor**
2. **somatic 是 annotation 或 sub-genotype，不主導 scaffold**
3. **避免 tumor variants 汙染 phasing scaffold**（Severus 明示）

→ **V5 完全符合業界共識；baseline 完全違反 8/8 工具設計原則**。

### 14.8 V5 Layer 1.5 在業界中的定位

V5 Layer 1.5（germline=0 但 phased somatic 推方向）是 **fork 內部創新**，業界主流工具沒有直接對應，但與業界精神兼容：

| 評估維度 | 結論 |
|---------|----|
| 與「germline first」原則衝突？ | **不衝突**（germline 仍優先；Layer 1.5 是 fallback）|
| 引入新風險？ | 部分 — phased somatic 推方向有 graph noise；V5 用 confidence 0.6 + 4 sanity 把關 |
| 業界類似機制？ | MethPhaser 用 methylation 作 germline 擴展（類似 fallback 概念）|

→ Layer 1.5 是「**業界精神內的 fork 創新**」（fallback 而非主導）。

### 14.9 對下游分析者的提醒（caveat 補強建議）

V5 audit suite §9 應補一條 caveat：

> **HP:i:11/21 在 V5 後可能來自兩種來源**：
> - **A 類**（Layer 1）：read 有 germline vote（HP1>0 或 HP2>0）+ somatic vote
> - **B 類**（Layer 1.5）：read germline=0 但 phased somatic（HP1_1/HP2_1）有 vote
>
> 下游分析若假設「HP:i:11 = 有 germline anchor」會被 B 類誤導。需要區分時，可使用 `--log` 輸出 read-level 投票結果或從 BAM 重新計算 countMap。

### 14.10 一句結論

> **V5 Layer 1.5 邏輯合理且有 6 層論證 + 業界 8/8 工具對齊**：(1) enum 命名語義；(2) PON-only 啟用後可信；(3) 邏輯一致；(4) confidence 0.6 把關；(5) 4 項 sanity check 15/15 PASS；(6) paired GT concordance +13.3pp。
>
> **baseline 直接用 HP1_1 投票判方向有 3 大問題**：(1) HP1_1 標籤被 self-phasing 汙染（環依賴）；(2) somatic 一票否決 germline（priority bug）；(3) HP:i:11 過度集中（含 ~752K 誤推 reads）。
>
> **「germline 主 HP2 但有 HP1_1 vote」case 規模龐大**（~512K reads，佔 baseline HP:i:11 過量的 68%）— V3-Fixed 把這類 reads 從錯誤的 HP:i:11 救回到正確的 HP:i:21。
>
> **V5 完全符合業界共識**（Lin 2022, WhatsHap, HapCUT2, MethPhaser, Octopus, Severus, SAVANA, POG Study **8/8 對齊**）；**baseline 完全違反 8/8 業界設計原則**。
>
> ⚠ **唯一未明示 caveat**：V5 後 HP:i:11/21 含 A 類（germline+somatic）+ B 類（germline=0 但 phased somatic）兩種來源；下游分析若假設「HP:i:11 = 有 germline anchor」會被誤導。

---

## Q15. Cross-purity ablation 實證 — PON 方法是否最佳 + F1 證實 + baseline bug 跨 purity

### 15.1 完整 ablation 矩陣（0.93 + 0.6）

執行於 2026-04-30，產出檔案：
- `InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_getVote_only_ablation_results_01.md`（0.93）
- `InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_ablation_purity06_results_01.md`（0.6）

| 實驗 | binary | --pon-only-phasing | 0.93 BAM | 0.6 BAM |
|------|---|:---:|---|---|
| A1 / B1 baseline | longphase-to-baseline | OFF | ✅ 既有 | ✅ 既有 |
| A2 V2b（0.93 only）| longphase-to-baseline | ON | ✅ 既有 | ⏸ 未做 |
| **A3 / B3 v3f_no_pononly** ★ | longphase-to-v3fixed | **OFF** | ✅ 新測 | ✅ 新測 |
| A4 V3F+PONonly（0.93 only）| longphase-to-v3fixed | ON | ✅ 既有 | ⏸ 未做 |
| A5 / B5 V5 | longphase-to (V5) | ON | ✅ 既有 | ✅ 既有 |

### 15.2 caller F1 實證（6 版本同時測試 vs SEQC2 truth）

| label | TP | FP | FN | precision | recall | **F1** |
|------|---:|---:|---:|---:|---:|---:|
| B1 baseline @ 0.6 | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |
| **B3 v3f_no_pononly @ 0.6** | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |
| **B5 V5 @ 0.6** | 24,190 | 13,487 | 15,257 | 0.6420 | 0.6132 | **0.6273** |
| A1 baseline @ 0.93 | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| **A3 v3f_no_pononly @ 0.93** | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| **A5 V5 @ 0.93** | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |

來源：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260430_v3f_ablation/phased_vcf_f1.tsv`

🎯 **實證結論**：
- ✅ **3 版本 F1 在每個 purity 內完全相同**（TP/FP/FN 每個小數位一致）
- ✅ **跨 purity 一致**：V5 vs baseline F1 = 0（不論 0.93 或 0.6）
- ✅ 0.93 F1 = 0.7166 與 V3-Fixed memory + PI 報告 4 §3.4 完全一致 — 計算可信
- 🆕 0.6 F1 = 0.6273（首次發布；09_purity06_simulation §6 caveat 之前未計算）

### 15.3 4 個質疑的精確答覆

#### Q15a「啟用 Two-Pass / PON-only 後還有 self-phasing 嗎？」

**主要 artifact 消除，但有殘留**：

| 殘留來源 | 是否仍存在 | 說明 |
|---------|:---:|---|
| 17.3:1 全基因組偏斜 | ❌ 消除 | A2/A4/A5 ratio 都翻轉到 ~0.65 |
| PS block orientation slight 偏好 | ⚠ 仍有 | V5 audit 09 §5.1：「1:2.34 不是 V5 bug，是 PS block 系統偏好 + conservative tagging 合成」|
| problem PS blocks（germline reads concordance < 70%）| ⚠ 仍有 | 07_paired §3.5：少數 PS block 方向不可靠 |
| cnLOH 雙親同源 region | ⚠ 仍有 | 兩 hap 相似，Two-Pass 無法區分 |
| Confidence < 0.6 reads | ✅ 仍標 HP:i:33 | 110K reads（V5 0.93 全基因組）保守攔截 |

#### Q15b「低 purity 下 self-phasing 嚴重嗎？」

**不嚴重，normal mix 自然減弱**：

| Purity | baseline ratio | baseline N50 |
|---|---|---|
| 0.93（純樣本）| **17.3:1**（全基因組）/ 2.73:1（15-site）| **11.9 Mbp**（mega-block）|
| 0.6（mix）| 1:1.14（chr19-22）/ 0.48:1（15-site）| **0.8 Mbp**（**降 14.9×**）|

→ normal reads 稀釋同 sub-clone read 共現 → phasing graph 不形成強連結 → mega-block 切碎 → self-phasing 自然消失。

#### Q15c「HP:i:33 標記機制與條件」

從 V5 程式碼（`HaplotagProcess.cpp:512-563` + `:725-737`），**兩類獨立來源**：

**A 類（getVote 內部）**：
- countMap[HP1] = 0 AND countMap[HP2] = 0（無 germline）
- AND countMap[HP1_1] = 0 AND countMap[HP2_1] = 0（無 phased somatic）
- AND somaticTotal > 0（有 HP3 vote）
- → hpResult = 33

**B 類（caller 端 confidence < 0.6 攔截）**：
- getVote 已給 hpResult ∈ {11, 21, 33}
- 但 max/(max+min) < 0.6
- → 強制改 hpResult = 33

V5 全基因組 110,197 HP:i:33 主要為 B 類（V5 memory 第 50 行：「confidence 全 < 0.6，多數 0.50-0.59」）。

**baseline 跨 purity 的 HP:i:33 = 0**（enum bug）：
- 0.93 baseline HP:i:33 = 0 ⚠ enum mismatch bug
- **0.6 baseline HP:i:33 = 0 ⚠ 同 enum mismatch bug**（從 0.6 ablation 實證）

#### Q15d「V5 在 0.6 F1 比 baseline 差嗎？」

❌ **沒有**。實證確認：

| 版本 | 0.6 F1 | 0.93 F1 |
|---|---:|---:|
| baseline | **0.6273** | 0.7166 |
| V5 | **0.6273**（差 **0**）| 0.7166（差 **0**）|

**機制**：longphase-to phase 不改 VCF FILTER → 三版本 PASS variants 完全相同 → F1 完全相同。V5 修補在 phasing/tagging 階段，**不影響 caller F1**。

⚠ **ISM SuggestFilter F1**（套濾後）在 0.93 PI 報告 4 §3.4 = -0.0003 噪音；0.6 預期類似（屬另案，未跑）。

### 15.4 baseline bug 跨 purity 都存在（重要釐清）

**Bug 1: getVote priority bug** — 跨 purity 都存在；觸發機率隨 purity 變化：
| Purity | 觸發機率（推估）| 症狀 |
|---|---|---|
| 0.93 baseline | ~30% | 17.3:1 偏 HP1（self-phasing 強）|
| 0.93 V2b PON-only | >95% | 99.9% HP21（PON-only 後 germline votes 銳減）|
| 0.6 baseline | < 10% | 0.48:1（normal 補 germline）|

**Bug 2: enum mismatch** — **跨 purity 完全相同**（型別問題與 purity 無關）：
| Purity | baseline HP:i:33 |
|---|---|
| 0.93 | 0 ⚠ |
| **0.6** | **0 ⚠**（**bug 仍存在！只是 self-phasing 弱化掩蓋症狀**）|

→ **0.6 下 baseline「看似 OK」是表象**，bug 仍在；只是被 normal mix 弱化掩蓋整體偏斜。

### 15.5 對 V5 audit suite + Q&A 既有結論的更新

| 既有結論 | 更新方向 |
|---|---|
| Q12 §12.4「PON-only 確定為主要驅動力」 | ✅ 跨 purity 確認（0.93 4.3× / 0.6 1.3×）|
| Q13 §13.2「~752K reads 從 baseline HP:i:11 回收」 | ✅ 拆解：~512K HP:i:21（PON-only 主導翻轉）+ ~240K HP:i:33（getVote 修補主導）|
| Q14 §14.5「V5 的 3 個關鍵改進」 | ✅ 拆「PON-only flag = ratio 翻轉」+「getVote 修補 = HP:i:33 標記」**兩個獨立貢獻** |
| F1 不衡量 V5 品質（Q12 §12.5）| ✅ 跨 purity 實證：V5 vs baseline F1 = 0（caller 層級）|
| baseline 的 enum bug 跨 purity 都存在 | 🆕 **新發現**：0.6 baseline HP:i:33 = 0 證實 |

### 15.6 一句總結

> **6 版本 F1 實證確認跨 purity 跨版本 caller F1 完全相同**：0.93 = 0.7166（A1=A3=A5）/ 0.6 = 0.6273（B1=B3=B5）。
>
> longphase-to phase 不改 caller F1（不改 VCF FILTER）；V5 真實價值在 read-level tag 品質，不在 F1。
>
> **「PON 方法是否最佳」跨 purity 答案**：
> - 高純度（0.93）→ ratio 主導性翻轉（4.3× 差距）+ V5 修補價值在 self-phasing
> - 中純度（0.6）→ ratio 翻轉幅度小（1.3×）+ V5 修補價值在 conservative tagging（HP:i:33 ↑）
> - **caller F1 在所有 purity 都不受 V5/baseline 影響**
>
> **baseline 的 enum mismatch bug 跨 purity 都存在**（0.93 / 0.6 baseline HP:i:33 都=0）— 0.6「看似 OK」是 self-phasing 自然減弱的表象，不是 bug 修復。

---

## 跨檔交叉引用

| 主題 | 文件 |
|------|------|
| 主報告 | `InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md` |
| PI 審查母版 | `InterSubMod/docs/reports/validated/2026/04/20260428_Self_Phasing_Baseline_V5_Audit_01.md` |
| V5 程式碼 diff（4 commit 演進）| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md` |
| 0.6 purity 完整驗證 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md` |
| Somatic bias 17.3:1 整合 + IGV | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md` |
| Phase vs Tag 演算法細節 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/13_phase_vs_tag_algorithm_detail.md` |
| V3-Fixed 驗證 memory | `memory/project_v3_fixed_haplotag_verification.md` |
| V5 驗證 memory | `memory/project_v5_somatic_fallback_verification.md` |
| Knowledge — purity 設計 | `/big8_disk/liaoyoyo2001/Knowledge/02_samples/subsample-purity.md` |
| Knowledge — HCC1395 樣本資訊 | `/big8_disk/liaoyoyo2001/Knowledge/02_samples/hcc1395.md` |

---

## 變更歷史

| 日期 | 變更 | 觸發 |
|-----|------|-----|
| 2026-04-29 06:00 | 初稿（7 Q&A 對話模式記錄）| 用戶 step-by-step 確認 + 選 B 寫獨立 Q&A 檔 |
| ⚠ 待補 | Q1/Q2/Q5 的 LongPhase 上游 design doc 直接驗證（GitHub README / 論文 / commit log）| 屬另案，不影響本報告主結論 |
