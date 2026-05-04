<!--
build_date: 2026-05-05
agent: fast-learning-coach skill (深度學習互動 → provenance audit)
status: validated
report_class: PI-audit / data-provenance
audience: PI / 自己未來 lab meeting / 論文撰寫者
parent_audit: InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md
inputs:
  - /big7_disk/liaoyoyo2001/longphase-to-mod/.git (commit history + reflog)
  - /big7_disk/liaoyoyo2001/longphase-to-mod/PhasingProcess.cpp:154-228
  - /big7_disk/liaoyoyo2001/longphase-to-mod/PhasingGraph.cpp:1854-1932
  - /big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp:484-563
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v2b/phasing.log
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/threshold_compare/{baseline_09,v5_flag}/run.log
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/v5_flag_force_path2only/run.log
  - InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md
  - InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md
outputs:
  - InterSubMod/docs/reports/validated/2026/05/20260505_self_phasing_V5_data_provenance_audit_01.md
related_memory:
  - memory/project_v5_somatic_fallback_verification.md (caveat 更新)
  - memory/project_self_phasing_causal_chain_confirmed.md
  - memory/project_pon_only_phasing_verification.md
verdict: CRITICAL FINDING — PI 報告 V5 數據為 Pass 1 only；4-30/5-01 已產出 Pass 2 觸發數據但未做對比分析；ploidy bug 已修但結論待重驗
last_verified: 2026-05-05
report_template: structured-tech-report v1.0 (provenance audit 變體)
-->

# Self-Phasing V5 Data Provenance Audit — 2026-05-05

## 0. TL;DR

> **PI 報告（2026-04-29）的全部 V5 數值（sanity 15/15 PASS、+13.3pp clean PS、AMB 17.5→8.0%、HP:i:33 −54%）都是「PON-only Pass 1 + tag layer fix」的結果。`d0bcd8c` (2026-04-30) 修了一個讓 ploidy 估值崩成 0 的關鍵 bug，使 Pass 2 second-round phasing 從未真正觸發。修補後（4-30 16:22 起）已重跑 baseline_09 / v5_flag / force_path2only 三組數據，purity 估值正常為 0.97-0.98，Pass 2 真實觸發 — 但這三組新 BAM 尚未做 sanity / concordance / ISM 對比分析。在新對比完成前，PI 報告 V5 結論需暫停為「Pass 1 only 條件下之觀察」，不能外推到完整 V5 (Pass 1 + Pass 2 + tag fix) 設計。**

四個必須立刻記住的事實：

| # | 事實 | 證據 |
|---|------|-----|
| 1 | PI 報告 V5 數據 = Pass 1 only | `output/pononly_v2b/phasing.log` 顯示 `purity: 0` + 無 "second round phasing" 字串 |
| 2 | V5 commit 鏈是 5 commits 不是 4 | `git log` 顯示 `8b8c1fd → 41ff147 → 380e8d2 → d0bcd8c → 938f0df` |
| 3 | working tree 未 commit caveat (R1) 已解決 | `git diff --stat HEAD` 為空（含 4-30 cherry-pick 完整化） |
| 4 | 4-30/5-01 新 BAM 尚未做對比分析 | `output/threshold_compare/`、`output/v5_flag_force_path2only/` 內無 .json/.tsv/.md 產出 |

---

## 1. 證據鏈（5 表 + 一個機制圖）

### 1.1 V5 commit 鏈最新狀態（5 commits，非 PI 報告所述 4 commits）

| 順序 | Git ref | Author date | **Commit date** | 修補在哪層 | 修什麼 |
|-----|---------|-----|-----|---------|-------|
| 1 | `8b8c1fd` | 2026-04-10 | 2026-04-10 | Phasing | 加 `--pon-only-phasing` flag → `convertNonGermlineToSomatic()` 降權 non-PON 變異 |
| 2 | `41ff147` | 2026-04-10 | 2026-04-10 | Tag | `getVote()` 重寫為兩層投票（germline first / somatic second）+ enum→int literal bug |
| 3 | `380e8d2` | 2026-04-10 | 2026-04-10 | Tag | `countINDELHaplotype()` 加 HAPLOTYPE_UNDEFINED guard |
| 4 | `d0bcd8c` | **2026-04-30** | **2026-04-30** | Phasing | **NEW**：fix ploidy collection timing — Pass 2 syncOrigins 後補回 ploidyRatio |
| 5 | `938f0df` | 2025-10-02 (上游 zhenyu 寫) | **2026-04-30** (本 fork cherry-pick) | Phasing | **NEW**：highPurity threshold `0.95` → `0.9` |

**關鍵釐清**：commit 4 與 commit 5 皆於 2026-04-30 15:27 本地 commit/cherry-pick 進入此 fork（`git reflog` 確認）。PI 報告 (4-29) 寫時，這兩個 commit 尚未存在；報告 caveat R1「working tree 未 commit 兩塊」實際對應的就是這兩個修補的雛形。

### 1.2 BAM/VCF 產出時間軸 vs commit 時間軸

```
2026-04-03 21:46  V2b phasing 跑 → output/pononly_v2b/tumor_phased.vcf
                  (binary 含 8b8c1fd；無 d0bcd8c → ploidy bug 觸發)
                  ↓
2026-04-10        commits 8b8c1fd, 41ff147, 380e8d2 完成（V2b/V3F/INDEL guard）
                  ↓
2026-04-12 13:43  V5 haplotag 跑 → output/pononly_v5_somatic_fallback/tumor_tagged.bam
                  (吃 V2b VCF，所以仍含 ploidy bug 後果 — Pass 2 從未觸發)
                  ↓
2026-04-29        PI 報告寫 (引用 pononly_v5_somatic_fallback BAM)
                  ↓  ⚠ 此時點所有 V5 結論 = Pass 1 only
                  ↓
2026-04-29 01:53  v5_93_purity_fix 跑（working tree 階段，purity=0.871，閾值 0.95→不觸發）
                  ↓
2026-04-30 04:09  v3f_no_pononly 跑（purity=0.927，觸發 Pass 2）
2026-04-30 12:00  v3f_ablation_ism_06 跑
                  ↓
2026-04-30 15:27  d0bcd8c commit + 938f0df cherry-pick + rebase abort 收尾
                  ↓
2026-04-30 16:22  threshold_compare/{baseline_09, v5_flag} phasing 跑
                  (purity = 0.977 / 0.984，**Pass 2 真實觸發**)
                  ↓
2026-04-30 17:10  threshold_compare 兩組 haplotag 完成
                  ↓
2026-05-01 03:37  binary 重編譯（穩定產物）
                  ↓
2026-05-01 11:46  v5_flag_force_path2only haplotag 完成
                  ↓
2026-05-05        本文件審核（尚無新 BAM 對比分析）
```

### 1.3 phasing.log purity=0 smoking gun

`/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v2b/phasing.log` 完整片段：

```
parsing BAM, somatic calling, and phasing
(chrY,10s)(chr21,344s)(chr13,428s)(chr19,443s)(chr22,556s)(chr17,579s)
... 24 chromosomes 共 832s 完成 Pass 1
purity: 0                          ← ploidy bug：ploidyRatioMap 為空 → q1=q3=0 → polynomial → clamp 到 0
export phasing result              ← 緊接 export，無 "second round phasing,"
(chrY,0s)(chr21,2s)(chr22,3s)...   ← 24 chr export 各 0-11 秒
parsing total: 844s
```

**機制驗證**（PhasingProcess.cpp:198-199）：
```cpp
bool highPurity = purity > 0.9;
if(highPurity){
    std::cerr << "second round phasing, ";   // ← 沒印出 = highPurity=false
}
std::cerr << "export phasing result" << std::endl;
```

→ 證據鏈 4 路全通：log 字串 / 時間特徵 / 程式碼邏輯 / 預期 cell line purity 矛盾，**HCC1395 5kHz V2b BAM 的 Pass 2 確認從未觸發**。

### 1.4 4-30/5-01 重跑後的 purity 結果（ploidy bug 已修證據）

| 實驗目錄 | mtime | purity | Pass 2 觸發 | log 來源 |
|---------|------|--------|-------------|----------|
| `pononly_v2b` | 4-03 | **0.0** ❌ | NO | `output/pononly_v2b/phasing.log` |
| `pononly_v5_somatic_fallback` | 4-12 | (吃 v2b VCF) | NO（繼承）| 同上 |
| `v5_93_purity_fix` | 4-29 01:53 | 0.871 | NO（< 0.95 當時閾值）| `output/v5_93_purity_fix/run.log` |
| `v3f_no_pononly` | 4-30 04:09 | 0.927 | **YES** | `output/v3f_no_pononly/phase.log` |
| `threshold_compare/baseline_09` | 4-30 16:22 | **0.977** | **YES** | `output/threshold_compare/baseline_09/run.log` |
| `threshold_compare/v5_flag` | 4-30 16:22 | **0.984** | **YES** | `output/threshold_compare/v5_flag/run.log` |
| `v5_flag_force_path2only` | 5-01 11:13 | 0.984 | (forced path) | `output/v5_flag_force_path2only/run.log` |

→ HCC1395 5kHz cell line 真實 purity ≈ 0.98（cell line 預期 >0.95），與 d0bcd8c 修補前的 0.0 估值形成**鐵證對照**。

### 1.5 d0bcd8c → cherry-pick → binary recompile 的完整時序

```
2026-04-30 15:27:07  d0bcd8c   commit (InterSubMod Research)
                                "fix(purity): collect ploidyRatio after PON-only Pass 2 syncOrigins"
                                ↓
2026-04-30 15:27:07  rebase aborted  ← reflog 顯示 "rebase (abort): updating HEAD"
                                ↓
2026-04-30 15:27:44  938f0df   cherry-pick (zhenyu upstream, author 2025-10-02)
                                "Update purity calculation in PurityCalculator
                                 to improve accuracy and adjust high purity threshold."
                                threshold: 0.95 → 0.9
                                ↓
2026-04-30 16:22     threshold_compare phasing 開始（用 working tree binary）
                                ↓
2026-05-01 03:37     重新編譯 binary（穩定版）
                                ↓
2026-05-01 11:13     v5_flag_force_path2only 用穩定 binary 跑
```

---

## 2. 對 PI 報告（2026-04-29）的 6 處修訂建議

| # | 原文段落 | 應改為 | 證據 |
|---|---------|--------|-----|
| 1 | §1 「修補在 longphase-to-mod fork 內以 4 commit 漸進完成」 | 「以 **5 commits** 漸進完成（4-10 完成 V2b/V3F/INDEL guard 三 commits；4-30 補 d0bcd8c ploidy fix + 938f0df threshold cherry-pick）」 | git log + reflog |
| 2 | §1 「V5 working tree 未 commit、Confidence threshold 0.6 未直接驗證」 | 「**working tree 已 commit (R1 解決)**；4-30/5-01 重跑數據已產出但對比分析未完成」 | `git diff --stat HEAD` 空；`output/threshold_compare/` 無 .json/.tsv |
| 3 | §6.1 ISM aggregate (HCC1395 5kHz) `version_summary.tsv` 引用 | 加註「以下 v5_new 數據為 4-12 BAM = Pass 1 only。完整 V5 (Pass 1 + Pass 2) 重跑後待回填」 | `pononly_v5_somatic_fallback` mtime 4-12；無 second round 字串 |
| 4 | §6.4 全基因組 paired ground-truth concordance「+13.3 pp / +8.3 pp」 | 加註「Pass 1 only 條件下測得；Pass 2 完整觸發後是否強化或持平待 4-30 BAM benchmark」 | 同 #3 |
| 5 | §6.5 「V5 修正後 HP1:HP2 ~1:1」「AMB% 17.5→8.0%」「HP:i:33 −54%」 | 加註「Pass 1 only 結論；Pass 2 觸發是否進一步改善需 threshold_compare/v5_flag BAM 重做 sanity check」 | 同 #3 |
| 6 | §8 caveat R1 「V5 working tree 未 commit」 | **刪除 R1**；新增 R11「Pass 2 verification gap — 4-30/5-01 已產 BAM 但 ISM benchmark / sanity 對比尚未完成；Pass 1-only V5 結論不可外推到完整 V5」 | git status 確認 + 4-30 後 output/ 內無對比檔 |

> **不直接 patch PI 報告**，PI 報告保留為 2026-04-29 的歷史快照（fact-check at that date）；本文件作為 follow-up + audit。

---

## 3. 三層因果鏈精準化（這次自學新理解，給 lab meeting 用）

### 3.1 phasing layer / tag layer / 修補兩層

```
[Phasing layer] PhasingGraph.cpp + PhasingProcess.cpp
  └─ self-phasing 根源：somatic 進 graph 反客為主
  └─ V2b 8b8c1fd 修：加 --pon-only-phasing flag (Pass 1 降權)
  └─ d0bcd8c 修：Pass 2 ploidy collection timing (Pass 2 真實觸發前提)
  └─ 938f0df 修：threshold 0.95 → 0.9 (Pass 2 觸發門檻)
  └─ 輸出：phased VCF (GT/PS/GT2) + LOH.bed
       │
       ▼
[Tag layer] HaplotagProcess.cpp::getVote()
  └─ self-phasing 放大：baseline getVote priority bug + enum mismatch
  └─ 41ff147 修：getVote 兩層投票（germline first / somatic second）
  └─ 380e8d2 修：countINDELHaplotype UNDEFINED guard
  └─ V5 working tree (已被 d0bcd8c 收尾)：Layer 1.5 somatic fallback
  └─ 輸出：tumor_tagged.bam (HP:i:1/2/11/21/33)
       │
       ▼
[ISM layer] InterSubMod ReadParser.cpp:120
  └─ HP:i tag → ISM 特徵 (HP_Ratio, HPFineNGroups, ...)
  └─ 本 repo 在此修補週期無 C++ 改動
```

**修補哲學**：必須兩層都動。
- 只修 phasing 不修 tag → tag 層 priority bug 仍會把對的 phase 解析錯
- 只修 tag 不修 phasing → tag 讀的是被 self-phasing 污染的 VCF，garbage in garbage out

### 3.2 共現 100% vs 50% 嚴格化（精通級才需要）

之前粗略說「somatic 共現 100% / germline 50%」 — 嚴格說法：

| 元素 | 在同 phase block 內 | 跨 phase block |
|-----|------------------|---------------|
| germline het | 100% 共現 | 50% 隨機（block 邊界打斷） |
| 同 sub-clone somatic | 100% 共現 | **100% 共現**（不受 block 邊界限制） |

→ 真正的不對稱：**somatic 跨整個 sub-clone 全基因組共現；germline 僅在 block 內共現** → phasing graph 看到 somatic-somatic edges 形成跨 block 的強連結 cluster，演算法把它當 phase block 骨幹。

### 3.3 r=0.001 攻防意義

n=288K pairs 下 r=0.001 比 17.3:1 偏移更具殺傷力：

| 指標 | 攻防意義 |
|------|---------|
| 17.3:1 偏移 | 可被反駁：「也許是 amplification bias / reference allele bias」 |
| **r = 0.001** | 任何 systematic bias 都會在 paired 也出現，不會破壞 r ≠ 0；**完全 random = TO 給的 HP_Ratio 跟正確答案無共同訊號** |

統計顯著性：n=288K 時 r=0.001 的 95% CI ≈ [-0.003, 0.005]，**統計上不顯著（p ≈ 0.6）**，正符合「完全無關」。
排除 paired 自身誤差：paired 內部一致性 r > 0.95（同 sample 不同 chr split 驗證）。

---

## 4. 「不存在的設計」修正 — Purity-aware two-pass phasing 真實存在

### 4.1 兩段式 phasing 結構（PhasingProcess.cpp:154-228）

```
┌─ Pass 1 (per-chromosome 平行) ──────────────────────────┐
│  if (params.ponOnlyPhasing) {  // V2b 起 = true        │
│      vGraph->convertNonGermlineToSomatic();            │
│      vGraph->phasingProcess(chrInfo.posPhasingResult,  │
│                              chrInfo.LOHSegments,      │
│                              nullptr);                 │ ← Pass 1 不收 ploidy
│      vGraph->resetNonPonOrigin();                      │
│      vGraph->somaticCalling(...);                      │
│      vGraph->syncPhasingResultOrigins(...);            │
│      vGraph->collectPloidyRatio(&ploidyRatioMap);      │ ← d0bcd8c NEW
│  } else {                                              │
│      vGraph->somaticCalling(...);                      │
│      vGraph->phasingProcess(...);                      │
│  }                                                     │
└────────────────────────────────────────────────────────┘
                       ↓
┌─ purity 計算 ───────────────────────────────────────────┐
│  double purity = PurityCalculator::getPurity(...);     │
│  bool highPurity = purity > 0.9;  // 938f0df: 0.95→0.9 │
└────────────────────────────────────────────────────────┘
                       ↓
┌─ Pass 2 (per-chromosome 平行，only if highPurity) ──────┐
│  if (highPurity) {                                     │
│      vGraph->convertNonGermlineToSomatic();            │
│      chrInfo.posPhasingResult = PosPhasingResult();    │
│      vGraph->phasingProcess(..., nullptr);             │
│  }                                                     │
│  vGraph->exportPhasingResult(...);                     │
└────────────────────────────────────────────────────────┘
```

### 4.2 PurityCalculator 多項式（PhasingGraph.cpp:1922-1928）

`purity = polynomial(q1, q3, lohRatio)` —— 用 SNP allele freq 分布的 quartile（q1, q3）+ LOH ratio 回歸出純度估值。**不是用戶設定，是演算法自己估**。多項式有兩組（依 caller 切換），結尾 clamp 到 [0, 1]。

### 4.3 highPurity threshold 0.95 → 0.9 變化

- **2025-10-02** 上游 zhenyu 在 longphase-to GitHub 主線把 threshold 從 `0.95` 改 `0.9`（commit 938f0df author date）
- **2026-04-30 15:27** 本 fork cherry-pick 進來
- 影響：之前 purity 介於 (0.9, 0.95) 的樣本 **不會** 觸發 Pass 2；現在會觸發

### 4.4 d0bcd8c ploidy bug 與影響

**Root cause**（commit message 摘要）：
- Pass 1 的 `phasingProcess(...)` 第三參數傳 `nullptr` → 跳過 ploidy 收集
- Pass 2 的 `somaticCalling + syncPhasingResultOrigins` 不會自己補回 `ploidyRatioMap`
- 結果：`mergedPloidyRatioMap` 為空 → `q1 = q3 = 0` → polynomial 算出負數 → clamp 到 0
- → `highPurity = (0 > 0.9) = false` → Pass 2 永遠不觸發

**影響範圍**：
- 4-30 之前產出的所有 V5 BAM（含 PI 報告引用的 `pononly_v5_somatic_fallback`）= Pass 1 only
- 4-30 16:22 之後的數據（threshold_compare）= Pass 2 真實觸發

---

## 5. 「V5 是否最佳」答案重述（依新事實）

| 路徑 | 機制 | 解決 self-phasing? | 在本專案的角色 |
|------|-----|-------------------|--------------|
| A: PON-only + tag fallback | 現行 V5 (5 commits) | ✅ Pass 1 already 解；Pass 2 對 high-purity 樣本進一步收斂 | **production**，但 Pass 2 效益待 4-30 數據對比驗證 |
| B: iterative phasing | 兩輪 phasing 中間過濾 somatic | 理論可行但工程量大 | 未實作 |
| C: long-read native phasing | 直接靠 long reads 物理連結 haplotype | ✅ 從根本繞過 variant graph | 需 ≥20kb HiFi/ONT；未在本專案路徑 |
| F: paired sequencing | matched normal 提供 germline scaffold | ✅ 完全消除 | 黃金標準但資料約束 |

**結論**：在 tumor-only + PON-only + ONT short-medium read 的設定下，V5 (5-commits) 仍是合理的修補路徑。**但 Pass 2 真實效益尚未獨立驗證** — 直到 4-30/5-01 數據完成 ISM benchmark 對比，「完整 V5 = Pass 1 + Pass 2 + tag fix」是否優於「Pass 1 + tag fix only」是 open question。

---

## 5.5 4-30/5-01 新數據盤點（**待補對比分析**）

### 5.5.1 6 個新實驗目錄

| 目錄 | mtime | purity | 用途 | 對比目標 |
|------|------|--------|-----|---------|
| `output/v5_93_purity_fix/` | 4-29 01:53 | 0.871 | working tree 階段驗證 | （閾值未過，比較不出 Pass 2）|
| `output/v3f_no_pononly/` | 4-30 04:09 | 0.927 | 不開 PON-only 也跑 | 對照 PON-only 是否必要 |
| `output/v3f_ablation_ism_06/` | 4-30 12:00 | TBD | 0.6 purity simulation | PI 報告 R7 caveat |
| `output/threshold_compare/baseline_09/` | 4-30 16:22 | 0.977 | baseline + threshold 0.9 | **vs v5_flag** |
| `output/threshold_compare/v5_flag/` | 4-30 16:22 | 0.984 | V5 flag + threshold 0.9 | **vs baseline_09** |
| `output/v5_flag_force_path2only/` | 5-01 11:13 | 0.984 | 強制 Pass 2 only | 對照「Pass 1+2」vs「Pass 2 only」 |

### 5.5.2 Pass 2 觸發的真實樣貌

`output/threshold_compare/v5_flag/run.log` 顯示：
```
purity: 0.983791
second round phasing, export phasing result
```

第二輪 phasing 真實觸發 — 與 `pononly_v2b/phasing.log` 的 `purity: 0` + 無 second round 形成鐵證對照。

### 5.5.3 待補對比分析（這是 PI 報告結論能否更新的關鍵）

**目前缺**（`find output/threshold_compare/ -name "*.json" -o -name "*.tsv" -o -name "*.md"` 為空）：

| 對比 | 預期意義 | 對 PI 報告影響 |
|------|---------|---------------|
| `threshold_compare/v5_flag/` BAM vs `pononly_v5_somatic_fallback/` BAM (4-12) | 量化 Pass 2 對 sanity check 15 點的貢獻 | 「+13.3pp clean PS」是 Pass 1 only / Pass 1+2 哪一版？ |
| `threshold_compare/v5_flag/` vs `threshold_compare/baseline_09/` | threshold 0.9 + Pass 2 觸發下 V5 flag 是否仍贏 baseline | 「V5 顯著優於 baseline」是否仍站得住 |
| `v5_flag_force_path2only/` vs `threshold_compare/v5_flag/` | Pass 1+2 vs Pass 2 only | Pass 1 (PON-only) 對最終結果是必要還是冗餘？ |

**潛在三種情境**：
- **情境 A**：Pass 2 後 V5 vs baseline 差距更大 → PI 報告結論強化
- **情境 B**：Pass 2 後差距不變 → V5 tag fix 才是主因，Pass 2 只在 phased VCF 層改善（不直接影響 ISM HP_Ratio）
- **情境 C**：Pass 2 後差距縮小 → V5 真實效益小於 PI 報告宣稱（需要修訂結論）

---

## 6. 後續行動

| ID | 動作 | 優先 | 對應 PI caveat |
|----|------|-----|--------------|
| **T1** (NEW) | 4-30 新 BAM 跑 ISM benchmark + sanity check 15 點對齊 PI 報告 §6.4 / §6.5 數字 | **P0 高** | PI R11 (新增) |
| **T2** | trace `pononly_v5_somatic_fallback/` (4-12) 是否要重產（現在有 d0bcd8c 應重跑取代）| **P1 中** | PI §6 全部 |
| **T3** | 7 樣本 V5 BAM 全量重跑（HCC1395_DORADO / HCC1937 / HCC1954 / H1437 / H2009 / COLO829）| **P1 中** | PI R3 |
| **T4** | PI 報告依本文件 §2 修訂建議 patch（或保留 4-29 為歷史快照，本文件作為 errata）| **P2 中** | PI 全文 |
| **T5** | manifest.yaml 加 `haplotag_version` + `binary_commit_hash` 欄位 | **P2 中** | PI R3 / F4 |
| **T6** | 0.6 purity simulation 用 d0bcd8c 後 binary 重做（已有 v3f_ablation_ism_06，待對齊）| **P2 中** | PI R7 |

**T1 詳細子步驟**（建議下一個 cycle 執行）：
1. 用 `compare_haplotag_ism.py` 對 `threshold_compare/v5_flag/tumor_tagged.bam` 跑 ISM
2. 對 `threshold_compare/baseline_09/tumor_tagged.bam` 跑 ISM 作 baseline
3. 對 4-12 的 `pononly_v5_somatic_fallback/tumor_tagged.bam` 跑 ISM 作對照（驗證舊 PI 數字）
4. 比對：HP_Ratio median / Potential_LOH% / HP_Ratio AUC / sanity 15 點 / clean PS pp
5. 產出 `output/threshold_compare/comparison_summary.json` + `comparison_report.md`
6. 依結果決定是否更新 PI 報告 / 是否觸發 T3 7 樣本擴展

---

## 7. 引用文件 + commit hash 清單

### 7.1 程式碼路徑（本次審核引用，longphase-to-mod fork）

| 路徑 | 用途 |
|------|-----|
| `/big7_disk/liaoyoyo2001/longphase-to-mod/PhasingProcess.cpp:154-228` | Pass 1 / Pass 2 核心邏輯 |
| `/big7_disk/liaoyoyo2001/longphase-to-mod/PhasingGraph.cpp:1854-1932` | PurityCalculator 多項式實作 |
| `/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp:484-563` | getVote 兩層投票 + Layer 1.5 |
| `/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.h:66-68` | Tag 層介面契約（5 commits 零變動）|

### 7.2 V5 commit 鏈完整 hash

| Hash | Commit date | 修補 |
|------|------------|------|
| `8b8c1fd` | 2026-04-10 | feat: add --pon-only-phasing mode |
| `41ff147` | 2026-04-10 | fix(haplotag): two-layer getVote |
| `380e8d2` | 2026-04-10 | fix(haplotag): countINDELHaplotype guard |
| `d0bcd8c` | 2026-04-30 | fix(purity): collect ploidyRatio after Pass 2 syncOrigins |
| `938f0df` | 2026-04-30 (cherry-pick) | Update purity calculation + threshold 0.9 |

### 7.3 InterSubMod 文件參照

| 類別 | 文件 |
|------|------|
| 父審核（4-29 PI 報告）| `InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md` |
| Self-phasing 根因 | `InterSubMod/docs/reports/research_landscape/02_Self_Phasing根因.md` |
| 4-28 V5 audit baseline | `InterSubMod/docs/reports/validated/2026/04/20260428_Self_Phasing_Baseline_V5_Audit_01.md` |
| V5 audit suite 母索引 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/00_INDEX.md` |
| 4-24 V5 vs Baseline 全基因組 | `InterSubMod/docs/reports/pi_reports/2026/04/20260424_V5_vs_Baseline_complete_comparison_01.md` |
| ReadParser HP demotion | `InterSubMod/src/core/ReadParser.cpp:120` |

### 7.4 Knowledge Base

- `/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-to.md`
- `/big8_disk/liaoyoyo2001/knowledge/06_workflows/phasing-workflow.md`
- `/big8_disk/liaoyoyo2001/knowledge/03_file_formats/vcf-longphase.md`

---

## 附錄 A：本次審核發現的歷程

| 時間（session） | 發現 | 觸發 |
|---|------|------|
| Step 1 | 用戶質疑「V5 是不是有更好的做法」+ 我答「沒有 purity-aware branching」 | 被用戶反詰「為什麼是這樣」 |
| Step 2 | 用戶硬要核查 → 我 grep `purity` 找到 `PurityCalculator` + `bool highPurity = purity > 0.9` | 用戶精通級拒絕模型權威壓過 source code |
| Step 3 | 找 phasing.log → 看到 `purity: 0` smoking gun | 我自己延伸驗證 |
| Step 4 | git log 發現 V5 是 5 commits 不是 4，d0bcd8c + 938f0df 是 4-30 才完成 | reflog + commit date |
| Step 5 | 找 4-30 後產出 → 6 個新實驗目錄 + 3 個 purity > 0.9 顯示 Pass 2 真實觸發 | 主動 find recent files |
| Step 6 | 確認新 BAM 沒有對比 .json/.tsv 產出 | 用戶要求確認結論時的 follow-up |

**教訓**：
1. **PI 報告引用的數據必須註明 binary commit hash + BAM mtime**（避免時間漂移）
2. **「working tree 未 commit」caveat 不該長期掛著** — 隔天就該 commit 並 promote 為穩定版
3. **修了關鍵 bug（如 ploidy timing）後，所有依賴下游分析都應該被視為 stale 直到重跑**
4. **AI 給 PI 級結論時必須先核查 source code 與最新 git log**，不能靠記憶或舊報告

---

## 附錄 B：變更歷史

| 日期 | 變更 | 觸發 | 作者 |
|-----|------|-----|-----|
| 2026-05-05 | 初稿建立 — 5 commits 鏈確認、smoking gun 驗證、4-30 新數據盤點 | fast-learning-coach skill 學習過程中用戶的精通級追問逼出 audit | fast-learning-coach skill (互動式 audit) |
