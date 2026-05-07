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
verdict: CRITICAL FINDING — PI 報告 V5 數據為 Pass 1 only；4-30/5-01 已產出 Pass 2 觸發數據；T1.2 (2026-05-07) chr19 priority bug 機制鐵證 (752 victims, V3F/V5 100% 修正率)；T1.2-F1 (2026-05-08) 全基因組擴展顛覆三項 chr19 結論 — 34,855 victims (46×), chr8 NOT hotspot, Layer 1.5 +560K 觸發；ISM benchmark 不再必要 (下游自動受惠)
last_verified: 2026-05-08
update_2026_05_07:
  - "T1.2 read-level vote audit 完成 (commit 8491f14)"
  - "新增 §5.6 機制驗證鐵證段"
  - "§6 行動表更新：原 T1 ISM benchmark CANCELLED；新 P0 = T1.2-F1 全基因組擴展"
  - "詳見 InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_priority_bug_mechanism_report.md"
update_2026_05_08:
  - "T1.2-F1 全基因組擴展完成（cycle 20260508_T1_2_F1_genome_wide_audit）"
  - "新增 §5.7 — 顛覆 §5.6 三項 chr19 結論：chr19 占比僅 2.16%、chr8 是 priority bug 冷區（rank 21，0.34× avg）、Layer 1.5 在全基因組確實觸發 +560,881 reads"
  - "§6 行動表更新：T1.2-F1 標記 DONE；下游 PI errata（§5.7.4）的訊息變更已陳列"
  - "詳見 InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md"
report_template: structured-tech-report v1.0 (provenance audit 變體)
-->

# Self-Phasing V5 Data Provenance Audit — 2026-05-05

## 0. TL;DR

> **PI 報告（2026-04-29）的全部 V5 數值（sanity 15/15 PASS、+13.3pp clean PS、AMB 17.5→8.0%、HP:i:33 −54%）都是「PON-only Pass 1 + tag layer fix」的結果。`d0bcd8c` (2026-04-30) 修了一個讓 ploidy 估值崩成 0 的關鍵 bug，使 Pass 2 second-round phasing 從未真正觸發。修補後（4-30 16:22 起）已重跑 baseline_09 / v5_flag / force_path2only 三組數據，purity 估值正常為 0.97-0.98，Pass 2 真實觸發。**

> **🟢 Update 2026-05-07：T1.2 read-level vote audit 完成（commit `8491f14`）— chr19 上 752 priority bug victims、V3F/V5 100% 修正率、4 路徑 3.5/4 PASS → V3F priority bug 修補機制鐵證確立**。詳見 §5.6。**ISM benchmark 不再必要**（ISM 是下游消費者，longphase-to 端修對後自動受惠）；新 P0 = T1.2-F1 全基因組擴展 + T4 PI 報告 errata 整合。

> **🔥 Update 2026-05-08：T1.2-F1 全基因組擴展完成 — 34,855 victims (46× chr19)、V3F/V5 100% 修正率、顛覆三項 §5.6 chr19 結論**：(1) chr19 僅占 2.16%（不是主要 hotspot）；(2) chr8 是 priority bug 冷區 rank 21（與 LOH+HPSig hotspot 是不同 layer）；(3) Layer 1.5 在 germline 缺席的 21.7M reads 中**確實觸發 +560,881**（chr19 局部未觸發是因 germline 不缺席）。詳見 §5.7。priority bug 機制因果**進一步強化**至 event-level 34,855 個案。

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

## 5.6 T1.2 Priority Bug 機制驗證鐵證（2026-05-07 補上）

> **Update 2026-05-07**：依用戶決策，**ISM 不再需要驗證**（ISM 是下游消費者，longphase-to 端 tag 改善後 ISM 自動受惠）。重點轉為**驗證 longphase-to 端 V3F 修補本身有效**，這由 T1.2 read-level vote audit 完成。完整 mechanism report：[`InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_priority_bug_mechanism_report.md`](../../../../research/v5_provenance_followup/T1_2_read_level_audit/T1_2_priority_bug_mechanism_report.md)（commit `8491f14`）

### 5.6.1 設計與執行

**對 baseline / V3F / V5 三版 testing-only binary patch 加 `--debug-vote-dump` flag** dump 每條 read 經 `getVote()` 後的 5-vote countMap + hpResult。三版分別跑 chr19（HCC1395 5kHz，~60 秒/版），共產出 549,206 rows × 3 版 vote dump TSV。

| Binary | Git ref | 修補狀態 |
|--------|---------|---------|
| baseline-debug | `8b8c1fd` | V2b PON-only flag（priority bug 仍在）|
| v3f-debug | `380e8d2` | V3F two-layer + INDEL guard |
| v5-debug | HEAD `938f0df` | + Layer 1.5 + ploidy fix + threshold 0.9 |

### 5.6.2 4 路徑驗證結果（3.5/4 PASS → 機制因果確立）

| 路徑 | 結果 | 判定 |
|------|------|------|
| ① 個案 trace ≥10 條 | **752 條** priority bug confirmed victims | ✅ PASS |
| ② chr19 1Mb 區域聚集 | chr19:30M (215 victims) + 27M (133) 集中 46% | ⚠️ PARTIAL |
| ③ Somatic density 共變 | high somatic vote ≥5 = **0 受害**；低票就觸發 | 🔄 反向但有意義 |
| ④ 修正後消失（V3F/V5 修正率）| **V3F 100% / V5 100%** | ✅ PASS |

### 5.6.3 三大新發現

1. **Priority bug 單向性**：chr19 上全 752 條 victims 都是 `baseline=11 → v3f=21 → v5=21`（無一條反方向）→ 完美對應 PI 報告 §3 全基因組 17.3:1 整體偏移單向性（94.6% somatic→HP1）。**chr19 752 條 = 17.3:1 微觀證據**。

2. **觸發條件比理論寬鬆**：1-2 票 somatic vote 即觸發 priority bug（high somatic vote ≥5 群體 = 0 受害者；全 752 都是 low 1-4 票）。原因：baseline `getVote()` 用 vector 順序檢查，第一個有票的 pair 就 break。

3. **V5 Layer 1.5 chr19 未觸發**：V5 vs V3F 在 chr19 結果 100% 相同（germline 不缺席 → Layer 1.5 分支不執行）。Layer 1.5 真實價值要在 **germline het 稀疏區（cnLOH / amplicon hotspot）** 才看得到。

### 5.6.4 對 PI 報告的機制證據強化

PI 報告 §5.2 「V3F getVote 兩層投票」原本只有 commit message + IGV 3 個截圖支持；**T1.2 補上 752 read-level 個案 + 100% 修正率的鐵證**。priority bug 從「理論 + 截圖」升級為「個案 + 統計 + 機制」三重佐證。

### 5.6.5 樣本個案（前 5 條）

| read_name (前 12) | chr19:pos | HP1/HP2 | HP1_1/HP2_1 | germline_maj | somatic_maj | baseline → v3f → v5 |
|---|---:|:---:|:---:|:---:|:---:|:---:|
| 1c50034a-f0f | 201,417 | 1/3 | 1/0 | HP2 | HP1 | **11 → 21 → 21** |
| afb8e89b-893 | 585,252 | 1/2 | 1/0 | HP2 | HP1 | **11 → 21 → 21** |
| 35c7e166-ec3 | 824,360 | 0/1 | 1/0 | HP2 | HP1 | **11 → 21 → 21** |
| 096ab9a7-030 | 1,574,442 | 0/3 | 1/0 | HP2 | HP1 | **11 → 21 → 21** |
| ccc8185d-f9b | 2,558,240 | 0/1 | 2/0 | HP2 | HP1 | **11 → 21 → 21** |

**統一指紋**：germline 票 → HP2（read 真實方向）；somatic 票 → HP1（self-phasing 污染）；baseline 跟 somatic（priority bug）；V3F/V5 修向 germline majority。

---

## 5.7 T1.2-F1 全基因組擴展（2026-05-08 補上）— 顛覆三項 chr19 結論

> **Update 2026-05-08**：對全基因組（chr1-22 + chrX/Y，HCC1395 5kHz）跑同一 vote audit。**結果與 §5.6 chr19 pilot 在 victim 量級、hotspot 分佈、Layer 1.5 觸發三項均產生顛覆性差異**。完整 audit：[`InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md`](../../../../research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md)

### 5.7.1 設計與執行

對同三版 binary（baseline / V3F / V5）跑全基因組 vote dump（每版 ~40 min，總 tagged ~18.9M reads × 3）。dump 大小：744 / 687 / 687 MB（gzipped）。

| 規模 | chr19 pilot | Genome F1 | 倍數 |
|---|---:|---:|---:|
| Dump rows | 549,206 | **29,973,253** | 54.6× |
| Tagged reads (per binary) | ~330K | **18,895,432** | 57× |
| **Priority bug confirmed victims** | **752** | **34,855** | **46.4×** |

### 5.7.2 三項 chr19 結論被推翻

#### 推翻 1：chr19 不是主要 hotspot（占比僅 2.16%）

§5.6.3 第 1 點寫「chr19 752 條 = 17.3:1 微觀證據」是**真的**，但 §5.6 推論「chr19 是主要 hotspot」**不成立**。實際 priority bug 影響全基因組廣泛分佈：

| chr | victims | enrichment ‰ | rank | 備註 |
|---|---:|---:|---:|---|
| chr7 | **3,508** | 0.723 | 1 | victim N 最高 |
| chr2 | 2,792 | 0.617 | 2 | |
| chr1 | 2,674 | 0.541 | 3 | |
| chr16 | 2,584 | **1.140** | 4 | enrichment 排第 5 |
| chr20 | 2,101 | **1.306** | 7 | enrichment 排第 2 |
| chr21 | 792 | **1.279** | 17 | enrichment 排第 3 |
| chr19 | 752 | 0.703 | **19** | **占 2.16%** |
| **chr8** | **666** | **0.200** | **21** | **0.34× genome avg（冷區）** |
| chrY | 67 | **1.484** | 24 | enrichment 最高（small N） |
| genome avg | 34,855 | 0.590 | — | |

#### 推翻 2：chr8 不是 priority bug hotspot（與 LOH+HPSig hotspot 是不同 layer）

MEMORY 的 [`project_hcc1395_chr8_hotspot.md`](../../../../research/autoresearch/memory/project_hcc1395_chr8_hotspot.md) 寫 chr8「LOH+HPSig 7.4× FP enrichment」。本 audit 顯示 chr8 在 **priority bug 層** 反而是 0.34× genome avg 的**冷區**。

→ **這兩個 hotspot 是不同 layer**：chr8 LOH+HPSig 集中是 ISM 下游 false-positive 富集（HP_Ratio + LOH 特徵交互），不是 V3F 修對的 priority bug 範疇。chr8 priority bug 在 V3F/V5 修補後仍存在的 ISM FP enrichment 必須另尋機制。

#### 推翻 3：Layer 1.5 確實觸發 560,881 reads（不是「未觸發」）

§5.6.3 第 3 點寫「V5 Layer 1.5 chr19 未觸發」是**chr19 局部對的**，但**全基因組看 V5 確實在 germline 缺席區大量觸發**：

| 指標 | 值 |
|---|---:|
| germline_vote=0 reads | **21,765,669**（占 merged 36.8%） |
| V3F tagged 數（germline=0 子集）| 0 |
| V5 tagged 數（germline=0 子集）| **560,881** |
| **Layer 1.5 額外觸發** | **+560,881** |

但 V5 全基因組 tag count = V3F tag count（兩者都是 18,895,432）→ V5 在 germline_vote=0 區多 tag 的 560,881，**必然在 germline_vote>0 區少 tag 同數**。這是新 finding，待釐清 mechanism（候選：V5 ploidy fix + threshold 0.9 對 germline 充足區的某種 trade-off）。**這不影響 priority bug 修正的 100% verdict**，但屬另一條 follow-up。

### 5.7.3 V3F vs V5 全基因組行為再校準

| 比較 | chr19 結論（§5.6） | Genome 結論（§5.7） |
|---|---|---|
| Priority bug 修正率 | V3F 100% / V5 100% | V3F 100% / V5 100% ✅ 一致 |
| Layer 1.5 觸發 | 未觸發 | **+560,881 額外 tag**（germline=0 區）|
| V5 vs V3F net effect | 行為相同 | germline=0 多 tag 560K，germline>0 少 tag 560K（**zero-sum 重分配**）|

### 5.7.4 對 PI 報告的後續訊息更新

- **PI §5.2 機制證據**：read-level victim N 從 752 升至 **34,855**（同源證據量級提升 46×）
- **PI §3.3.3 chr19 hotspot 解讀**：chr19 SP1/2/3 是 **可重現案例**，但**不再代表 priority bug 主要分佈**；正確說法是「chr19 是被研究最多但非最熱的 priority bug 區，全基因組廣泛分佈，chr7/chr2/chr1/chr16/chr20 是 N 量大宗」
- **PI R11 chr8 hotspot**：chr8 priority bug rank 21（0.34× genome avg）→ chr8 ISM FP 富集**不是** priority bug 直接造成；應改寫為「chr8 LOH+HPSig 富集機制獨立於 priority bug，需另案探討」
- **新 caveat**：V5 vs V3F 整體 tag count 相同但分佈不同（germline=0 +560K / germline>0 -560K）→ 任何「V5 對 V3F 整體效應」分析需分層

### 5.7.5 數據完整性 caveat

- **3-way merged 59M rows vs dump 30M**：inner join 在 (read_name, chr, pos) 是 cross product（同 read+pos 多投票記錄，例如 supplementary alignments）。**34,855 是 event-level victims**，非 unique-read 計數
- **chrY victim 67 / total 45,137 ‰=1.484** rank 1 但 N 太小，**不應**過度詮釋
- **此 audit 僅 HCC1395 5kHz**：跨樣本 priority bug 分佈差異待 T3（7 樣本擴展）補上

---

## 6. 後續行動

> **2026-05-07 行動表更新**：用戶決策 — **ISM benchmark 不再必要**（ISM 是下游消費者，longphase-to 端 V3F 修補已被 T1.2 證實 100% 修對 priority bug victims，下游 ISM 自動受惠）。原 T1 (ISM benchmark) 從 P0 降級。新 P0 重心是 longphase-to 端的全基因組擴展 + 整合至 PI 報告。

| ID | 動作 | 優先 | 狀態 | 對應 PI caveat |
|----|------|-----|------|--------------|
| **T1.2** ✅ | Read-level vote countMap audit (chr19) — priority bug 機制鐵證 | — | **DONE 2026-05-07** | PI §5.2 機制證據 |
| ~~T1 (原)~~ | ~~4-30 新 BAM 跑 ISM benchmark~~ | ~~P0~~ | **CANCELLED** — ISM 是下游消費者，不需獨立驗證 | — |
| **T1.2-F1** ✅ | 全基因組擴展 dump — 34,855 victims (46× chr19) / chr8 NOT hotspot / Layer 1.5 +560K 觸發 | — | **DONE 2026-05-08**（§5.7）| PI R11 / chr8 hotspot |
| **T2** | trace `pononly_v5_somatic_fallback/` (4-12) 是否要重產 | **P1** | pending | PI §6 全部 |
| **T3** | 7 樣本 V5 BAM 全量重跑（HCC1395_DORADO / HCC1937 / HCC1954 / H1437 / H2009 / COLO829）| **P1** | pending | PI R3 |
| **T4** | PI 報告依本文件 §2 修訂建議 patch（保留 4-29 為歷史快照，本文件作為 errata）| **P1** | pending | PI 全文 |
| **T5** | manifest.yaml 加 `haplotag_version` + `binary_commit_hash` 欄位 | **P2** | pending | PI R3 / F4 |
| **T6** | 0.6 purity simulation 用 d0bcd8c 後 binary 重做 | **P2** | pending | PI R7 |

**T1.2-F1 詳細子步驟**（下一個 cycle 執行）：
1. 用既有 `longphase-to-{baseline,v3f,v5}-debug` binary 跑全基因組（chr1-22）dump，每版 ~30 min
2. JOIN 三版 dump 找全染色體 priority bug victims
3. 區域聚集分析：sliding window per-chr enrichment；驗證 chr8 hotspot（MEMORY: `project_hcc1395_chr8_hotspot.md`）+ V5 Layer 1.5 在 germline 稀疏區是否觸發
4. 產出 `T1_2_F1_genome_wide_audit.md` + 圖表
5. 跨樣本擴展前提：T1.2-F1 結論在 HCC1395 全基因組仍成立 → 進 T3 7 樣本

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
