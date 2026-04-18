<!--
建立時間: 2026-04-03 22:00
目標: 記錄 LongPhase-TO self-phasing circular dependency 修正方案的完整推論、觀察、驗證計畫與方法學論證
處理範圍: LongPhase-TO phasing pipeline 修改設計、haplotag HP tag 影響分析、ISM 驗證框架、生物資訊方法學論證
關聯檔案:
  - /big7_disk/liaoyoyo2001/longphase-to-mod/PhasingProcess.cpp (核心修改)
  - /big7_disk/liaoyoyo2001/longphase-to-mod/PhasingGraph.cpp (graph 邏輯)
  - /big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp (haplotag 邏輯)
  - research/loh_investigation/reports/20260403_pon_only_phasing_verification_report.md
  - research/loh_investigation/reports/20260403_pon_only_haplotag_ism_verification_report.md
-->

# LongPhase-TO 二次 Phasing 方案：設計、推論與驗證計畫

## 摘要

本文件記錄 LongPhase-TO self-phasing circular dependency 的完整分析過程，包含：問題發現、根因確認、第一版修正（PON-only phasing）的優缺點、第二版修正（二次 phasing + 正確分類）的設計推論、以及用 sSNV benchmark 驗證的完整框架。

---

## 1. 問題發現與根因

### 1.1 觀察到的現象

ISM 在 TO 模式下的 HP_Ratio 出現系統性偏移：

| 指標 | Paired（金標準） | TO-Baseline | 差距 |
|------|---------------|-------------|------|
| HP_Ratio median (TP) | 0.5000 | 0.8358 | 偏向 HP1 |
| HP_Ratio extreme rate | 46.8% | 60.0% | +13.2pp |
| ISM-only LOH excess | 3.9% | 15.4% | +11.5pp |

### 1.2 根因確認：Self-Phasing Circular Dependency

**程式碼層面**（PhasingProcess.cpp）：

```
Line 153: vGraph->addEdge(readVariantVec)            ← 所有 variants 加入 graph
Line 156: vGraph->somaticCalling(...)                ← somatic 識別在 addEdge 之後
Line 161: vGraph->phasingProcess(...)                ← phasing 使用含 somatic 的 graph
```

**因果鏈**：
1. Somatic variants 在 `addEdge()` 時以 edges 加入 phasing graph
2. 因為 somatic mutation 是 subclonal，ALT allele 一致出現在同一組 reads
3. Phasing algorithm 將這些 reads 統一分配到 HP1
4. `somaticCalling()` 之後執行，但 graph edges 已被汙染
5. 結果：94.6% somatic variants → HP1，bias = 17.3:1

**論文是否承認？** 論文（碩士論文）未討論此 circular dependency。論文將 somatic 參與 phasing 視為 feature（建構 somatic haplotype branches），但未分析其對下游工具（如 ISM）的 side effect。

### 1.3 Self-Phasing 對 ISM 的影響路徑

```
Self-phasing bias (17.3:1)
  → Haplotag: somatic reads 偏向 HP1 → HP:i:1/HP:i:11 佔絕對多數
  → ISM: HP_Ratio 偏向 HP1（median 0.8358 vs paired 0.5000）
  → Potential_LOH 判定過度（60% vs 47%）
  → VerificationClass 退化（Strong/Weak 降為 Noise）
  → QS LOH penalty 反向作用 → TO mode 被迫停用 LOH penalty
  → HP 維度分析效力大幅降低
```

---

## 2. 第一版修正：PON-Only Phasing

### 2.1 設計

在 phasing 前呼叫 `convertNonGermlineToSomatic()`，將所有非 PON variants 標記為 SOMATIC，使 phasing graph 只由 PON germline 主導。

### 2.2 驗證結果

| 指標 | Baseline | PON-only | 判定 |
|------|----------|----------|------|
| LOH.bed Jaccard | — | **1.0000** | 完全不變 |
| Somatic bias | **17.3:1** | 消除 | Circular dependency 修正 |
| Phase block N50 | 4,061 | **8,109 (+99.7%)** | 品質翻倍 |
| Phased rate | 54.9% | **78.5%** | +23.6pp |
| SEQC2 LOH Jaccard | 0.8470 | **0.8470** | 不變 |
| 執行時間 | 2693s | **1976s (1.36x)** | 加速 |
| Purity | 0.927 | **0（失效）** | 副作用 |

### 2.3 Haplotag + ISM 驗證：發現新問題

| 指標 | Paired | TO-Baseline | TO-PON-Only |
|------|--------|-------------|-------------|
| HP_Ratio extreme rate (TP) | 46.8% | 60.0% | **99.9%** |
| Both HP > 0 rate | 67.1% | 55.9% | **0.15%** |
| HPMergedSig rate | 37.4% | 29.3% | **0.003%** |
| ISM-only LOH excess | 3.9% | 15.4% | **54.8%** |
| F1 (SuggestFilter) | 0.9282 | 0.9284 | **0.9291** |

**HP_Ratio 反而更差！** 原因：

1. 所有非 PON variants 被標為 SOMATIC → VCF 中 GT=`0|0`, GT2=`.|1`
2. Haplotag 解析 GT=`0|0` → `refHaplotype = UNDEFINED`
3. REF reads 無法投票（被忽略）
4. ALT reads 投票 → HP:i:21
5. getVote() 優先級：somatic > germline → 即使有 germline 投票也被覆蓋
6. 結果：somatic 位點所有 reads → 同一 HP → HP 二分群崩壞

### 2.4 第一版優缺點

| 優點 | 缺點 |
|------|------|
| Phasing 品質大幅改善（N50 翻倍） | HP tag 全部錯誤（HP:i:21 覆蓋所有 reads） |
| Somatic bias 消除 | HP 二分群完全崩壞（0.15%） |
| LOH.bed 不受影響 | Purity 估計失效 |
| F1 不退步 | ISM HP 維度完全不可用 |
| 執行更快 | 私有 germline variants 被錯誤分類為 SOMATIC |

### 2.5 第一版的根本問題

**二元強制分類不合理**：將所有非 PON variants 統一標為 SOMATIC，忽略了三種實際類別：
1. **PON germline**：在 PON 資料庫中 → 正確
2. **True somatic**：真正的腫瘤突變 → 應該被識別，而非預設
3. **Non-PON germline**：私有罕見 germline variants → 被錯誤強制為 SOMATIC

---

## 3. 第二版修正：二次 Phasing + 正確分類

### 3.1 核心概念

**三層 Variant 分類**：

| 類別 | 定義 | Phasing 角色 | VCF GT 格式 | Haplotag HP |
|------|------|-------------|-------------|-------------|
| PON Germline | PON 資料庫確認 | 驅動 phasing | `0\|1` / `1\|0` | HP:i:1 / HP:i:2 |
| True Somatic | somaticCalling() 識別 | 被動 phased | `0\|0` + GT2 | HP:i:11 / HP:i:21 |
| Non-PON Unknown | 不在 PON，非 somatic | 不干擾 phasing | `0/1`（unphased） | 由附近 germline 決定 |

### 3.2 執行流程

```
Step 1: addEdge()                      → 所有 variants 建構 graph
Step 2: convertNonGermlineToSomatic()  → 暫時標記非 PON 為 SOMATIC
Step 3: phasingProcess()               → 第一次 phasing（PON germline 主導）
Step 4: resetNonPonOrigin()            → [新函數] 重置非 PON 回 UNDEFINED
Step 5: somaticCalling()               → 用 clean graph 正確識別 true somatic
Step 6: 清空 posPhasingResult          → 清除第一次結果
Step 7: phasingProcess()               → 第二次 phasing（正確 origin）
Step 8: exportPhasingResult()          → 正確 GT/GT2 格式
```

**Haplotag 不需要改動** — VCF 格式正確後，haplotag 自然產出正確 HP tags。

### 3.3 為何此流程可行（程式碼驗證）

| 條件 | 驗證結果 | 依據 |
|------|---------|------|
| Edge 資訊在 phasing 後保留 | ✅ | phasingProcess() 不修改 edge counts（PhasingGraph.cpp） |
| somaticCalling() 不依賴 origin | ✅ | 三點邊投票基於 edge pattern，不檢查 origin（line 861-989） |
| 二次 phasing 已有先例 | ✅ | highPurity 路徑（PhasingProcess.cpp:178-197）已實現相同模式 |
| posPhasingResult 可清空重建 | ✅ | highPurity 路徑 line 194 示範清空後重跑 |

### 3.4 需要修改的檔案

| 檔案 | 修改 | 行數估計 |
|------|------|---------|
| `PhasingGraph.cpp` | 新增 `resetNonPonOrigin()` 函數 | ~7 行 |
| `PhasingGraph.h` | 宣告 `resetNonPonOrigin()` | ~1 行 |
| `PhasingProcess.cpp` | 修改 PON-only 流程，加入二次 phasing | ~15 行 |

**Haplotag 相關檔案（Haplotag.cpp, HaplotagProcess.cpp）不需要修改。**

### 3.5 設計推論

**為何不直接在第一次 phasing 後 exportPhasingResult？**

因為第一次 phasing 的 `posPhasingResult.somatic` flag 是基於 `convertNonGermlineToSomatic()` 設定的——所有非 PON 都標為 somatic。`exportPhasingResult()` 根據此 flag 決定 GT vs GT2 格式。如果不重新分類就 export，所有非 PON variants 都會得到 GT=`0|0`+GT2 格式，haplotag 仍然會產出錯誤的 HP tags。

**為何需要第二次 phasingProcess？**

因為 `edgeConnectResult()`（phasing 核心）根據 variant origin 計算不同的 edge weight（`findBestEdgePair()` line 191-203）。第一次 phasing 時所有非 PON 都是 SOMATIC（低 weight），第二次 phasing 時 somaticCalling() 已正確分類——non-PON germline 恢復正常 weight，true somatic 保持低 weight。第二次 phasing 結果更準確。

**為何不在 somaticCalling() 之後只更新 somatic flag？**

理論上可以（簡化方案），但 phasing 結果會略有偏差——被重分類為 germline 的 variants 在第一次 phasing 中 edge weight 被低估。二次 phasing 確保結果完全正確。

---

## 4. 為何原始 HP Tag 對 ISM 和生物資訊研究不好

### 4.1 ISM 的 HP 機制

**ReadParser.cpp 的 HP 映射**：
```
HP:i:1  → "1"   → HP1-family
HP:i:2  → "2"   → HP2-family
HP:i:11 → "1-1" → HP1-family
HP:i:21 → "2-1" → HP2-family
HP:i:33 → "3"   → 排除
```

**HP_Ratio 計算**（RegionProcessor.cpp）：
```
HP_Ratio = (HP1Family_N + 0.001) / (HP1Family_N + HP2Family_N + 0.002)
```

**問題**：ISM 映射本身正確，但 haplotag 在 somatic 位點把**所有 reads 推到同一 HP family**（因為 GT=`0|0` 的 refHaplotype=UNDEFINED 使 REF reads 無法投票）。

### 4.2 影響量化

| ISM 指標 | 正確值（Paired） | 錯誤值（TO-PON-Only） | 影響 |
|----------|---------------|---------------------|------|
| HP 二分群 | 67.1% | **0.15%** | 完全崩壞 |
| HPMergedSig | 37.4% | **0.003%** | HP 統計檢驗失效 |
| HPFineSig | 61.2% | **0.02%** | HP 精細檢驗失效 |
| Potential_LOH | ~47% | **99.9%** | 假 LOH 膨脹 |
| VerificationClass 改變 | — | 7.8% 位點降級 | Strong/Weak → Noise |

**影響鏈**：
1. HP 二分群崩壞 → VerificationClass 系統性降級
2. Potential_LOH 膨脹 → QS LOH penalty 失效（TO mode 已被迫停用）
3. HP 維度分析不可用 → ISM 失去核心分析維度
4. 未來改進被封死 → 任何依賴 HP 的新特徵都無法在 TO mode 使用

### 4.3 對其他工具的影響

| 工具 | 問題 |
|------|------|
| WhatsHap | 只認 HP:i:1/2，HP:i:11/21/33 被忽略 → phased rate 低估 |
| samtools split | HP:i:21 被分到獨立 group，非標準分組 |
| IGV | HP:i:21 無法正確顏色編碼 |
| Allele-specific methylation | 無法用 HP 分割兩個 haplotype 的甲基化模式 |
| Subclone analysis | 無法用 HP 區分不同亞系 |

### 4.4 生物學角度

**Germline haplotype 和 somatic status 是正交屬性**：

- 一條 read 來自哪條父母染色體（HP:i:1 or 2）由遺傳決定
- 這條 read 是否攜帶 somatic mutation 由腫瘤演化決定
- 兩者**獨立**，不應該用 somatic tag（HP:i:21）覆蓋 germline tag（HP:i:2）

**例子**：在一個 somatic variant 位點（GT2=`.|1`，mutation on HP2）：
- HP1 上的 reads：全部是 REF（HP1 沒有突變）→ 應該標 HP:i:1
- HP2 上的 reads：部分 ALT（腫瘤細胞），部分 REF（正常細胞）→ 應該標 HP:i:2
- 當前行為：所有 ALT reads → HP:i:21，所有 REF reads → 無 HP（UNDEFINED）

---

## 5. sSNV Benchmark 驗證框架

### 5.1 是否可以用 sSNV Benchmark 驗證 F1？

**可以，但需要注意 F1 的局限性。**

**可用的 benchmark 基礎設施**：
- Truth set: SEQC2 `high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`
- TP/FP VCFs: 已用 `bcftools isec` 分割（30,490 TP + 4,822 FP）
- ISM pipeline: `run_batch_vcf_analysis.sh` → `significance_summary.csv`
- F1 計算: `scripts/pipeline/steps/03_filter_analysis.py`

**F1 的測量方式**：
```
TP = SEQC2 truth variants 中被 caller 呼叫的（30,490）
FP = caller 呼叫但不在 truth 中的（4,822）
ISM SuggestFilter=True → 建議移除的 variants
Precision = 未被 SuggestFilter 標記的 TP / (未被標記的 TP + 未被標記的 FP)
Recall = 未被 SuggestFilter 標記的 TP / 30,490
F1 = 2 * Precision * Recall / (Precision + Recall)
```

### 5.2 F1 預期不會大幅變化的原因

**SuggestFilter 不直接使用 HP 統計量**。SuggestFilter 的決策依據：
1. `Significant` = GlobalP < 0.05 AND CramersV > threshold
2. `LocalBestCluster` = methylation clustering 的最佳分群
3. `AlleleDelta` = allele frequency difference between clusters

這些指標都不依賴 HP tag。因此修正 HP tag 對 SuggestFilter F1 的影響 < 0.005。

**但這不代表修正無意義**（見 5.3）。

### 5.3 F1 以外的決定性指標

修正的價值不在 SuggestFilter F1，而在以下指標：

| 指標 | 測量方式 | 預期改善 | 為何重要 |
|------|---------|---------|---------|
| **HP 二分群恢復** | Both HP1>0 AND HP2>0 rate | 0.15% → ~65% | HP 維度分析的前提條件 |
| **HPMergedSig 恢復** | HPMergedSig=True rate | 0.003% → ~30% | HP 統計檢驗效力恢復 |
| **HP_Ratio 正常化** | extreme rate (<0.1 or >0.9) | 99.9% → ~50% | 接近 Paired 的正確分佈 |
| **Potential_LOH 修正** | Potential_LOH=True rate | 99.9% → ~50% | 消除假 LOH 膨脹 |
| **VerificationClass 修正** | Strong + Weak rate | ~50% → ~60% | 恢復正確的位點分類 |
| **QS LOH penalty 可重啟** | QS AUC with LOH penalty | N/A → > 0.63 | 打通 TO mode QS 改進空間 |

### 5.4 十項充要條件

| # | 條件 | 閾值 | 現狀 | 修正後預期 |
|---|------|------|------|-----------|
| C1 | HP 二分群恢復 | Both HP > 50% | **FAIL** (0.15%) | PASS (~65%) |
| C2 | HPMergedSig 恢復 | > 20% | **FAIL** (0.003%) | PASS (~30%) |
| C3 | HP_Ratio 正常化 | extreme < 65% | **FAIL** (99.9%) | PASS (~50%) |
| C4 | Potential_LOH 不膨脹 | < 65% | **FAIL** (99.9%) | PASS (~50%) |
| C5 | F1 不退步 | >= Baseline - 0.005 | PASS (0.9291) | PASS |
| C6 | TP loss 不增加 | <= 0.5% | PASS (0.2%) | PASS |
| C7 | FP removal 不減少 | >= 4.7% | PASS (5.2%) | PASS |
| C8 | LOH.bed 不變 | Jaccard >= 0.99 | PASS (1.0) | PASS |
| C9 | Phasing N50 不退步 | >= 4,061 | PASS (8,109) | PASS |
| C10 | Allele 統計不變 | r > 0.99 | PASS (1.0) | PASS |

**判定規則**：
- **必要條件**（C1-C4）：HP 修正有效
- **安全條件**（C5-C7）：不引入退步
- **穩定條件**（C8-C10）：不影響其他模組
- **10/10 PASS** = 方案「更好」

### 5.5 sSNV Benchmark 執行計畫

```bash
# Phase 1: 編譯修改版 LongPhase-TO
cd /big7_disk/liaoyoyo2001/longphase-to-mod && make -j$(nproc)

# Phase 2: 執行二次 Phasing
./longphase-to phase --pon-only-phasing \
    --snp-file /big7_disk/liaoyoyo2001/data/HCC1395_5kHz/snv.vcf.gz \
    --bam-file /big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam \
    --reference /big7_disk/liaoyoyo2001/data/ref/GRCh38_no_alt_analysis_set.fasta \
    --pon-file [PON] --strict-pon-file [STRICT_PON] \
    -t 48 --ont -o output/pononly_v2/tumor_phased

# Phase 3: Haplotag（不改動 haplotag 程式碼）
./longphase-to haplotag \
    -s output/pononly_v2/tumor_phased.vcf \
    -b /big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam \
    -r /big7_disk/liaoyoyo2001/data/ref/GRCh38_no_alt_analysis_set.fasta \
    --tagSupplementary -t 48 -o output/pononly_v2/tumor_tagged

# Phase 4: ISM (TP + FP)
ISM_BIN=/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod
${ISM_BIN} --tumor-bam output/pononly_v2/tumor_tagged.bam \
    --normal-bam [NORMAL_BAM] --reference [REF] \
    --vcf [TP_VCF] --output-dir output/ism_pononly_v2_tp \
    --threads 48 --window-size 5000 --distance-metric BERNOULLI
# (同上但用 FP_VCF 到 ism_pononly_v2_fp)

# Phase 5: 計算所有指標（C1-C10）
python3 scripts/pipeline/steps/03_filter_analysis.py \
    --tp output/ism_pononly_v2_tp --fp output/ism_pononly_v2_fp
```

---

## 6. 自我論證：質疑與反駁

### 質疑 1：「F1 幾乎不變（0.9282 vs 0.9291），為什麼要改？」

**反駁**：

F1 不變是因為 SuggestFilter 不直接使用 HP 統計量。但 HP 是 ISM 架構的核心支柱之一：

- HPMergedSig 從 29.3% 崩壞到 0.003%（3,000 倍退化）
- HP 二分群從 55.9% 崩壞到 0.15%（370 倍退化）
- 不修正 = **放棄整個 HP 分析維度**

**生物資訊佐證**：在 allele-specific analysis 中，haplotype balance 是核心指標。HP_Ratio 偏移代表我們觀測到的「表觀遺傳 × 基因型」關聯是被 phasing artifact 扭曲的，而非真實的生物信號。任何基於 HP 的下游分析（ASM, LOH, subclone）都需要正確的 HP tag。

### 質疑 2：「為何不改 ISM 而要改 LongPhase-TO？」

**反駁**：

ISM 的 HP 映射已經正確（HP:i:21 → HP2-family）。問題在 haplotag 把 somatic 位點**所有 reads 推到同一 HP**。

```
即使修改 ISM 的映射邏輯：
  HP:i:21 → "2"（直接當 germline HP2）
結果：
  somatic 位點仍然所有 reads 都是 HP:i:21
  → HP_Ratio 仍然 = 0.0 或 1.0
  → HP 二分群仍然崩壞
```

**根因在上游**：haplotag 的 HP assignment 邏輯決定了 reads 的 HP tag。ISM 無論怎麼改映射都無法修復「所有 reads 同 HP」的問題。

**生物資訊佐證**：data pipeline 的最佳實踐是「在產生錯誤的地方修正」（fix at source），而非讓每個下游工具各自 workaround。修正 phasing 讓所有下游工具（ISM, WhatsHap, methylation analysis, subclone analysis）同時受益。

### 質疑 3：「N50 翻倍，但怎麼知道 phasing 是正確的？」

**反駁**：

多個獨立指標一致指向 phasing 品質改善：

1. **N50 +99.7%**（4,061 → 8,109）：phase block 更長，代表 phasing 一致性更好
2. **Phased rate +23.6pp**（54.9% → 78.5%）：更多 variants 成功 phased
3. **Somatic bias 消除**（17.3:1 → 1:1）：circular dependency 被修正
4. **LOH.bed Jaccard = 1.0**：LOH region 判定完全不變 → phasing 改動不影響 LOH
5. **SEQC2 LOH Jaccard = 0.8470**：外部 truth set 驗證一致

**進一步驗證方法**（如需）：
- Switch Error Rate（需 truth phasing set）
- Germline variant phasing concordance vs paired mode
- Haplotype-specific read depth balance in non-LOH regions

**生物資訊佐證**：在 short-read phasing 中，N50 和 phased rate 是評估 phasing 品質的標準指標（WhatsHap paper）。Long-read phasing 的 N50 翻倍代表 somatic contamination 曾嚴重破壞 phase continuity。

### 質疑 4：「只在 HCC1395 一個樣本測試」

**承認**：這是目前最大的弱點。

**緩解計畫**：
- **P0**（必須）：COLO829（另一 paired 樣本）+ H1437（純 TO）
- **P1**（建議）：H2009（高 FP rate）+ HCC1937（不同組織）
- **P2**（理想）：HCC1395 低純度系列（25%, 50%, 75%）+ HG002 正常對照

**判定**：至少 3 個樣本全部通過 C1-C10。

### 質疑 5：「二次 phasing 增加執行時間」

**評估**：

| 模式 | Phasing 時間 | 相對 baseline |
|------|-------------|-------------|
| Baseline | 2,693s | 1.00x |
| PON-only（單次） | 1,976s | 0.73x |
| 二次 phasing（預估） | ~3,952s | 1.47x |

增加 ~20 分鐘（全基因組），對臨床/研究場景可接受。

**生物資訊佐證**：在 somatic variant calling 中，迭代式 refinement（如 GATK 的 Base Quality Score Recalibration 多次 pass）是常見做法。二次 phasing 的概念類似：第一次建立 scaffold，第二次用正確分類精煉。

### 質疑 6：「二次 phasing 是否保證比 baseline 更好？可能引入新問題？」

**分析**：

二次 phasing 的每個步驟都有明確的改善邏輯：
1. **第一次 phasing**（clean graph）：N50 翻倍 → phasing scaffold 更準確
2. **somaticCalling() on clean graph**：比在 contaminated graph 上更準確地識別 somatic
3. **第二次 phasing**（correct origin）：germline + correct somatic → 最終結果接近理想

**可能的問題**：
- somaticCalling() 在 clean phasing 後可能識別出不同的 somatic set
- 第二次 phasing 的 N50 可能略低於第一次（因為 germline leakthrough 不再被強制為 SOMATIC）
- 但 phasing 結果的**正確性**更高（正確的 origin classification）

**風險緩解**：C1-C10 十項條件是全面的驗證框架。如果任何條件 FAIL，可以回溯分析具體原因。

---

## 7. 結論

### 7.1 方案選擇

| 方案 | Phasing 品質 | HP Tag 正確性 | 功能保留 | 推薦 |
|------|-------------|-------------|---------|------|
| Baseline（原始） | 差（N50=4,061） | 有 bias（17.3:1） | 完整 | ✗ |
| PON-only v1（單次 phasing） | 好（N50=8,109） | 崩壞（99.9% extreme） | 部分失效 | ✗ |
| **二次 Phasing v2** | **好** | **預期正確** | **完整** | **✓** |

### 7.2 驗證路線圖

```
Phase 1: 實作二次 phasing（修改 3 個檔案）
Phase 2: HCC1395 驗證（C1-C10）
Phase 3: 跨樣本驗證（COLO829 + H1437 + H2009）
Phase 4: 低純度驗證（HCC1395 25/50/75%）
```

### 7.3 成功定義

**10/10 充要條件全部 PASS** + **至少 3 個樣本一致** = 方案確認「更好」。

---

## 附錄 A：關鍵檔案路徑

| 檔案 | 說明 |
|------|------|
| `/big7_disk/liaoyoyo2001/longphase-to-mod/PhasingProcess.cpp` | 核心流程修改 |
| `/big7_disk/liaoyoyo2001/longphase-to-mod/PhasingGraph.cpp` | graph 邏輯 + resetNonPonOrigin() |
| `/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp` | haplotag（不修改） |
| `/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/` | Baseline 結果 |
| `/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly/` | PON-only v1 結果 |
| `/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/` | TP/FP VCFs |
| `/big8_disk/data/HCC1395/SEQC2/` | SEQC2 truth set |
| `research/loh_investigation/reports/20260403_pon_only_phasing_verification_report.md` | Step 1 報告 |
| `research/loh_investigation/reports/20260403_pon_only_haplotag_ism_verification_report.md` | Step 2 報告 |

## 附錄 B：現有數據完整統計

### B.1 三組 F1 比較

| 指標 | Paired | TO-Baseline | TO-PON-Only v1 |
|------|--------|-------------|----------------|
| TP 總數 | 30,476 | 30,476 | 30,476 |
| FP 總數 | 4,822 | 4,822 | 4,822 |
| SuggestFilter TP | 113 | 112 | 73 |
| SuggestFilter FP | 235 | 250 | 252 |
| Precision | 0.8688 | 0.8691 | 0.8693 |
| Recall | 0.9963 | 0.9963 | 0.9976 |
| **F1** | **0.9282** | **0.9284** | **0.9291** |

### B.2 三組 HP 分佈比較

| 指標 | Paired | TO-Baseline | TO-PON-Only v1 |
|------|--------|-------------|----------------|
| HP_Ratio median | 0.5000 | 0.8358 | 0.0000 |
| HP_Ratio extreme rate | 46.8% | 60.0% | 99.9% |
| Both HP > 0 | 67.1% | 55.9% | 0.15% |
| HPMergedSig | 37.4% | 29.3% | 0.003% |
| HPFineSig | 61.2% | 60.1% | 0.02% |
| Potential_LOH | ~47% | ~60% | 99.9% |

### B.3 Phasing 品質比較

| 指標 | Baseline | PON-only v1 |
|------|----------|-------------|
| Phase block N50 | 4,061 | 8,109 |
| Phased rate | 54.9% | 78.5% |
| Somatic HP bias | 17.3:1 | 消除 |
| LOH.bed Jaccard | — | 1.0000 |
| SEQC2 LOH Jaccard | 0.8470 | 0.8470 |
