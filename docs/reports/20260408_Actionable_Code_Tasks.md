<!--
建立時間: 2026-04-08 04:40
目標: 不依賴 PON-only 重跑、現在可執行的程式改進清單
處理範圍: C++ 核心 + Python 分析腳本 + 基礎設施
關聯檔案:
  - docs/reports/20260408_InterSubMod_Stage_Report_v1.md
  - src/core/ReadParser.cpp
  - src/core/RegionProcessor.cpp (QualityScorer)
  - src/core/FisherExact.cpp
-->

# 可執行程式工作清單

**不依賴 PON-only 全量重跑，現在就能開始**

---

## P0: 高優先（直接影響下游分析品質）

### P0-1: ReadParser somatic HP tag 過濾選項

**檔案**: `src/core/ReadParser.cpp:122-139`

**現狀**:
```cpp
case 11: info.hp_tag = "1-1"; break;  // somatic from HP1
case 21: info.hp_tag = "2-1"; break;  // somatic from HP2
case 33: info.hp_tag = "3";   break;  // unphased
```

**問題**: 在 self-phasing 下，HP:i:11/21 是循環標記（somatic variant phase 自己）。下游 LabelTest 用這些做分組會引入 circular bias。

**修改方案**: 新增 CLI flag `--ignore-somatic-hp`，啟用時：
- HP:i:11 → 視為 HP:i:1（保留 germline haplotype）
- HP:i:21 → 視為 HP:i:2
- HP:i:33 → 維持 HP:i:0（unphased）

**影響範圍**: ReadParser → RegionProcessor → LabelTest → FisherExact → SignificanceAnalyzer

**預估工時**: 2-3 小時（含測試）

---

### P0-2: QualityScore TO 模式根本重設計

**檔案**: `src/core/RegionProcessor.cpp:119-183`

**現狀**: TO mode 只停用 LOH penalty 和 verify bonus（`loh_penalty=0, verify_bonus=0`），但整體公式仍基於 coverage/CpG 等與 TP/FP 無關的指標。

**問題**: TO QS AUC=0.497（random），FN vs TP AUC=0.338（反轉）。

**修改方案 A（保守）**: 在現有框架內加入 caller_af 權重
```cpp
QualityScoreWeights get_tumor_only_weights() {
    QualityScoreWeights w{};
    w.loh_penalty = 0.0f;
    w.verify_bonus = 0.0f;
    // NEW: AF-based scoring for TO
    w.use_caller_af = true;
    w.af_bonus_high = 15.0f;   // AF > 0.3
    w.af_bonus_medium = 8.0f;  // AF > 0.15
    return w;
}
```

**修改方案 B（激進）**: 完全替換 TO QS 為 caller-feature-only 公式

**建議**: 先做方案 A，等 PON-only 重跑後再評估是否需要方案 B。

**前提**: 需要 VCF INFO 中的 AF 傳遞到 RegionProcessor。目前 VCF 的 caller 特徵是在 manifest building 時提取的，ISM C++ 本身不讀取。需要確認是否要在 C++ 層面加入 VCF INFO 解析。

**預估工時**: 4-6 小時（方案 A），需要設計決策。

---

## P1: 中優先（改善分析品質但不阻塞主線）

### P1-1: CramersV 多群框架修正

**檔案**: `src/core/FisherExact.cpp`

**現狀**: CramersV 基於 2×2 contingency table（HP1 vs HP2），93% 的 regions CramersV=0。

**問題**: 對 HP 分配不均或多群結構的 regions 完全失效。

**修改方案**: 
- 改用 R×C contingency table（R = cluster groups, C = methylation states）
- 或改用 HPFineNGroups 分群後的多群 CramersV

**預估工時**: 3-4 小時

---

### P1-2: PON-only haplotag 自動化 pipeline

**檔案**: 新建 `scripts/pipeline/steps/haplotag_pon_only.sh`

**功能**:
```bash
# For each sample:
longphase haplotag \
    --bam tumor.bam \
    --snp-file pon_only_phased.vcf.gz \
    --reference ref.fasta \
    --output tumor_pon_haplotagged.bam \
    --threads 112

# Then run ISM
inter_sub_mod \
    --tumor-bam tumor_pon_haplotagged.bam \
    [--normal-bam normal.bam] \
    --reference ref.fasta \
    --vcf caller_output.vcf.gz \
    --output-dir output/ \
    ...
```

**已有**: PON-only phased VCF 存在於 HCC1395 的 step03 中。需要確認其他 6 samples 是否也有。

**預估工時**: 2-3 小時（腳本 + 7 samples 配置）

---

### P1-3: Before/After 自動比較框架

**檔案**: 新建 `scripts/analysis/compare_before_after_hp_fix.py`

**功能**: 
- 載入重跑前後的 significance_summary.csv
- 計算每個 HP-dependent 特徵的 AUC 變化
- 生成 before/after heatmap
- 自動判定哪些特徵有顯著改善

**預估工時**: 2 小時

---

### P1-4: Gating 邏輯 TO 適配

**檔案**: `src/core/SignificanceAnalyzer.cpp:89-99`

**現狀**:
```cpp
result.passed_gate = expanded.alt.passed_gate || expanded.hp_family.passed_gate;
// ...
if (result.passed_gate) {
    // Phase 3: Structure Tests (PERMANOVA)
}
```

**問題**: TO mode 只有 22% regions passed_gate → ClusterPermanovaF 只在 22% 上計算。

**修改方案**: TO mode 下放寬 gating（如降低 hp_family 最低 read 要求），或新增 cluster-only gate（不依賴 HP）。

**注意**: ClusterPermanovaF AUC=0.512（即使在有效子集上也是 random），放寬 gate 不會改善區分力，但會讓更多 regions 有 PERMANOVA 值（對下游分析有用）。

**預估工時**: 1-2 小時

---

## P2: 低優先（nice-to-have）

### P2-1: Output 格式擴充 — caller feature passthrough

讓 ISM 直接輸出 VCF 中的 caller features (AF, DP, GQ)，避免需要在 Python 層面 merge。需要修改 VCF parser 和 output writer。

### P2-2: H2009 異質性診斷

Phase 1A 中唯一負向 sample。需要深度分析其 ISM 特徵分布、FP 來源、coverage pattern。

### P2-3: Platform normalization 初步分析

HCC1395（5kHz）vs HCC1395_DORADO（DORADO re-basecall）已有數據。比較兩者 ISM 特徵差異，為跨平台校準建立 baseline。

---

## 執行順序建議

```
Week 1:
  [P0-1] ReadParser somatic HP filter → 2-3h
  [P1-2] PON-only haplotag pipeline → 2-3h
  [P1-3] Before/After comparison framework → 2h
  → 啟動 7 samples haplotag + ISM 全量重跑 (background, ~8-12h)

Week 2:
  [P0-2] QS TO redesign (方案 A) → 4-6h
  [P1-1] CramersV 多群框架 → 3-4h
  → 等重跑完成 → 用 P1-3 自動比較

Week 3:
  [P1-4] Gating TO 適配 → 1-2h
  → 根據重跑結果決定後續
  [P2-x] 低優先項目按需處理
```

**總預估**: ~20-25 小時程式工作 + ~12 小時計算時間（重跑）
