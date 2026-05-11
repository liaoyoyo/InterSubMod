---
id: ism-kb-03-pipelines-index
name: "Pipelines 目錄索引"
description: "InterSubMod 三條分析 pipeline 索引：paired_full（canonical benchmark）、paired_pileup（模型變體）、tumor_only（有 self-phasing bias）。含選擇指引與對比。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "pipeline directory structure and selection guidance"
related_ids:
  - ism-kb-03-pipelines-paired-full
  - ism-kb-03-pipelines-paired-pileup
  - ism-kb-03-pipelines-tumor-only
  - ism-kb-03-pipelines-pipeline-comparison
  - ism-kb-05-data-formats-output-directory-layout
  - ism-kb-03-pipelines-f1-baseline-canonical
tags: [pipeline, index, paired, tumor-only, benchmark]
canonical_paths: [03_pipelines/00_index.md]
alias_paths: []
---

# Pipelines 目錄索引

- 一句結論：ISM 有 3 條 pipeline；**canonical benchmark = paired_full**，TO 模式結果不可直接對比 paired
- 適用對象：初次使用 ISM 者、選擇 pipeline 前
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big7_disk/liaoyoyo2001/InterSubMod/knowledge/03_pipelines/
  ```

---

## 三條 Pipeline 速查

| Pipeline | VCF 來源 | Normal BAM | Haplotag 工具 | 可信度 | 典型用途 |
|----------|---------|-----------|---------------|--------|----------|
| **paired_full** | ClairS paired full | ✅ 需要 | LongPhase-S | 🟢 Canonical | 論文 benchmark、F1 主表 |
| **paired_pileup** | ClairS paired pileup | ✅ 需要 | LongPhase-S | 🟡 模型評估 | 評估 pileup 模型效能 |
| **tumor_only** (TO) | ClairS-TO | ❌ 不用 | LongPhase-TO (self-phasing) | 🔴 有偏差 | TO 場景探索；結論須謹慎 |

---

## 選擇決策樹

```
  有 Normal BAM 配對嗎？
   │
   ├── 是 ──┬── 要 canonical benchmark ──→ paired_full
   │        └── 要評估 pileup 模型 ──→ paired_pileup
   │
   └── 否 ──→ tumor_only（注意 self-phasing bias，非 canonical）
```

---

## 文件列表

| 檔案 | 主題 | 關鍵結論 |
|------|------|----------|
| [01_paired_full.md](01_paired_full.md) | Canonical pipeline（paired + full ClairS VCF） | Phase 1A ΔF1 locked（見 05 SoT） |
| [02_paired_pileup.md](02_paired_pileup.md) | Pileup 模型變體 | 變異數 < full；評估用 |
| [03_tumor_only.md](03_tumor_only.md) | TO 模式（無 Normal） | 94.6% HP1 偏差；ΔF1=-0.0206（NEGATIVE） |
| [04_pipeline_comparison.md](04_pipeline_comparison.md) | 三條對比表 | 何時用哪條、delta F1 對照 |
| [05_f1_baseline_canonical.md](05_f1_baseline_canonical.md) | ★ **F1 / ΔF1 SoT** | 所有 ΔF1 數字的唯一權威（含完整 provenance） |

---

## 核心差異

### VCF 來源差異
- **paired_full**：ClairS paired mode 的 full model 輸出（含 pileup + full 複合）
- **paired_pileup**：只用 pileup model（模型較簡單、變異數較少）
- **TO**：ClairS-TO（tumor-only caller，依賴 PON 過濾 germline）

### Haplotag 差異
- **paired (S)**：LongPhase-S 用 normal BAM 作 phasing reference → 精確
- **TO**：LongPhase-TO self-phasing → HP1/HP2 分配 94.6% 偏向 HP1（已知 bias）

### 輸出目錄差異
- `output/canonical/<sample>/paired_full/`
- `output/canonical/<sample>/paired_pileup/`
- `output/canonical/<sample>/TO/`

---

## 常見誤區

- ❌ **誤把 TO F1 與 paired F1 直接比較** — 不同 VCF 來源 + phasing 方式，差異大部分來自 pipeline 差異
- ❌ **誤把 paired_pileup 當 paired_full 的子集** — pileup model 變異數較少但**不是**完全包含於 full
- ❌ **假設 ISM 對 TO 有正 F1 增益** — 實測 ΔF1=-0.0206（Phase 1A 已凍結）

---

## 相關資源

- 執行流程：[../06_workflows/02_full_vcf_analysis.md](../06_workflows/02_full_vcf_analysis.md)
- F1 計算：[../08_truth_and_benchmark/02_f1_calculation.md](../08_truth_and_benchmark/02_f1_calculation.md)
- Phase 1A locked 結論：[../../docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md)
