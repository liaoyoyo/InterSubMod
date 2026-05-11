---
id: ism-kb-03-pipelines-pipeline-comparison
name: "三條 Pipeline 對比表"
description: "paired_full / paired_pileup / tumor_only 全面對照；包含 VCF 來源、haplotag 方式、ΔF1、輸出差異、使用建議。選擇 pipeline 前必讀。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "cross-pipeline comparison against individual pipeline docs"
related_ids:
  - ism-kb-03-pipelines-index
  - ism-kb-03-pipelines-paired-full
  - ism-kb-03-pipelines-paired-pileup
  - ism-kb-03-pipelines-tumor-only
  - ism-kb-03-pipelines-f1-baseline-canonical
tags: [pipeline, comparison, selection-guide, benchmark]
canonical_paths: [03_pipelines/04_pipeline_comparison.md]
alias_paths: []
---

# 三條 Pipeline 對比表

- 一句結論：paired_full 是 canonical，paired_pileup 是 full 的子集驗證，TO 不可與前兩者直接比 F1
- 適用對象：選擇 pipeline 前、寫論文 methods 時
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  # 快速看三條 canonical run 存在
  ls /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/
  ```

---

## 全面對照表

| 面向 | paired_full | paired_pileup | tumor_only (TO) |
|------|-------------|---------------|-----------------|
| **VCF 來源** | ClairS paired full | ClairS paired pileup | ClairS-TO |
| **VCF 變異數** | 多（full + pileup） | 少（僅 pileup） | 中 |
| **Normal BAM** | ✅ 必要 | ✅ 必要 | ❌ 不用 |
| **Haplotag 工具** | LongPhase-S | LongPhase-S | LongPhase-TO |
| **Haplotag bias** | 低 | 低 | 🔴 高（94.6% HP1） |
| **LOH BED 來源** | LongPhase-S phased.bed | 同 | LongPhase-TO phased.bed |
| **Self-phasing** | 否 | 否 | ✅ 是（依 tumor 自身） |
| **ISM ΔF1 (locked)** | **+0.0112** | ~類似 paired_full | **-0.0206** |
| **Canonical status** | 🟢 論文標準 | 🟡 副次 | 🔴 非 canonical |
| **輸出目錄** | `canonical/<sample>/paired_full/` | `canonical/<sample>/paired_pileup/` | `canonical/<sample>/TO/` |
| **適用場景** | Benchmark 主表 | pileup 模型評估 | 無 normal 場景探索 |

---

## 何時用哪條？

```
┌─ 論文 benchmark 主表 ────→ paired_full
├─ Pileup 模型 vs full 對照 ─→ paired_pileup
├─ 評估 ClairS 本身的效能 ──→ paired_full + paired_pileup 並跑
├─ 無 normal BAM 的場景 ────→ tumor_only（結果需加註 caveat）
└─ Self-phasing 研究 ───────→ tumor_only（自身就是研究對象）
```

---

## F1 / ΔF1 對照（Phase 1A locked）

👉 **完整 SoT（含 CI、樣本、方法、運行條件、限制、lock 日期）**：[05_f1_baseline_canonical.md](05_f1_baseline_canonical.md)

**速查**：

| Pipeline | ΔF1 | Status |
|----------|------|--------|
| paired_full | **+0.0112** [CI +0.0044, +0.0188] | 🟢 locked |
| paired_pileup | ~+0.01 範圍 | 🟡 未正式 lock |
| TO | **-0.0206** | 🔴 NEGATIVE locked |

**重要**：ΔF1 不跨 pipeline 比較；只能在同一 pipeline 內比較加與不加 ISM filter 的差異。Baseline F1 / method / reproducibility 規格見 SoT 文件。

---

## 輸出結構一致性

所有 3 條 pipeline 的輸出**結構相同**（差異只在目錄名與 VCF 來源）：

```
<pipeline_output>/
├── significance_summary.csv         # 59 欄，結構同
├── TP/
│   └── region_<N>/
│       ├── reads.tsv
│       ├── methylation/
│       └── distance_<METRIC>/
└── FP/
    └── <same>
```

欄位字典通用：[../05_data_formats/01_significance_summary_schema.md](../05_data_formats/01_significance_summary_schema.md)

---

## 跨 pipeline 欄位意義差異

雖然 CSV schema 相同，某些欄位在不同 pipeline 的「可信度」不同：

| 欄位 | paired_full | paired_pileup | TO |
|------|-------------|---------------|-----|
| HP_Ratio | 🟢 精確 | 🟢 精確 | 🔴 偏向 HP1 |
| HPFineNGroups | 🟢 characterization OK | 🟢 同左 | ⚠ flag=on 下全消失 |
| AlleleDelta | 🟢 | 🟢 | 🟡 受 HP bias 影響 |
| Potential_LOH | 🟢 | 🟢 | 🟡 62% LOH 消失 |
| Coverage_Multiple | 🟢（r≈0.83） | 🟢 | 🟢 |

---

## 批次對比腳本

```bash
# 所有 7 樣本 × 3 pipeline 的 F1 對比
python3 scripts/analysis/compare_benchmark_f1.py \
  --root output/canonical/ \
  --output-tsv benchmark_all.tsv

# 跨 pipeline 的 region 欄位對比
python3 scripts/analysis/build_loh_round1_cross_sample_audit.py
```

輸出：`all_region_rows.tsv.gz`（116 欄）— 見 [../05_data_formats/02_master_dataset_schema.md](../05_data_formats/02_master_dataset_schema.md)

---

## 論文寫作建議

**Methods 章節**應清楚陳述：
1. Primary benchmark：paired_full
2. Secondary analysis：paired_pileup（for pileup model ablation）
3. TO：only discussed as characterization study with caveats noted

**Results 章節不可**：
- 把 TO F1 與 paired F1 放同一表（會誤導讀者）
- 宣稱 TO 有「ISM F1 gain」（實測為負）

---

## 相關

- 各 pipeline 詳細：[01_paired_full.md](01_paired_full.md)、[02_paired_pileup.md](02_paired_pileup.md)、[03_tumor_only.md](03_tumor_only.md)
- CURRENT_FOCUS 權威：[../../docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md)
- Self-phasing 根因：[../../docs/reports/research_landscape/02_Self_Phasing根因.md](../../docs/reports/research_landscape/02_Self_Phasing根因.md)
