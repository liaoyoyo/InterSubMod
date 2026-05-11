---
id: ism-kb-03-pipelines-tumor-only
name: "Tumor-Only (TO) Pipeline"
description: "ClairS-TO + LongPhase-TO + ISM；無需 Normal BAM。⚠️ 已知 self-phasing bias（94.6% HP1）；ΔF1 為負（-0.0206）。非 canonical；探索用途。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "pipeline against scripts/run_vcf_all_snv.sh; self-phasing bias against docs/reports/research_landscape/02_Self_Phasing根因.md"
related_ids:
  - ism-kb-03-pipelines-index
  - ism-kb-03-pipelines-paired-full
  - ism-kb-03-pipelines-pipeline-comparison
  - ism-kb-09-conclusions-concluded-negative
  - ism-kb-05-data-formats-vcf-sources
  - ism-kb-06-workflows-full-vcf-analysis
  - ism-kb-03-pipelines-f1-baseline-canonical
tags: [pipeline, tumor-only, to, clairs-to, longphase-to, self-phasing]
canonical_paths: [03_pipelines/03_tumor_only.md]
alias_paths: []
---

# Tumor-Only (TO) Pipeline

- 一句結論：無 Normal BAM 的 ClairS-TO + LongPhase-TO + ISM；⚠️ self-phasing bias 嚴重（94.6% reads 偏向 HP1），ΔF1=-0.0206（負增益），**非 canonical**
- 適用對象：無 Normal BAM 場景探索、Self-phasing bias 研究
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big8_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v3_fixed/
  ```

---

## ⚠️ 核心警告

**不要把 TO 結果與 paired 直接比較**。TO 有三大結構性差異：

1. **Self-phasing bias**：無 normal reference，LongPhase-TO 依 tumor 自身推 haplotype → 94.6% reads 標為 HP1
2. **Germline-somatic confound**：ClairS-TO 用 PON 過濾 germline，但無法完全區分，導致 HP tag 意義不穩
3. **ΔF1 為負**：全 7 樣本平均 ISM 增益 -0.0206（Phase 1A 已凍結）

詳細因果鏈：[../../docs/reports/research_landscape/02_Self_Phasing根因.md](../../docs/reports/research_landscape/02_Self_Phasing根因.md)

---

## 輸入鏈路

```
┌─────────────────┐       ┌──────────────┐      ┌──────────────┐
│ Tumor BAM       │       │ Reference     │      │ PON database  │
│ (unphased)      │       │ hg38.fa       │      │               │
└────────┬────────┘       └──────┬───────┘      └──────┬───────┘
         │                        │                    │
         ▼                        ▼                    ▼
┌──────────────────────────────────────────────────────────┐
│  ClairS-TO (tumor-only caller)                            │
│  → clairsto_tp.vcf.gz (uses PON to filter germline)       │
└──────────────────────────┬───────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│  LongPhase-TO (self-phasing, no normal reference)         │
│  → tumor_tagged.bam (HP:i:1/2/11/21/33)                   │
│  → phased.bed (LOH, note: 62% LOH lost in some settings)  │
└──────────────────────────┬───────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│  InterSubMod                                              │
│  ⚠ 94.6% HP1 bias → distance matrix 結構偏誤              │
└──────────────────────────────────────────────────────────┘
```

---

## 典型命令

```bash
./scripts/run_vcf_all_snv.sh \
  --caller-mode to \
  --mode all-with-w5000 \
  -v /big8_disk/liaoyoyo2001/longphase-to-mod/output/clairsto_v3fixed_work/clairsto_tp.vcf.gz \
  -t /big8_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v3_fixed/tumor_tagged.bam \
  -o output/canonical/<sample>/TO/<timestamp>
```

**建議加上**：
```bash
--germline-hp-only   # 降級體細胞 HP:i:11/21/33 為純 germline HP:i:1/2
                     # （已驗證機制正確但 HCC1395 TO 無 F1 改善）
```

---

## `--germline-hp-only` 旗標

- **位置**：`include/utils/ArgParser.hpp:127-128`
- **用途**：忽略 LongPhase-TO 產生的 somatic HP tags（HP:i:11/21/33），僅用 germline HP:i:1/2
- **Phase 1 結論（2026-04-21）**：機制正確，HCC1395 TO 全量無 TSV 特徵 AUC ≥0.02 改善；不進 Phase 2
- **⚠️ 副作用**：HPFineNGroups N≥3 在 flag=on 下消失，需重新驗證 subclone marker 結論

詳見：`09_conclusions/03_concluded_negative.md` 的 `ReadParser --germline-hp-only Phase 1 CONDITIONAL NEGATIVE`

---

## 輸出結構

同 paired_full，但目錄名為 `TO/`：
```
output/canonical/<sample>/TO/<timestamp>/
├── significance_summary.csv
├── TP/
└── FP/
```

欄位結構不變，但**部分欄位意義被 bias 污染**（尤其 HPMergedDelta、HPFineNGroups、Allele*）。

---

## 資料信任度

- **status**：🔴 非 canonical；characterization 用途
- **ClairS-TO Verdict 欄位（2026-04-20）**：內部校準 POSITIVE（Germline FP 96.1%、Somatic TP 91.8%），但 `Verdict_Germline 100%` 落 LowQual，對 ISM F1 無影響
- **歷史問題（2026-04-04 前）**：TO VCF 來源曾誤用 pileup symlink，已矯正（MEMORY: project_vcf_source_error_correction）

---

## 已知 TO canonical 現況

| 樣本 | TO canonical Run | 有效性 |
|------|-----------------|---------|
| HCC1395 | ✅ 20260314 baseline + 2026-04-04 VCF 矯正 | 有效（HP tag fix 後） |
| HCC1395_DORADO | ⚠ ΔF1 only +0.000476 | 低信噪比 |
| COLO829 | ⚠ 無 methylation (ONT_R10) | 部分結果 |
| H1437 | ✅ 5/6 zone TP 0.61-0.94 | 有效 |
| H2009 | ✅ 負向例（caller 已完美） | 有效 |
| HCC1937 | ✅ 標準情況 | 有效 |
| HCC1954 | ⚠ Amplicon artifact | 需額外過濾 |

詳見 [../02_samples/](../02_samples/) 各樣本文件。

---

## 相關

- 索引：[00_index.md](00_index.md)
- Self-phasing 根因：[../../docs/reports/research_landscape/02_Self_Phasing根因.md](../../docs/reports/research_landscape/02_Self_Phasing根因.md)
- TO FP 全貌：[../../docs/reports/research_landscape/01_TO_FP問題全貌.md](../../docs/reports/research_landscape/01_TO_FP問題全貌.md)
- 對照表：[04_pipeline_comparison.md](04_pipeline_comparison.md)
