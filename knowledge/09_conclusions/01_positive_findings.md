---
id: ism-kb-09-conclusions-positive-findings
name: "Positive Findings 索引"
description: "已驗證為真實信號的 positive findings 索引：HPFineNGroups subclone marker、LOH×AF×Methylation、Self-Phasing 根因、Zone-Aware H1、V5 Somatic Fallback 等。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "positive findings against MEMORY Active Research and Concluded"
related_ids:
  - ism-kb-09-conclusions-index
  - ism-kb-07-derived-features-index
  - ism-kb-07-derived-features-hpfinengroups
  - ism-kb-07-derived-features-loh-af-methylation
  - ism-kb-09-conclusions-research-landscape-index
  - ism-kb-03-pipelines-f1-baseline-canonical
tags: [conclusions, positive, findings, index]
canonical_paths: [09_conclusions/01_positive_findings.md]
alias_paths: []
---

# Positive Findings 索引

- 一句結論：10+ 個已驗證為真實的 positive 結論；跳轉到 docs/reports/research_landscape/ 或 research/ 原始證據
- 適用對象：論文 results 章節、新進者了解可靠發現
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/research_landscape/
  ```

---

## Positive Findings 清單

### 🟢 P1：LOH × AF × Methylation（最強發現）
- **結論**：Inter AF→NGroups +0.705；7/7 樣本 p<10^-39；文獻空白
- **KB 詳細**：[../07_derived_features/04_loh_af_methylation.md](../07_derived_features/04_loh_af_methylation.md)
- **Research**：`research/loh_subclone_af/`（TO）, `research/loh_subclone_af_paired/`（Paired）
- **MEMORY**：`project_loh_subclone_af_methylation_positive`

### 🟢 P2：HPFineNGroups Subclone Marker
- **結論**：**新 canonical filter `NG=4 + AF<0.4 + NR≥80 NonLOH`** TP rate 92.81%（vs 舊 `N≥4+NR≥80 NonLOH` = 89.12%）；+3.7pp（HCC1954 +21pp 挽救）
- **⚠️ 注意**：舊 filter（僅兩條件 N≥4+NR≥80）TP rate 為 89.12%；新 filter 需加 `AF<0.4` 才達 92.81%
- **KB 詳細**：[../07_derived_features/01_hpfinengroups.md](../07_derived_features/01_hpfinengroups.md)
- **⚠️ 警告**：2026-04-21 flag=on 下 N≥3 消失；需重驗
- **Research**：`research/F_hpfinengroups_deepening/`
- **Landscape**：`docs/reports/research_landscape/09_Part_B.md`
- **MEMORY**：`project_hpfinengroups_subclone_marker`

### 🟢 P3：Self-Phasing Circular Dependency
- **結論**：TO self-phasing 導致 62% LOH 消失
- **Landscape**：`docs/reports/research_landscape/02_Self_Phasing根因.md`
- **MEMORY**：`project_self_phasing_causal_chain_confirmed`

### 🟢 P4：V5 Somatic Fallback Haplotag
- **結論**：當前最佳；AMB% 17.5→8.0%；F1=0.7154；clean blocks 95%
- **取代**：V3-Fixed 已被 V5 取代
- **MEMORY**：`project_v5_somatic_fallback_verification`

### 🟢 P5：Zone-Aware H1（TP rate 差異真實）
- **結論**：Zone TP rate 有統計顯著差異（7 樣本一致）
- **但**：H3 QS 調整 NEGATIVE（characterization only）
- **KB 詳細**：[../07_derived_features/03_zone_aware_framework.md](../07_derived_features/03_zone_aware_framework.md)
- **Landscape**：`docs/reports/research_landscape/08_Zone_Aware.md`

### 🟢 P6：SNV-Methylation ASM
- **結論**：ASM 32-66% 但 FP>>TP 重疊大
- **MEMORY**：`project_snv_methylation_association`

### 🟢 P7：PON-Only Phasing 驗證
- **結論**：LOH.bed 不變；somatic bias 消除
- **MEMORY**：`project_pon_only_phasing_verification`

### 🟢 P8：Normal BAM 整合 Phase 2A（進行中但 POSITIVE signal）
- **結論**：HCC1395 pilot 97.3% sig；pending 7 樣本
- **狀態**：Ongoing（Phase 2 active）
- **MEMORY**：`project_normal_bam_progressive_copy`

---

## 與 Phase 1A locked F1 結論

⚠️ **注意**：上述 positive findings 多為 **characterization** 性質；F1 filter 主表結論（完整 provenance 見 [03_pipelines/05_f1_baseline_canonical.md](../03_pipelines/05_f1_baseline_canonical.md)）：
- paired_full ΔF1 = **+0.0112**（locked）
- TO ΔF1 = **-0.0206**（NEGATIVE locked）

詳見 [../03_pipelines/04_pipeline_comparison.md](../03_pipelines/04_pipeline_comparison.md)

---

## Research Landscape 權威文件對照

| Positive | Landscape 文件 |
|----------|---------------|
| LOH × AF | `07_LOH_CN_AF_研究總整理.md` |
| HPFineNGroups | `09_Part_B.md` |
| Self-Phasing | `02_Self_Phasing根因.md` |
| Zone-Aware | `08_Zone_Aware.md` |
| ISM 分析價值 | `03_ISM分析價值界定.md` |

全索引：[04_research_landscape_index.md](04_research_landscape_index.md)

---

## 相關

- Characterization：[02_characterization_only.md](02_characterization_only.md)
- NEGATIVE：[03_concluded_negative.md](03_concluded_negative.md)
- Derived features：[../07_derived_features/](../07_derived_features/)
- Landscape 全索引：[04_research_landscape_index.md](04_research_landscape_index.md)
