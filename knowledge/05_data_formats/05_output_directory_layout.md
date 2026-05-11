---
id: ism-kb-05-data-formats-output-directory-layout
name: "Output Directory Layout"
description: "output/ 頂層結構：canonical/（基線 7×3）、synthesis/（觀察）、archive/（歷史）、multilayer_hp_benchmark/ 等；信任度分級。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "output dir against output/OUTPUT_INDEX.md and docs/data_specs/20260414_output資料信任度與生命週期_01.md"
related_ids:
  - ism-kb-05-data-formats-index
  - ism-kb-03-pipelines-index
  - ism-kb-03-pipelines-paired-full
tags: [data-format, output, layout, canonical, archive]
canonical_paths: [05_data_formats/05_output_directory_layout.md]
alias_paths: []
---

# Output Directory Layout

- 一句結論：`output/` 下 canonical/ 是 7×3 基線（信任度 A），synthesis/ 是研究探索（信任度 B），archive/ 是歷史（信任度 C-）
- 適用對象：定位 canonical run、區分實驗與基線、理解信任度
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big7_disk/liaoyoyo2001/big7_disk_output/
  cat /big7_disk/liaoyoyo2001/big7_disk_output/OUTPUT_INDEX.md 2>/dev/null | head -40
  ```

---

## 頂層結構

| 目錄 | 分類 | 說明 | 信任度 |
|------|------|------|--------|
| `canonical/` | Canonical | 7 樣本 × 3 模式 ISM baseline（19 runs） | 🟢 A |
| `synthesis/` | Synthesis | 觀察、研究 rounds、批次、concluded、manifests | 🟡 B |
| `bip8_output_archive/` | Archive | bip8 歷史（39 dirs） | 🔴 C- |
| `big8_output_archive/` | Archive | big8 歷史（5 dirs） | 🔴 C- |
| `multilayer_hp_benchmark/` | Standalone | HP 比較 benchmark | 🟡 B |
| `hcc1395_normal_pilot*/` | Pilot | Normal BAM 整合試驗 | 🟡 B |
| `kde_smoke_test/` | Diagnostic | KDE 修正驗證 | 🟡 B |

**完整索引**：`output/OUTPUT_INDEX.md`
**信任度規範**：`docs/data_specs/20260414_output資料信任度與生命週期_01.md`

---

## Canonical/ 結構（信任度 A）

```
canonical/
├── HCC1395/
│   ├── paired_full/<timestamp>/
│   ├── paired_pileup/<timestamp>/
│   └── TO/<timestamp>/
├── HCC1395_DORADO/
│   └── <same>
├── COLO829/
│   └── <same>
├── H1437/
│   └── <same>
├── H2009/
│   └── <same>
├── HCC1937/
│   └── <same>
└── HCC1954/
    └── <same>
```

**總計**：7 樣本 × 3 模式 = 21 目錄（實際 19 個因 DORADO 未全跑）

每個 `<timestamp>/` 包含：
- `longphase_s/`（LongPhase 輸出）
- `ism_final/`（ISM 輸出含 distance matrix, heatmaps, cluster）
- `figures/`（視覺化）

---

## Synthesis/ 結構（信任度 B）

```
synthesis/
├── observation_workspaces/        # 系統性觀察（O-系列 .md）
├── research_rounds/               # 自動化研究迴圈輸出
├── batches/                       # 批次 benchmark 結果
├── concluded/                     # 已結案的研究 artifact
├── legacy_partials/               # 部分完成舊跑
└── manifests/                     # 實驗 manifests
```

**典型用途**：
- 跨樣本分析（all_region_rows.tsv.gz 居於 `observation_workspaces/20260327_loh_round1_cross_sample_audit/`）
- Round-based 假說驗證產出

---

## 信任度分級與生命週期

### 🟢 A - Canonical（論文標準）
- 目錄：`canonical/`
- 特性：
  - 參數穩定、版本固定
  - 可重現性 100%
  - 論文 benchmark 主表唯一來源
- 保留：永久

### 🟡 B - Synthesis / Pilot（研究探索）
- 目錄：`synthesis/`, `multilayer_hp_benchmark/`, `*_pilot*/`, `kde_smoke_test/`
- 特性：
  - 特定研究假說驗證
  - 參數可能變動
  - 可能需要重跑
- 保留：relevant 研究結案後 6-12 個月

### 🔴 C- Archive（歷史）
- 目錄：`*_output_archive/`
- 特性：
  - 舊版 binary 或已矯正錯誤
  - 僅作 reproducibility trace
  - 不應用於現行 benchmark
- 保留：至少 1 年（compliance）

---

## 典型 canonical run 路徑範例

**HCC1395 paired_full**：
```
/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/
  20260314_HCC1395_paired_full_full_complete_matrix/
    ├── longphase_s/
    │   ├── germline_phased_merged.vcf.gz
    │   ├── per_region_features.tsv
    │   ├── significance_summary.csv       # ← 59 欄 ISM 輸出
    │   └── master_dataset.csv
    ├── ism_final/
    │   ├── distance_matrix.h5
    │   ├── heatmaps/
    │   └── cluster labels
    └── figures/
```

---

## 查找 canonical run 的快速方法

```bash
# 找特定樣本 × pipeline 的最新 canonical
ls -td /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/*/ | head -1

# 找 7 樣本 × 3 模式 × 最新 run
for s in HCC1395 HCC1395_DORADO COLO829 H1437 H2009 HCC1937 HCC1954; do
  for m in paired_full paired_pileup TO; do
    d=$(ls -td /big7_disk/liaoyoyo2001/big7_disk_output/canonical/$s/$m/*/ 2>/dev/null | head -1)
    echo "$s/$m: $d"
  done
done
```

---

## 符號連結

專案根 `output -> /big7_disk/liaoyoyo2001/big7_disk_output/`

```bash
readlink /big7_disk/liaoyoyo2001/InterSubMod/output
# → /big7_disk/liaoyoyo2001/big7_disk_output/
```

---

## 相關

- 信任度規範：[../../docs/data_specs/20260414_output資料信任度與生命週期_01.md](../../docs/data_specs/20260414_output資料信任度與生命週期_01.md)
- Output 頂層索引：`output/OUTPUT_INDEX.md`
- VCF 來源：[04_vcf_sources.md](04_vcf_sources.md)
- Per-region 結構：[03_per_region_files.md](03_per_region_files.md)
