---
id: ism-kb-10-research-status-blockers-and-risks
name: "Blockers & Risks"
description: "當前阻塞：haplotag 重跑、expected_coverage hardcoded bug、COLO829 TO 無 methylation；merged 檔 AF 欄位陷阱 + HCC1395 phase1_new LOH 殘缺 (2026-04-23 新加)；風險清單與緩解。⚠ 2 週有效。"
status: active
last_verified: 2026-04-23
content_nature: runtime-fact
doc_type: reference
verified_scope: "blockers against MEMORY pending items"
related_ids:
  - ism-kb-10-research-status-index
  - ism-kb-10-research-status-next-milestones
  - ism-kb-07-derived-features-coverage-multiple
  - ism-kb-10-research-status-current-focus-snapshot
  - ism-kb-05-data-formats-merged-dataset-pitfalls
tags: [status, blockers, risks, pending]
canonical_paths: [10_research_status/04_blockers_and_risks.md]
alias_paths: []
---

# Blockers & Risks

> ⚠️ **此為 2026-04-22 快照**。阻塞狀態變動快；決策前先核對 `docs/CURRENT_FOCUS.md`

- 一句結論：3 個主要阻塞（haplotag 重跑、expected_coverage bug、COLO829 TO）；風險可控
- 適用對象：決策前、排除阻塞優先級排序
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  grep -A 2 "Pending\|Blocker\|阻塞" /big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md
  ```

---

## 當前主要阻塞

### B0：Merged 合成檔資料陷阱（2026-04-23 新增）⚠
- **影響**：跨樣本 AF 視覺化、LOH-aware 分析全受影響
- **陷阱 1**：`merged_7samples_paired_full_plus_hcc1395_to.tsv.gz` 的 `AF` 欄位**非 caller_af**（非 HCC1395 的 5 樣本 p75 皆 <0.06，視覺全堆底）
- **陷阱 2**：HCC1395 `kde_status='phase1_new'` rows 的 LOH annotation 殘缺（Potential_LOH=True 僅 5.7% vs 正常 58.8%）
- **偵測**：`df.groupby('sample')['AF'].quantile(0.75)` 若 <0.1 即可疑；`df.groupby('sample')['Potential_LOH'].mean()` 若單樣本偏離 median >30pp 即可疑
- **緩解**：改從 stale master archive 讀取（`all_region_rows.tsv.gz`），用 `caller_af` + `NumReads/bx` 重算
- **文件**：[05_data_formats/06_merged_dataset_pitfalls.md](../05_data_formats/06_merged_dataset_pitfalls.md)
- **預估解除**：下週 P0 master × 兩 flag × 7 樣本重跑完成
- **優先級**：🔴 高（已影響 2026-04-23 S5 PPT 圖製作）

### B1：Haplotag + ISM 全量重跑待執行
- **影響**：7 樣本 Phase 2A 無法啟動
- **依賴**：PON-only phasing 執行完成
- **預估解除**：2026-05（pending schedule）
- **優先級**：🔴 高

### B2：expected_coverage hardcoded 75.0 infra bug
- **現象**：master dataset 全 7 樣本共用 default
- **影響**：COLO829 (9x) 等低覆蓋樣本的 Coverage_Multiple 失真
- **解法**：`/cpp-change` 修 KDE 啟用並重跑
- **MEMORY**：`project_expected_coverage_baseline_bug`
- **優先級**：🔴 高

### B3：COLO829 TO 無 methylation
- **現象**：ONT_R10 data 缺 MM/ML tag
- **影響**：TO pipeline 無法跨 7 樣本完整驗證
- **解法**：Accept limitation 或 skip TO for COLO829
- **優先級**：🟡 中

### B4：HCC1395 chr8 熱點待分析
- **現象**：LOH+HPSig 7.4× FP enrichment
- **影響**：Zone-Aware characterization 需額外考量
- **MEMORY**：`project_hcc1395_chr8_hotspot`
- **優先級**：🟡 中

### B5：KDE Fix 下游 COLO829 衝擊
- **現象**：9 筆 variant 待重跑
- **MEMORY**：`project_kde_fix_downstream_quantification`
- **優先級**：🟡 中

---

## 風險清單

### R1：Phase 2 Normal Reference 跨樣本不一致
- **風險**：HCC1395 pilot POSITIVE 可能不複製到其他樣本
- **緩解**：7 樣本全量驗證（待 haplotag）
- **機率**：中

### R2：HPFineNGroups subclone marker 重新驗證失敗
- **風險**：flag=on 下 N≥3 消失，可能 somatic HP tag artifact
- **緩解**：Part B flag=on 重驗
- **機率**：中

### R3：TO 跨樣本 archive 需重跑
- **現象**：5 樣本 TO pipeline archive (post-HP-fix) 已彙整 master_extended 但缺 LOH/CN
- **建議**：Path 2 重跑 ISM（COLO829+DORADO+H2009）
- **MEMORY**：`project_to_cross_sample_archive_data_exists`
- **機率**：高（須執行）

---

## 緩解策略

| 阻塞 | 最短解決路徑 |
|------|-------------|
| B1 | 跑 PON-only phasing → 7 樣本 haplotag → ISM 批次 |
| B2 | `/cpp-change` 修 Config.hpp expected_coverage auto-estimate → rebuild → 重跑 7 樣本 |
| B3 | Decision：跳過 COLO829 TO 或 accept limitation |
| B4 | 獨立 research sub-project for chr8 |
| B5 | 等 KDE Fix 下游結果後重跑受影響 cycle |

---

## 相關

- Current Focus：[01_current_focus_snapshot.md](01_current_focus_snapshot.md)
- Milestones：[05_next_milestones.md](05_next_milestones.md)
- Coverage_Multiple：[../07_derived_features/02_coverage_multiple.md](../07_derived_features/02_coverage_multiple.md)
