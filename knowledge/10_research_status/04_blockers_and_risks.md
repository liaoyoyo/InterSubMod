---
id: ism-kb-10-research-status-blockers-and-risks
name: "Blockers & Risks"
description: "當前阻塞：T1.2 V6 production tag Hard Gate、Tier 2 證據強化 backlog、~63% FP framework gap、CURRENT_FOCUS goal drift；風險清單與緩解。⚠ 2 週有效。"
status: active
last_verified: 2026-05-18
content_nature: runtime-fact
doc_type: reference
verified_scope: "blockers against docs/CURRENT_FOCUS.md §2026-05-17 Tier 1-4 + MEMORY pending items"
related_ids:
  - ism-kb-10-research-status-index
  - ism-kb-10-research-status-current-focus-snapshot
  - ism-kb-10-research-status-next-milestones
  - ism-kb-07-derived-features-coverage-multiple
  - ism-kb-05-data-formats-merged-dataset-pitfalls
tags: [status, blockers, risks, hard-gate, tier-1-4, thread-d, v6-production]
canonical_paths: [10_research_status/04_blockers_and_risks.md]
alias_paths: []
---

# Blockers & Risks

> ⚠️ **此為 2026-05-18 快照（基於 2026-05-17 plan tender-pondering-blossom）**
> 決策前先核對 `docs/CURRENT_FOCUS.md §2026-05-17`
>
> **本快照已從 2026-04-22 B1/B2/B3 鏡像深度更新到 Tier 1-4 雙軌（Thread D paper + V6 production）阻塞階層。**

- 一句結論：當前主阻塞為 T1.2 V6 production tag Hard Gate 5-day workflow + Tier 2 證據強化 backlog（Z-AUTO KDE 4-sample + Strategy A 章節骨架）；風險集中於 framework 63% FP gap + cancer-only limitation reviewer 質疑機率
- 適用對象：決策前、阻塞解除優先級排序、Hard Gate 評估前
- 可直接執行命令：
  ```bash
  # Tier 1-4 阻塞狀態
  sed -n '15,80p' /big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md
  # Hard Gate workflow
  grep -A 10 "Tier 1.2 V6 Production Tag" /big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md
  ```

---

## §1 當前主要阻塞（按 Hard Gate / 工作流序列化）

### B-T1.2：V6 Production Tag Hard Gate（🔴 最高優先）
- **Track**：V6 production
- **影響**：解鎖 Thread D paper Tier 2 Archive TO 7-sample rerun + T4.3 PI errata package
- **依賴**：5-day workflow（Day 1-2 COLO829 V6 ISM 補完 → Day 3 7-sample marker → Day 4 manifest + git tag → Day 5 PI errata email）
- **Gate Level**：Day 4 `git tag v6-prod-{YYYYMMDD}` + Day 5 email send 兩個 🔴 Hard Gate
- **預估解除**：W3 結束（2026-05-22）
- **計劃**：`InterSubMod/research/selfphasing_v6_production/00_PLAN.md`

### B-T2.1：Z-AUTO KDE 跨 4 樣本擴展待執行（🟠 高）
- **Track**：Thread D paper
- **影響**：⭐3 → ⭐4 升級必要條件；未通過則 Strategy A §3 HCC1395 primary 章節缺 cross-sample 證據
- **依賴**：B-T1.2 完成（V6 binary tag）→ 4 樣本 (H1437/H2009/HCC1954/HCC1937) KDE 各自做
- **預估解除**：W4 中

### B-T2.2/T2.3：Strategy A 章節骨架待寫（🟠 高）
- **Track**：Thread D paper
- **影響**：Tier 3 paper draft 啟動 prerequisite
- **依賴**：T2.1 部分證據 + HCC1395 primary discovery + 6-sample replication cohort
- **預估解除**：W3-W4 並行

### B-FW-GAP：Framework 63% FP unexplained（🟡 中，conceptual）
- **影響**：Reviewer 可能質疑「framework 不完整」
- **緩解策略**：
  - 接受 limitation（決策 #7）+ Discussion 章節 (T3.5) 明示
  - Reactive T4.2 GC/mappability/repeat 新軸 pilot（reviewer 強質疑時觸發）
- **預估解除**：Paper submission 後 reviewer feedback 階段

### B0：Merged 合成檔資料陷阱（2026-04-23 持續）⚠
- **持續中**（未解除）
- **影響**：跨樣本 AF 視覺化、LOH-aware 分析全受影響
- **陷阱 1**：`merged_7samples_paired_full_plus_hcc1395_to.tsv.gz` 的 `AF` 欄位**非 caller_af**（非 HCC1395 的 5 樣本 p75 皆 <0.06）
- **陷阱 2**：HCC1395 `kde_status='phase1_new'` rows 的 LOH annotation 殘缺（Potential_LOH=True 僅 5.7% vs 正常 58.8%）
- **緩解**：改從 stale master archive 讀取（`all_region_rows.tsv.gz`），用 `caller_af` + `NumReads/bx` 重算
- **文件**：[../05_data_formats/06_merged_dataset_pitfalls.md](../05_data_formats/06_merged_dataset_pitfalls.md)

### B2：expected_coverage hardcoded 75.0 infra bug（持續）
- **現象**：master dataset 全 7 樣本共用 default
- **影響**：COLO829 (9x) 等低覆蓋樣本的 Coverage_Multiple 失真
- **解法**：`/cpp-change` 修 KDE 啟用並重跑（受 B-T1.2 影響需先確定 V6 production binary）
- **MEMORY**：`project_expected_coverage_baseline_bug`

### B3：COLO829 TO 無 methylation（接受 limitation）
- **現象**：ONT_R10 data 缺 MM/ML tag
- **影響**：TO pipeline 無法跨 7 樣本完整驗證
- **決策**：accept limitation（決策 #8 cancer-only acceptance 框架下）
- **後續**：Strategy A 樣本階層改為 HCC1395 primary + 6 replication（不含 COLO829 TO methylation 軸）

### B4：HCC1395 chr8 sample-specific hotspot 限制（characterization-only）
- **現象**：chr8 FP enrichment 2.31× highest of 23 chr，但 sample-specific 無法跨樣本 generalize
- **影響**：H4 POSITIVE 但無法升級為 cross-sample marker
- **MEMORY**：`project_hcc1395_chr8_hotspot`

---

## §2 已解除 / 解除中阻塞

| ID | 名稱 | 解除日期 | 證據 |
|----|------|---------|------|
| T1.1 | Thread D 主軸正名 banner 缺失 | 2026-05-16 | commit + 報告 338→381 行 |
| T1.3 | init-research scaffolding 缺失 | 2026-05-16 | `thread_d_paper/` + `selfphasing_v6_production/` 兩目錄 |
| Q5 | biorxiv/ensembl MCP「僵屍」誤判 | 2026-05-17 | commit `f3611a7` erratum + memory feedback |
| Goverance | CLAUDE.md/AGENTS.md 內容 drift | 2026-05-17 | v3 D2 分流 commit `696c7c1` |

---

## §3 風險清單（依 plan tender-pondering-blossom）

### R-PAPER-1：Cancer-only limitation reviewer 強質疑
- **風險**：審稿質疑 framework 對 non-cancer benchmark 適用性
- **機率**：中（接受 limitation 已標明於 Discussion）
- **緩解**：T4.4 HG002 non-cancer pilot reactive

### R-PAPER-2：Framework 37% coverage 不夠
- **風險**：審稿質疑 framework 不完整
- **機率**：中
- **緩解**：T3.5 Discussion 章節 + T4.2 reactive 新軸 pilot

### R-V6-1：V6 marker coverage 退步
- **風險**：V6 binary 比 V3F marker coverage 下降，導致 tag 不可發
- **機率**：低（pilot 顯示 +9% over V3F）
- **緩解**：Day 3 7-sample marker 檢查 → 通過才進 Day 4 tag

### R-V6-2：PI errata email 內容遺漏
- **風險**：T4.3 errata 5 條缺項
- **機率**：中
- **緩解**：Day 5 user review 強制 + PI Report 4-29 原文逐條對照

### R-KB-1：CURRENT_FOCUS goal drift > 25 天
- **風險**：KB snapshot 與 docs/CURRENT_FOCUS.md 內容偏離超過 2 週
- **機率**：高（hook 警告已觸發）
- **緩解**：本次（2026-05-18）深度更新 5 snapshot；之後 W4 重新評估

### R2：HPFineNGroups subclone marker 重新驗證失敗
- **持續**（未解除）
- **風險**：flag=on 下 N≥3 消失，可能 somatic HP tag artifact
- **緩解**：Part B flag=on 重驗

### R3：TO 跨樣本 archive Path 2 重跑
- **現象**：5 樣本 TO pipeline archive (post-HP-fix) 已彙整 master_extended 但缺 LOH/CN
- **建議**：依 V6 production tag 完成後 Path 2 重跑 ISM
- **MEMORY**：`project_to_cross_sample_archive_data_exists`

---

## §4 緩解策略矩陣

| 阻塞 | 最短解決路徑 |
|------|-------------|
| B-T1.2 | 5-day workflow 按 Day 1-5 順序執行；Day 4/Day 5 兩 Hard Gate 不可繞過 |
| B-T2.1 | B-T1.2 通過後並行 4 樣本 KDE（V6 binary） |
| B-T2.2/T2.3 | T2.1 部分證據 + HCC1395 primary 章節骨架先寫 |
| B-FW-GAP | T3.5 Discussion + T4.2 reactive |
| B0 | 改從 stale master archive 讀取 + `caller_af` 重算 |
| B2 | `/cpp-change` 修 Config.hpp expected_coverage auto-estimate → rebuild → 重跑（依賴 B-T1.2 V6 binary） |
| B3 | Accept limitation（決策 #8） |
| B4 | Independent research sub-project for chr8（characterization-only） |

---

## §5 Hard Gate 階層提醒（依 CLAUDE.md §1）

依專案 Hard Gate 不可逆操作清單：
1. **刪除檔案** 🔴
2. **C++ commit**（pre-commit hook 強制編譯檢查）🔴
3. **研究方向 NO-GO 判定**（依 Pre-registration 條件）🔴
4. **V6 production `git tag`**（不可逆）🔴 — T1.2 Day 4
5. **PI errata email send**（送出後不可逆）🔴 — T1.2 Day 5

所有 🔴 Hard Gate 觸發前必須用戶確認（非 auto-pass）。

---

## §6 相關

- 當前主軸：[01_current_focus_snapshot.md](01_current_focus_snapshot.md)
- 活躍假說：[02_active_hypotheses.md](02_active_hypotheses.md)
- Milestones：[05_next_milestones.md](05_next_milestones.md)
- Coverage_Multiple：[../07_derived_features/02_coverage_multiple.md](../07_derived_features/02_coverage_multiple.md)
- Merged dataset pitfalls：[../05_data_formats/06_merged_dataset_pitfalls.md](../05_data_formats/06_merged_dataset_pitfalls.md)
- /scientific-rigor §7 (Pre-registration 強制)：`InterSubMod/.claude/skills/scientific-rigor/SKILL.md`
- T1.2 Workflow：`InterSubMod/research/selfphasing_v6_production/00_PLAN.md`
