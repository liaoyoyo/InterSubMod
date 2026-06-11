---
name: known-pitfalls
description: InterSubMod 已知 AI 陷阱清單。每條記錄具體錯誤、正確做法、觸發場景。避免重複犯錯。USE WHEN「known pitfalls」「踩雷」「avoid mistake」「common bug」「之前怎麼錯的」「歷史教訓」、涉及 OLS/residualization、VCF 來源、特徵設計、AUC 分析、binary commit / KDE fix / working tree、證據鏈 / single-track / ⭐4-5 升級、**外部工具 F1/AUC/metric claim、論文 claim 對照、業界標準口徑（bcftools isec / hap.py / som.py）、tool README 行為敘述（longphase / clairs-to / claude code 內部演算法）— 對應 P-14 outside-claim 必查 KB**。SKIP WHEN 純 build / commit / docs 寫作、無對應陷阱類別的場景（先查表，無命中即略）、新 feature 探索初期（先設計再查陷阱）、純 UI / 視覺 / PPTX 製作（但 PPTX 涉及外部工具 metric 敘述時仍需查）。
allowed-tools: Read, Write, Edit
user-invocable: true
---

# InterSubMod Known Pitfalls（AI 已知陷阱）

> 每條記錄來自過去對話中 AI 犯的具體錯誤。新陷阱發現時追加到對應分類下。

## 使用規則

| 觸發場景 | 必讀陷阱 |
|----------|---------|
| 涉及 OLS / residualization / confound control | P-01, P-02 |
| 涉及 VCF 來源識別 / 數據溯源 | P-03, P-04, P-12 |
| 涉及特徵設計 / AUC 分析 | P-05, P-06, P-09, P-10 |
| 涉及證據鏈 / single-track / ⭐4-5 升級 | P-07 (anchor #1) |
| 涉及 binary commit / KDE fix / working tree | P-08, P-13 |
| 涉及跨樣本 n_passed / saturation | P-11 |
| **涉及外部工具 F1/AUC/metric claim、論文 claim 對照、業界口徑** | **P-14** ★ |

---

## 統計方法陷阱

### P-01: L2 Collider Bias

**錯誤**：對 near-constant 特徵（如 AlleleDelta in LOH regions）做 OLS residualization on AF，產生虛假 AUC 信號（表面 AUC 從 0.50 跳到 0.59）。

**正確做法**：L2 residualized AUC 必須用 L3 AF-bin 交叉驗證。若 L2 與 L3 差距 > 0.10，即為 collider bias，該特徵應判定 CONFOUND。

**來源**：O12 LOH 甲基化場景分析。Memory: `feedback_L2_collider_bias.md`

### P-02: Pooled OLS Residualization Trap

**錯誤**：Pooled OLS（TP+FP 合併後 fit）殘差仍保留分組信息，因為 TP/FP 在特徵空間中佔據不同位置，殘差 = 組間差 + 組內差。

**正確做法**：必須使用 within-group OLS（TP/FP 分別 fit），殘差才真正移除 confound。

**來源**：Beyond-AUC M2 驗證。Memory: `feedback_pooled_ols_residualization_trap.md`

---

## 數據來源陷阱

### P-03: VCF 來源錯誤歸因

**錯誤**：將 canonical VCF 錯誤歸因為 "chenhan112 pipeline"。實際上 canonical TO pipeline VCF 是 liaoyoyo2001 使用 ONT_5kHz BAM（有 5mCG+5hmCG MM/ML）執行 ClairS-TO 產生的。

**正確做法**：確認 VCF 來源時必須追蹤：(1) 誰執行了 caller，(2) 使用哪個 BAM，(3) BAM 是否有 MM/ML tags。查閱 Knowledge/02_samples/ 和 Knowledge/03_file_formats/ 交叉驗證。

**來源**：2026-04-14 TO pipeline staging v2 修正。

### P-04: pileup Symlink 指向錯誤 Caller

**錯誤**：pileup 模式的 output symlink 實際指向 ClairS paired（非 TO），導致 TO 分析使用了 paired caller 的輸出。

**正確做法**：追蹤 symlink 實際目標（`readlink -f`），確認 caller pipeline 與分析模式匹配（TO 分析必須用 ClairS-TO VCF，paired 分析用 ClairS paired VCF）。

**來源**：Memory: `project_vcf_source_error_correction.md`

---

## 特徵分析陷阱

### P-05: CramersV 93% Zero Artifact

**錯誤**：將 CramersV 視為連續區分特徵使用。實際上 CramersV 在 2×2 contingency table（ISM 的標準框架）中只有 {0, 1} 兩個值，93% 的 regions 為 0。

**正確做法**：CramersV 不適合作為連續特徵使用。使用 HPFineNGroups（已克服此限制，AUC 提升 +0.125）作為替代。

**來源**：R1-R5 特徵設計研究。Memory: `project_feature_design_limitations_r1r5.md`

### P-06: n_reads / NumReads Confound

**錯誤**：忽略 read count 對所有統計量的系統性影響。較多 reads → 較高統計功效 → PERMANOVA p-value 更小、HPFineNGroups 更大，但這反映檢測力而非生物效應。

**正確做法**：所有特徵分析必須控制 n_reads（residualize 或分層）。任何 AUC > 0.58 的特徵都需排除 read count confound 後才能判定。

**來源**：O11 heterogeneity 分析。Memory: `project_O11_heterogeneity_negative.md`

---

## 證據鏈陷阱（v1.6 anchor #1 hardening）

### P-07: Single-Track Validated Cycle (Missing Orthogonal Evidence)

**錯誤**：cycle 標 ⭐4 / ⭐5 但缺第四軌「Orthogonal」證據（archive comparison / replicate run / alternate caller）。L4 mandatory 要求 4-track：(i) Statistical (ii) Cross-sample (iii) Mechanism (iv) Orthogonal — 缺 (iv) 是最常見的撤回原因（如 04-26 thread B whitelist）。

**正確做法**：宣告 ⭐4 / ⭐5 之前，至少加 1 個 orthogonal-track artifact 並寫入 `plan.preconditions.upstream_reports`。否則降為 ⭐3（described, single-track）。

**來源**：plan v1.6 §4.5.4-G batch 5a anchor #1；validation-protocol L4 mandatory；Drill 1 retro hpfinengroups + thread_b 案例驗證

---

## 數據新鮮度與整合陷阱（v1.8 T1-5 從 2026-04 retract events 編寫）

### P-08: KDE-fix Stale Binary Downstream

**錯誤**：binary commit (KDE fix) 之後，下游 master_*.tsv.gz dataset 仍由舊版 binary 重建，bias 持續到下游分析（如 HCC1395 S3 TP 95.5%→58.3% post-fix）。

**正確做法**：使用 master dataset 前必驗 dataset 重建時間 ≥ binary fix commit time。`/check-staleness` 已自動檢查 `precheck.checks.binary.stale_distance`；≥1 時 BLOCK，要求重建。

**來源**：20260420_KDE_Fix_Acceptance_Validation。Memory: `project_kde_fix_downstream_quantification.md`

### P-09: Spatial Autocorrelation Confound (chr+pos aggregation)

**錯誤**：以 chr+pos 視窗聚合的特徵會帶入 spatial autocorrelation（linkage / hotspot density / replication timing），AUC 看似顯著但 mid-TP-rate window 分層後消失。

**正確做法**：所有 chr+pos 聚合特徵必須跑 mid-TP-rate window stratification；若 AUC 下降 >0.05 → 判定 spatial autocorr confound，要求 within-window control。

**來源**：Memory: `feedback_spatial_autocorrelation_confound.md`

### P-10: Feature Name Literal Interpretation

**錯誤**：把含生物學語意的特徵名（例如 HPFineNGroups → "methylation subclone marker"）當成生物學意涵直接論證。實際 C++ 實作可能只是 phasing-derived 統計量，與 methylation 無直接關係（HPFineNGroups 04-22 撤回事件）。

**正確做法**：宣告生物學語意前必讀 `src/include/` 對應 feature 定義；plan.preconditions 應加 `source_code_refs` 欄列出 C++ 路徑。

**來源**：04-22 HPFineNGroups 撤回事件。Memory: `feedback_feature_name_vs_definition_rule.md`

### P-11: Saturation Artifact (1/N Samples Drives All)

**錯誤**：cross-sample n_passed=1/N（單一樣本飽和如 H2009 thread B 04-26）但 pooled / mean aggregation 隱藏這個 1-vs-rest pattern。pooled metric 看似 promising，實則由 1 樣本拉動。

**正確做法**：若 `generalize.consistency.n_samples_passed / n_samples_total ≤ 0.2` → BLOCK 不論 pooled metric 多好。深入分析哪個樣本驅動信號（saturation? batch effect?）。

**來源**：20260426_Thread_B_Whitelist_Retraction

### P-12: Merged Dataset AF Schema Drift

**錯誤**：merged_*.tsv.gz dataset 的 `AF` 欄位語意可能不等於 caller_af（5/7 樣本 archive: AF p75 < 0.06）。下游用 `AF` 當 caller_af 處理→錯誤結論。

**正確做法**：用 merged dataset 前驗 schema：`AF.p75 > 0.10` AND `dataset_id contains "caller_af_separate"`。否則用 canonical/ 拉 caller_af join 替換。

**來源**：20260424_X6_Caller_AF_S3S5_CrossSample。Memory: `feedback_merged_dataset_af_and_loh_pitfalls.md`

### P-13: Working Tree Dirty / Uncommitted Binary

**錯誤**：pipeline 跑時 git working tree dirty（modified/uncommitted）。結果不可重現；effect size 通常 near-noise（如 longphase 04-29 ΔF1=-0.0003）卻被當訊號。

**正確做法**：`/check-staleness` 已加 `precheck.checks.git.working_tree_clean` 檢查。dirty → BLOCK 直到 commit OR `binary_version` 寫成顯式 SHA + diff snapshot 歸檔。

**來源**：20260429_longphase_TO_vs_V5_Somatic_Fallback

---

## 外部工具 / 論文 claim 陷阱

### P-14: outside-claim 必先查 KB（不可從本專案 report 推論外部行為）

**錯誤**：討論外部工具 (longphase / clairs-to / claude code) F1/AUC/metric claim 時，直接從本專案內部 report (§8.6.x 等) 推論外部行為，未先查 KB 對應 `Knowledge/05_tools/<tool>.md` 或論文 §N。

**症狀**：本專案 report 與外部工具 README/論文 claim 看似衝突（如「longphase-to 論文 §4.3 寫 F1 改善」vs「本專案 §8.6.2 寫 F1 三版相同」），實則**口徑不同**（V_H/V_L post-filter F1 vs caller-level FILTER=PASS F1），但因未查 KB → 浪費 3-4 輪來回才釐清。

**正確做法**：
1. 用戶問題或我即將回答時，若觸發關鍵字「longphase / clairs-to / claude code / 論文 / README / paper §N / 改善 F1 / metric claim / 業界標準」→ **先查**：
   - `mcp__knowledge__knowledge_search "<tool> F1"`
   - `Read /big8_disk/liaoyoyo2001/Knowledge/05_tools/<tool>.md`
2. 引用 KB 段落或 paper §N 再下結論
3. 若 KB 與本專案 report 衝突 → KB + 論文優先（除非本專案 report 明寫「我們的修補 vs 論文此處」）

**反例**：「§8.6.2 證據 ② 寫 F1 三版相同 → 結論 F1 不變」（沒查論文，沒釐清「哪個 F1 口徑」）

**正例**：「先 Read `Knowledge/05_tools/longphase-to.md` 確認論文 §4.3 F1 口徑為 V_H/V_L post-filter；本專案 §8.6.2 是 caller-level F1。兩者不衝突，是不同 metric。」

**來源**：2026-05-13 slide 14 F1 雙口徑事件。Memory: `feedback_outside_claim_must_query_kb.md`

---

## 數據誠信陷阱

### P-15: 捏造 metric / 報告搶先分析（同一 batch 平行）

**錯誤**：含數字報告把「預期值」當真值寫入；或報告 Write 與產生數字的分析 Bash 放**同一 tool-call batch** → Write 拿不到當批未回傳的數字 → 用記憶/預期補。

**症狀**：報告數字 grep 不到任何來源檔；或方向與真值相反（2026-06-01：H19 捏造 0.985；BRCA2 寫 0.572「弱」真值 0.866「強」）。同 session 30 分鐘違反 3 次 → 純文字規則必失敗，需機械層。

**正確做法**：CLAUDE.md §13.0 鐵則 — `分析跑完 → 落檔 → Read 驗證 → 才撰寫`；Write 與分析 Bash **永不放同一 batch**；機械後盾 §13 三層（`fill_report.py` by-construction / `number_provenance_check.sh` gate / `audit` 溯源表）。

**來源**：`20260601_fabricated_metric_in_html_preview_postmortem.md`；memory `feedback_no_fabricated_numbers_in_reports`。

### P-16: 修 Bug 不先 root-cause / 多變數齊改 / 失敗硬疊

**錯誤**：「先 quick fix 之後再查」；一次改多個變數無法隔離哪個有效；修 2-3 次不好仍硬加 fix 不質疑架構。

**正確做法**：`/cpp-change` 4-phase root-cause（穩定重現 → 多組件邊界 instrument → 單一假說一次一變數 → 失敗測試先行）；**2-3 次失敗 → 停下質疑架構**；回歸測試用 revert-must-fail-again 驗有效性。

**來源**：2026-06-02 借鑑 superpowers systematic-debugging。

### P-17: 盤點 / 狀態 / 索引類事實憑記憶寫 status（§13.0 非數字版）+ cross-branch 幻覺

**錯誤**：做「盤點 / audit / 現況整理」時，把某檔「存在 / current / stale / 有 N 個 run」當已知直接寫進報告，**沒實際 ls/grep 驗證**；或不標自己所在 git branch/worktree → 把「只在某 branch / 某 worktree 看到的檔」當全域 current。

**症狀**（2026-06-12 盤點稽核實例，連批判層都犯）：宣稱 HCC1395「7 個 complete_matrix run、版本混亂」實際 **1 個**；ISM region 路徑少一層（漏 `filtered_snv_tp/`）；`data_sources` 填充率「0%」實際 1 檔；批判反指「LAUNCH_READINESS 當前 branch 不存在」實際存在。**根因 = §13.0 同源：寫進去的事實沒回檔案 grep**，只是對象從「數字」換成「檔案/狀態/路徑」。

**正確做法**：
1. 任何 `status=current/stale` 或「某檔/路徑/數量」斷言，**旁邊必附「驗證指令 + 一行輸出摘要」**；無附 → 自動降 `status=unverified`，不可寫 current。
2. 盤點/audit/status 報告**首行標 provenance stamp**（branch + 短 HEAD + worktree path + 驗證日期）；助手：`bash scripts/provenance_stamp.sh`。
3. 跨 5 worktree 環境：宣稱「唯一 SoT」前先 `grep -Fxvf` 對帳，確認非各 branch 物理分叉（2026-06-12 實測 evidence_ledger MAIN=superset、無分叉，但 line 數 99/49/15/43/15 不同會誤導）。

**來源**：2026-06-12 數據準確度×可尋性盤點稽核（`docs/data_specs/20260612_data_accuracy_findability_improvement_audit_01.md`）；§13.7 完成宣稱 gate 延伸到盤點類斷言。

---

## 與 /scientific-rigor 元方法論的關係

本 skill 為 `/scientific-rigor` 元方法論層的**反例 / 陷阱知識庫**:
- `/scientific-rigor §4 DAG 因果審計` 引用 **P-01 (L2 Collider Bias)** + **P-02 (Pooled OLS Trap)** 作為「殘差 over collider → 自動降級 characterization only」的反例依據
- `/scientific-rigor §3 Effect Size + §7 文獻追溯` 引用 **P-14 (外部工具 claim 必查 KB)** 作為 effect size 宣告必跑 KB 對照的 trigger
- `/scientific-rigor §8.3 Reflexion buffer` 的 reopen threshold 設計呼應本 skill 的「陷阱可重啟條件」邏輯

**級聯觸發**: 任務啟動 → 本 skill 查歷史陷阱 → `/scientific-rigor §0.5 最小可用子集`決定下游 → 視場景進入 `§4 DAG` 或 `§9.2 Postmortem`

---

## Phase Chain Position & Dependencies

- **Phase**: 元方法論層（防重複犯錯，任何研究啟動前查詢）
- **Upstream**: 新研究方向 / 重大改動 / 用戶 query
- **Downstream**: 依 pitfalls 對齊到對應 skill 修正
- **Reads**: 本 skill 內 P-01 至 P-N pitfall 條目（references/）
- **Writes**: Inline 警告引用

## Failure Mode & Diagnostics

| 症狀 | 可能原因 | 修法 |
|------|---------|------|
| Pitfall 條目 >3 個月 stale | 新教訓未寫入 | 走 §9.2 postmortem 後同步條目 |
| 用戶忽略警告強推 | 缺 Hard Gate 對齊 | 依 §8.3.1 reopen threshold 拒絕，C1/C2/C3 一條方可繞 pitfall |

