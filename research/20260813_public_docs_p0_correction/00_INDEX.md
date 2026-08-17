<!--
建立時間: 2026-08-13 09:52
目標: 把 2026-08-12 公開文件全面稽核轉成可驗證的 P0 修正循環
處理範圍: Task B — 全部 34 個 InterSubMod P0 claim families，以及 CCU lab-tutorial 全部 OLD-P0 與 REGRESSED findings；P1/P2 保留在後續 queue
build_branch: chore/handoff-20260813
build_commit_at_start: 83741469
current_head_at_last_verify: 95d420f6
status: VALIDATED_LOCAL_SOURCE_WITH_EXTERNAL_PUBLICATION_REQUIRED
artifact_tracking: UNCOMMITTED_WORKSPACE_DELIVERABLE_WITH_SHA256_MANIFEST
worktree: /big7_disk/liaoyoyo2001/InterSubMod
data_sources: research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv,research/20260812_intersubmod_github_public_docs_full_validation/lab_tutorial_delta_reaudit.tsv,docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json
驗證方式: 逐 claim occurrence 對照、殘餘字串 gate、文件／HTML 結構檢查、桌機／手機／列印 browser QA；remote GitHub About、merge、Wiki publish 與 CCU deploy 明確列為外部 action
證據等級: L2 ⭐⭐⭐⭐（完整文件母體與 fresh runtime 已驗；未執行 remote publication）
關聯檔案:
  - InterSubMod/research/20260813_public_docs_p0_correction/pre-decision-audit.md
  - InterSubMod/research/20260813_public_docs_p0_correction/implementation-notes.md
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/00_INDEX.md
-->

# 2026-08-13 公開文件 P0 修正循環

> **Task B / P0 correction scope**：本輪完整處理既定 P0 母體，不把 P1/P2 或尚未發布的 remote 狀態算成已解決。

## 服務目標

- G2：避免把 molecular candidate structure 誤稱 cellular subclone。
- G3：把 methylation 限定為 pattern-conditioned association，不升格成 lineage/function truth。
- G4：同步 EN/ZH、Wiki、Pages 與 CCU 教材的版本與分母口徑。
- G5：建立可由 CI 重跑的 public-claim guard 與 correction receipt。

## Pre-registration（Confirmatory）

| 預測 H | 否證條件 | Decision threshold |
|---|---|---|
| H-P0-1：34 個 InterSubMod P0 claim families 都能映射到本機可控來源或明示 external action | 任一 P0 claim 無 occurrence／owner／修正文案，或把 About／remote merge 誤標成 local fixed | 34/34 有 disposition；local-fixed 與 external-action 分開計數 |
| H-P0-2：更正後所有本機公開入口都遵守 molecular-not-cell claim ceiling | target files 仍出現「confirmed clone／same cell lineage／methylation confirms clone／automatic heatmap」等禁語，或候選／local／model-conditional 限定缺失 | P0 residual gate 0 個未豁免命中；所有允許例外逐條列 receipt |
| H-P0-3：CCU OLD-P0 與 REGRESSED claims 可在 commit `9eb1618` 上形成可套用 patch | patch 無法乾淨套用、generated print-all 未同步、或 source page 仍保留 regression wording | 所有 OLD-P0 + REGRESSED finding 有 fixed／external disposition；patch apply-check PASS；site QA PASS |
| H-P0-4：修正不破壞既有公開文件與 HTML 互動 | link、TOC、details、responsive、print 或 no-JS 任一失敗 | 結構檢查與 Chromium desktop/mobile/print 全 PASS；0 console/page error；0 horizontal overflow |

## Step → Verify

1. 鎖定 dirty-tree 邊界與 audited versions → 驗證：列 branch／commit／target diff，非目標檔 0 次修改。
2. 逐項修 InterSubMod 34 個 P0 claim → 驗證：每個 ID 有 occurrence-before／after／evidence／disposition。
3. 修 CCU 全部 OLD-P0 與 REGRESSED → 驗證：基於 `9eb1618` 的 patch、generated aggregate 同步、finding-level receipt。
4. 新增 claim guard → 驗證：正向 contract、禁語、分母、CLI/schema assertions 全部自動檢核。
5. 重跑文件與 HTML QA → 驗證：命令 exit 0、實際輸出片段、browser viewport 與 print receipt。
6. 產生修正報告與 standalone HTML／SVG → 驗證：claim mapping、external actions、remaining P1/P2、provenance 皆可見。
7. 執行獨立 fresh-reader gate → 驗證：數量閉合、首屏狀態、SVG 幾何、未 commit provenance 與資料優先序可由只讀者重建。

## Scope boundary

- 本輪不修改生物資訊核心演算法、不重跑 7 樣本 BAM pipeline。
- 本輪不 push、merge、改 GitHub About、發布 Wiki、部署 Pages 或 CCU live site。
- P1/P2 不會被本輪的 P0 成功率吞掉；另列 remaining queue。
- 既有 2026-08-12 validated audit 保持 immutable；修正另建 receipt。

## 輸入與預期輸出

- 輸入：`InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv`
- 輸入：`InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/lab_tutorial_delta_reaudit.tsv`
- 輸入：`InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/{authority_manifest.json,denominator_registry.tsv}`
- InterSubMod 輸出：本機修正檔、claim receipt、validator 與 QA receipt。
- CCU 輸出：基於 `9eb1618d359e602d9c528675952b20d051fb2346` 的可套用 patch 與 QA receipt。
- 報告輸出：`InterSubMod/docs/reports/validated/2026/08/20260813_公開文件P0修正與驗證_01.{md,standalone.html}`。

## 終局輸出

- Markdown：`InterSubMod/docs/reports/validated/2026/08/20260813_公開文件P0修正與驗證_01.md`
- Standalone HTML／inline SVG：`InterSubMod/docs/reports/validated/2026/08/20260813_公開文件P0修正與驗證_01.standalone.html`
- Browser QA：`InterSubMod/research/20260813_public_docs_p0_correction/html_qa_receipt.md`
- Validated correction bytes：`InterSubMod/research/20260813_public_docs_p0_correction/validated_target_manifest.sha256`（33 entries；SHA-256 `4c7b0e...f6e1`）
- Final report bytes：`InterSubMod/research/20260813_public_docs_p0_correction/final_report_manifest.sha256`
- Fresh-reader gate：`InterSubMod/research/20260813_public_docs_p0_correction/reader_test_receipt.md`（第三輪 PASS）
- Evidence ledger：`InterSubMod/research/autoresearch/evidence_ledger.jsonl`（cycle `20260813_public_docs_p0_correction`）
- P0 guard：34/34 inventory、27 documents、77 target rules、140 required、79 forbidden、errors 0。
- External state：main／Wiki／Pages／About／CCU live 仍未發布；只允許標 `LOCAL_SOURCE_CORRECTED`／`PATCH_VALIDATED_ON_PINNED_CLONE`，CCU 正式狀態為 `NOT_APPLIED / NOT_DEPLOYED`。
