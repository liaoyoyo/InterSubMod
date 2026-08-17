<!--
建立時間: 2026-08-13 09:52
目標: 公開文件 P0 修正過程 living document
處理範圍: 34 個 InterSubMod P0 claim families + CCU OLD-P0/REGRESSED findings
cycle_id: 20260813_public_docs_p0_correction
spec_id: public_docs_p0_correction
status: validated_local_source_remote_publication_pending
advisory: on
關聯檔案:
  - InterSubMod/research/20260813_public_docs_p0_correction/00_INDEX.md
  - InterSubMod/research/20260813_public_docs_p0_correction/pre-decision-audit.md
  - InterSubMod/.claude/skills/implementation-notes/SKILL.md
-->

# Implementation Notes：公開文件 P0 修正

> 本文件只記錄實作期設計決定、偏離、折衷、未決與 gotcha；最終逐 claim 結果另寫 correction receipt。

## 🔵 設計決定

### [2026-08-13 09:52] 以 frozen authority 統一公開詞彙

- **Status**：Accepted。
- **背景**：公開面把 molecule→cell、candidate topology→clone tree、association→causal confirmation 升格。
- **決定**：We will standardize on “local recurrence-allowed mutation-state candidate structure”, “physical-molecule co-occurrence”, and “pattern-conditioned methylation association”.
- **理由**：與 current authority machine contract 一致；代價是標題與圖說需較長限定語。
- **影響範圍**：README EN/ZH、Wiki、Pages、CCU OLD-P0／REGRESSED source。
- **Revisit if**：新增 L1 cellular truth、CN/LOH/purity/CCF calibration 且 authority manifest正式升版。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐。

### [2026-08-13 09:52] 分離 local correction 與 remote publication

- **Status**：Accepted。
- **背景**：GitHub About、main merge、Wiki publish、Pages deploy、CCU deploy 都是外部 state change。
- **決定**：We will report local source fixes and patch readiness separately from external publication; no push or deploy will be performed.
- **理由**：避免把 working-tree state 誤報成 live public state。
- **影響範圍**：correction receipt 的 disposition 欄與 final report。
- **Revisit if**：用戶明確授權 remote actions。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐。

### [2026-08-13 09:52] 保持原稽核 immutable，新增 correction cycle

- **Status**：Accepted。
- **背景**：2026-08-12 報告是當時 public snapshot 的證據，不應被事後改寫。
- **決定**：We will create new correction artifacts and link back to the immutable audit rather than edit its verdict/counts.
- **理由**：保留 before/after 可稽核性。
- **影響範圍**：`InterSubMod/research/20260813_public_docs_p0_correction/` 與新 validated report。
- **Revisit if**：僅在原報告有事實錯誤時以 erratum companion 處理。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐。

## 🟠 偏離之處

### [2026-08-13 10:16] 共享 branch 在執行中前進兩個既有 commit

- **Status**：Accepted after read-only conflict scan。
- **觀察**：`HEAD` 從本輪啟動時的 `83741469` 前進至 `95d420f6`；新增 commit 為
  `d3bcf0c8` 與 `95d420f6`，不是本輪 correction agents 建立。
- **衝突檢查**：`git diff 83741469..95d420f6 --` 對 README、Quickstart、Summary、
  `docs/wiki`、`docs/explain`、CMake、`include/`、`src/`、`tests/` 皆為空。
- **決定**：保留新 HEAD，不 reset；以 target-specific diff 與最終 HEAD 重寫 provenance。
- **風險評估**：低；兩 commit 涉及文件標準與資料資產追蹤，未覆蓋本輪 target 或 core build input。
- **Revisit if**：最終 conflict scan 發現上述 target 在本輪後又被修改。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐。

### [2026-08-13] Fresh-reader 阻擋首版報告

- **Status**：Resolved after revision。
- **觀察**：獨立讀者指出 P0 摘要表只列 33 個唯一 ID（漏 C136）、HTML 首屏易把 disposition／patch validated 誤讀成已發布修完、圖 1 SVG bar 互相重疊，且 dirty `HEAD` 不能識別 validated bytes。
- **決定**：補列 C136 與 34-ID 閉合式；全站統一 `PATCH_VALIDATED_ON_PINNED_CLONE / NOT_APPLIED / NOT_DEPLOYED`；修正 SVG 幾何；建立 33-entry correction SHA-256 manifest 與獨立 final-report sidecar；MD／HTML／footer 明示 uncommitted。
- **驗證**：correction manifest 33/33 OK；browser QA 最終 desktop 1440/1440、mobile 390/390、4/4 SVG a11y、exit 0；final reader gate 另存 receipt。
- **Evidence tier**：L1 ⭐⭐⭐⭐⭐。

### [2026-08-13 09:52] CCU source 以 clean temporary clone 修補

- **Status**：Accepted。
- **規範要求**：正式 source 通常在共享 workspace 內修改。
- **實作偏離**：CCU repo 不在 workspace；以 `/tmp/lab_tutorial_p0.GIEknc/repo` 的 clean `9eb1618` clone 修改並輸出 durable patch／receipt到 InterSubMod。
- **理由**：沒有 remote push/deploy authority，也不應任意在 `/big7_disk/liaoyoyo2001/` 建新 sibling repo。
- **風險評估**：temp clone 可消失；以 patch、commit SHA、checksums與 apply-check消除交付風險。
- **回退路徑**：maintainer 在正式 clone 套用 patch；若 upstream HEAD 變更，先三方 rebase/re-audit。
- **Revisit if**：用戶指定正式 CCU repo 工作路徑。
- **Evidence tier**：L2 ⭐⭐⭐⭐。

## 🟡 折衷考量

### [2026-08-13 09:52] 先完整 P0，不把 P1/P2 假裝解完

- **Status**：Accepted。
- **方案 A**：We will fix the complete predefined P0 scope and preserve a separate P1/P2 queue.
- **方案 B**：一次修改全部 58 個問題 claim + 32 CCU findings；風險是範圍過大、難以辨識回歸。
- **方案 C**：只修 README；不符合使用者要求與 Task B scope。
- **採用 A 理由**：P0 是完整且可驗收的 bounded unit；仍會在報告顯著標註 overall remediation 未完成。
- **Tier check**：34/34 P0 disposition，CCU OLD-P0 + REGRESSED 全 disposition；不以 P0 分母宣稱全站修復率。
- **Revisit if**：P0 全綠後用戶要求進入 P1/P2 cycle。
- **Evidence tier**：L2 ⭐⭐⭐⭐。

### [2026-08-13 09:52／終局 revisit] 直接 tracked HTML 與 generator-first 原則

- **Status**：Accepted；終局驗證已觸發 revisit 並完成。
- **方案 A**：We will locate and edit a generator when one exists; otherwise treat tracked standalone HTML as the source and protect it with residual tests.
- **方案 B**：一律直接全域替換 HTML；容易破壞 SVG、script與語意。
- **方案 C**：只加頁首警告、不修正文；無法消除錯誤 claim。
- **採用 A 理由**：兼顧維護性與現有 repo reality。
- **Tier check**：每頁 diff review + DOM/browser QA。
- **Revisit 結果**：找到 page-07 canonical generator；已同步 `CLEAN`／CCF／subclone 語意、補 SVG title/desc、實際重建，並把 generator 納入 C125/C127/C128 claim guard。
- **Evidence tier**：L3 ⭐⭐⭐。

## 🔴 未決問題

### [2026-08-13 09:52] GitHub About 與遠端發布由誰執行

- **Status**：open。
- **Question**：P0 local-ready 後，由 repo owner 直接套用建議 About 文案並 merge/publish，或另授權本 agent？
- **Context**：本輪沒有外部寫入授權。
- **Options**：owner 手動；後續明確授權 agent；只保留 patch。
- **Default if no answer**：標 `external_action_required`，不發布。
- **Revisit if**：本機 QA 全綠且使用者要求上線。
- **Priority**：critical。
- **Evidence tier**：L5。

## 📚 Lore

### [2026-08-13] Generated HTML 與 generator 必須同時受 claim guard

- **Constraint**：只修 tracked standalone HTML 會在下次重建時復發。
- **Resolution**：page-07 generator 與 generated HTML 同步；終局 guard 為 27 documents／77 target rules／140 required／79 forbidden，errors=0。
- **Evidence**：`InterSubMod/research/20260813_public_docs_p0_correction/claim_guard_receipt.md`。

### [2026-08-13] 66.52% 的 grain 是 strict components

- **Constraint**：`170,131/255,752=66.5219%`，不是 `170,131/469,849` 的 all-mutation rate。
- **Why it matters**：同一精確數字換錯母體會形成看似合理但科學錯誤的 claim。
- **Evidence**：`InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv:2-5`。

### [2026-08-13] “read = cell” 是本輪最高風險抽象錯誤

- **Constraint**：read/alignment record 提供受測序、比對、phasing 與 tag-calling error 影響的 molecular recorded/called evidence；originating cell 未 barcode、未知，不能把 HP/PS 或 MM/ML 當無誤差物理真值。
- **Why it matters**：所有 clone、lineage、cell proportion 與 correct cell tag claim 都會被這個錯誤放大。
- **Evidence**：`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md:39-49`。

## Provenance Footer

- Commit at start：`83741469`；base HEAD：`95d420f6`；本輪 deliverables 未 commit，不能以 HEAD 識別
- Validated correction bytes：33-entry manifest `4c7b0e6004c565dd55070793231dc280b8755edc5965fb8a0b1d9df05ecaf6e1`
- Final report bytes：`InterSubMod/research/20260813_public_docs_p0_correction/final_report_manifest.sha256`
- Branch：`chore/handoff-20260813`
- Build time：2026-08-13 09:52 +08:00
- Cycle：`20260813_public_docs_p0_correction`
- Line count：<400。
