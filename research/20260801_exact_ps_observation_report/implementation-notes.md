<!--
建立時間: 2026-08-01 17:01 +0800
目標: Exact-PS 全資料最終 HTML observation report 實作過程 living document
處理範圍: Python builder、report-data schema、standalone HTML、receipt、pipeline final observation hook、browser QA
cycle_id: cycle_20260801-exact-ps-observation-report
spec_id: exact_ps_observation_report
status: validated
advisory: on
關聯檔案:
  - InterSubMod/research/20260801_exact_ps_observation_report/pre-decision-audit.md
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json
  - InterSubMod/.claude/skills/implementation-notes/SKILL.md
-->

# Implementation Notes: Exact-PS Final Observation HTML

> **Purpose**: 記錄從 frozen authorities 生成美觀、可重生且不越過 claim ceiling 的
> standalone HTML observation report。

## Frontmatter

- **Spec source**: 使用者 2026-08-01 要求「所有數據出來後，自動輸出美觀 HTML 觀察表」
- **AI session**: 2026-08-01
- **Last updated**: 2026-08-01 17:40 +0800
- **Line count**: <200 / 400

## 🔵 設計決定

### 2026-08-01 17:01 — Authority-driven derived view

- **Status**: Accepted
- **決定**: We will derive every displayed number from frozen JSON/YAML/TSV authorities
  identified by `authority_manifest.json`; HTML will never be a numerical SoT.
- **理由**: 防止漂亮頁面與 canonical output 漂移。
- **影響範圍**: builder、report data、receipt、HTML provenance drawer。
- **Revisit if**: 新 release manifest 取代現行 handoff manifest。
- **Evidence tier**: L1

### 2026-08-01 17:01 — Standalone progressive enhancement

- **Status**: Accepted
- **決定**: We will generate a single offline HTML with inline CSS/SVG and no external CDN;
  core tables and conclusions will exist before JavaScript runs.
- **理由**: 支援 workstation、交付、no-JS 與列印。
- **影響範圍**: HTML renderer、Playwright QA。
- **Revisit if**: 正式部署需要多人協作後端。
- **Evidence tier**: L2

### 2026-08-01 17:01 — CN/LOH absence is a visible state

- **Status**: Accepted
- **決定**: We will render CN/LOH as `NOT_INTEGRATED` with its required fields and next gate;
  we will not infer missing values.
- **理由**: 避免 AF winner 被誤讀為 CN-adjusted lineage。
- **影響範圍**: method-status matrix、tree-decision section。
- **Revisit if**: P0-2/P0-3 產生 hashed allele-specific CN source。
- **Evidence tier**: L1

### 2026-08-01 17:26 — Independent finalizer, not a frozen-runner mutation

- **Status**: Accepted
- **決定**: 使用獨立 post-run finalizer；不修改 2026-07-24 frozen cohort runner。
- **理由**: HTML 同時依賴 runner 後才完成的 exact topology census 與 methyl sidecar；
  runner 結束只是必要條件，不是完整 readiness。
- **影響範圍**:
  `InterSubMod/research/20260801_exact_ps_observation_report/scripts/finalize_exact_ps_observation_report.py`。
- **驗證**: finalizer 依序執行 builder、browser QA，最後才寫 release marker。
- **Evidence tier**: L1

### 2026-08-01 17:34 — Receipt 必須重新綁定 identity

- **Status**: Accepted
- **決定**: Finalizer 不信任單一 `all_pass=true`；必須驗 exact receipt schema、
  manifest／registry／13 authorities／9 strict nested bundles 與 report outputs SHA。
- **理由**: 紅隊證明舊版 gate 可由空 manifest 與偽造 receipt 通過。
- **影響範圍**: builder receipt、QA receipt、finalizer、negative tests、release reservation。
- **驗證**: forged build receipt 與 forged QA schema 測試均 fail closed。
- **Evidence tier**: L1

## 🟠 偏離之處

### 2026-08-01 17:34 — 增加 release reservation

- **偏離**: 原規格只規劃 staging→release；紅隊發現同 release ID 的 race condition。
- **處置**: 加入 atomic exclusive reservation file；失敗時與 staging 一起搬入
  `failed_attempts`，成功時永久保留 release ID reservation。
- **結果**: 第二個同 ID 執行不會移動或覆寫第一個 release。

## 🟡 折衷考量

### 2026-08-01 17:01 — 新 builder 與舊工作站

- **Status**: Accepted
- **方案 A**: We will build a small authority-driven report builder and reuse only verified
  visual patterns from prior workstations。
- **方案 B**: 修改大型 layered workstation；拒絕，耦合高且資料契約不同。
- **方案 C**: 手工寫一份 HTML；拒絕，無法在新數據後重生。
- **採用 A 理由**: 最小耦合、可測試、可作 final observation hook。
- **Revisit if**: 現有 builder 已完全涵蓋新 handoff decision policy。
- **Evidence tier**: L2

## 🔴 未決問題

### 2026-08-01 17:01 — Pipeline hook 位置

- **Status**: resolved
- **Question**: 應修改 20260724 research runner，或新增 stable wrapper／post-run command？
- **Context**: 舊 runner 對 frozen 研究結果有 authority 意義，直接修改可能破壞重現性。
- **Resolution**: 新增獨立 finalizer CLI，不改 frozen runner；正式 `all7_v1` 已由此入口發布。
- **Revisit if**: 盤點發現既有 runner 已有正式 report hook。
- **Priority**: critical
- **Evidence tier**: L5

## ✅ 驗證與正式輸出

### 2026-08-01 17:40 — Contract tests

```text
Ran 12 tests in 9.735s
OK
```

- authority artifact SHA drift：阻擋。
- duplicate artifact ID：阻擋。
- denominator row loss：阻擋。
- non-empty output／existing release：阻擋。
- forged build／QA receipts：阻擋。
- HTML／report-data 數值不一致與 schema drift：阻擋。
- QA failure：保留 failed attempt，不發布 release marker。

### 2026-08-01 17:40 — Formal all7 release

- **Input manifest**:
  `InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json`
- **Output**:
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260801_exact_ps_observation_report/all7_v1/`
- **Browser QA**: 1440／1024／390／320、no-JS、A4 全 PASS。
- **Integrity**: report data、HTML、build receipt、release receipt sidecars 全 PASS。
- **Release status**: `VALIDATED_DERIVED_OBSERVATION`。
- **Scientific status**: `INCOMPLETE_WITH_ABSTAIN`，因 10,717 resource abstains
  且 CN/LOH 尚未整合。

### 2026-08-01 17:45 — Report-data schema

- **Path**:
  `InterSubMod/research/20260801_exact_ps_observation_report/schemas/report_data.schema.json`
- **Contract**: draft-07、22 個頂層必填欄位、7 個 exact samples、23 個
  named conservation checks。
- **Positive validation**: 正式 `all7_v1/report_data.json` exit 0。
- **Negative validation**: duplicate/missing sample、CN status 升格、methyl topology
  rescoring、tree-decision materialization 升格與 check 缺列均被拒絕。

## 📚 Lore

### 2026-08-01 — Ranked-only denominator

- **Constraint**: 63,506/71,955 是 ranked-only rooted-unlabeled topology resolution。
- **Why it matters**: 不可畫成 63,506/98,955 或稱全區域 biological topology。
- **Evidence**: handoff `denominator_registry.tsv`。

### 2026-08-01 — Local blocks are not a whole-sample tree

- **Constraint**: 不同 PS×HP blocks 沒有 bridging evidence／union re-solve 時不可串樹。
- **Why it matters**: HTML 可列多個 local candidate trees，但不可畫假的全樣本 lineage。
- **Evidence**: handoff read-AF/CNV/single-tree decision spec。

## Provenance

- Commit: `387a101e6a3292e0d7f230ba8d20271c7434972a`
- Skill: `/implementation-notes` v0.1
- Pre-decision audit: `InterSubMod/research/20260801_exact_ps_observation_report/pre-decision-audit.md`
