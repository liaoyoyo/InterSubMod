<!--
建立時間: 2026-08-13 10:40
目標: 完整研究交接與 GitHub 公開準備的 process-time living record
處理範圍: InterSubMod handoff snapshot + LongLineage private public-preview preparation
cycle_id: cycle_20260813-1032-complete-research-handoff
spec_id: complete_research_handoff
status: in_progress
advisory: on
build_branch: agent/research-handoff-audit-evidence-20260813
build_commit: ddd8909a838318d8a77969313e9561c8ff9d01c2
worktree: /big7_disk/liaoyoyo2001/worktrees/ism-handoff-20260813
data_sources: InterSubMod/research/20260813_complete_research_handoff/pre-decision-audit.md
驗證方式: 每個 Accepted/Closed entry 都需對應 receipt 或 commit；否則維持 Proposed/open
證據等級: L2 ⭐⭐⭐⭐（process record；不本身授權 release）
關聯檔案:
  - InterSubMod/research/20260813_complete_research_handoff/00_INDEX.md
  - InterSubMod/research/20260813_complete_research_handoff/pre-decision-audit.md
  - InterSubMod/.claude/skills/implementation-notes/SKILL.md
-->

# Implementation Notes：完整研究交接與 GitHub 公開準備

> **Purpose**：記錄規格解讀、偏離、折衷、未決與不明顯約束；本檔不是 release receipt。

## Frontmatter

- **Spec source**: 使用者 2026-08-13 鎖定之完整實作計畫
- **AI session**: 2026-08-13 Codex
- **Last updated**: 2026-08-13 13:08 +08:00
- **Cycle**: `cycle_20260813-1032-complete-research-handoff`

## 設計決定（Design Decisions）

### [2026-08-13 10:20] 以兩個獨立 clean worktree 隔離 release

- **Status**: Accepted
- **背景**: 原 InterSubMod workspace 有 tracked/untracked changes；LongLineage 本機 HEAD 晚於 candidate。
- **決定**: We will implement only in dedicated worktrees rooted at the two user-locked commits.
- **理由**: 不 stage/reset 原研究工作；release diff 可相對固定 baseline 審核。
- **影響範圍**: 所有 commit、test、receipt 與 PR。
- **Revisit if**: 使用者明確變更 baseline。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐（HEAD/dirty checks 可重現）

<!-- BEGIN USER-SPECIFIED -->

### [2026-08-13] 科學與發布 ceiling

- **Status**: Accepted / Protected
- **Decision**: authority/denominator 不改寫；0 confirmed cellular subclone、0 linear ancestry、association-only、CN/LOH not integrated；88.2579% 不得升格。
- **DO NOT change**: 不以文件包裝、clean build 或 preview CI 提高 biological claim。
- **Rationale**: 使用者固定邊界 + frozen authority。

### [2026-08-13] Exclusion lock

- **Status**: Accepted / Protected
- **Decision**: `73afaeac`、08-13 drilldown/CNV active research、LongLineage `9ad976b`／`6ce62b2` 與其他後續 dirty work 不進 release baseline。
- **DO NOT change**: 只可登錄 `IN_PROGRESS/PARTIAL`，不可 cherry-pick 進 snapshot。
- **Rationale**: 防止 active science 與 handoff release 混線。

<!-- END USER-SPECIFIED -->

### [2026-08-13 10:40] 分離 frozen replay 與 source reproduction

- **Status**: Accepted
- **背景**: 19/19 hash replay 只證明登錄路徑未漂移；不證明 `ddd8909a` 重生該 science。
- **決定**: We will maintain separate fields for authority artifact hash, producer/source snapshot, baseline build, and regeneration status.
- **理由**: 防止 clean build 被誤讀為 science rerun。
- **影響範圍**: artifact registry、acceptance receipt、README、HTML。
- **Revisit if**: 有 byte-reproducible clean-source rebuild receipt。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-08-13 13:08] 可攜介面採 fail-closed site profile 1.1

- **Status**: Accepted
- **背景**: 只驗證檔案存在無法證明 genome build、index、工具 bytes 或 pipeline clone 一致。
- **決定**: `site-profile` 1.1 強制 reference contig contract、必要 index、required-tool SHA-256、project-root binding 與 real-data preflight；執行時鎖定 profile bytes，不讓 child step 重讀可變設定。
- **理由**: 將 bip7／bip8 的路徑差異轉成可驗證 contract，並阻止 inherited environment 或 `latest` 猜測偷偷改變輸入。
- **影響範圍**: `scripts/site/*`、portable pipeline、machine receipt 與 CI。
- **驗證**: 222 Python tests + 8 subtests、7/7 JSON Schema、13/13 handoff package checks；真正雙機 receipt 仍為 publication blocker。
- **Revisit if**: 需要支援第二 genome build，或 profile schema 進入 2.x migration。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

## 偏離之處（Deviations）

### [2026-08-13 10:40] pre-decision verdict 為 PROBE 仍開啟 living record

- **Status**: Accepted
- **規範要求**: implementation-notes 通常接在 verdict GO 後。
- **實作偏離**: 在 scoped GO／publication NO-GO 的 PROBE 狀態啟動 notes。
- **理由**: 使用者明示實作；可逆的 registry、docs、tests 與 private PR 可安全推進，publishing action 已鎖。
- **風險評估**: 若讀者忽略 gate，可能把工作進度誤當 release readiness。
- **回退路徑**: 保持本檔 in_progress；任何公開 action 由 acceptance receipt 另行授權。
- **Revisit if**: 全部 blocking gates PASS 或任一核心 falsifier 觸發。
- **Evidence tier**: L2 ⭐⭐⭐⭐

### [2026-08-13 10:40] 本輪不重跑 TiB science

- **Status**: Accepted
- **規範要求**: Comprehensive validation 通常要求 full data rerun。
- **實作偏離**: full data roots 全盤點，但 science 以 frozen hash replay 驗證。
- **理由**: 使用者明示 scope；現有 frozen authority 是本輪輸入，不是待重算輸出。
- **風險評估**: 只能證明 frozen snapshot integrity，不證明新 commit science parity。
- **回退路徑**: 後續另開 compute cycle，絕不在 handoff receipt 補寫成已跑。
- **Revisit if**: authority mismatch、source changes 或新 science claim。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

## 折衷考量（Trade-offs）

### [2026-08-13 10:40] 機器可讀 registry 優先於單篇長文

- **Status**: Accepted
- **方案 A**: We will make registries/schemas/receipts the SoT and render Markdown/HTML from them.
- **方案 B**: 手寫完整 narrative；rejected，容易 drift 且 AI 無法精確查詢。
- **方案 C**: 只保留路徑清單；deferred，缺 semantic/finality/provenance。
- **採用 A 理由**: 同時服務新人、AI、CI 與 release audit。
- **Tier check**: 工程 claim；不涉及 biological effect size。
- **Revisit if**: registry validator 無法涵蓋必要語意。
- **Evidence tier**: L2 ⭐⭐⭐⭐

### [2026-08-13 10:40] 實際路徑與 alias 分 repo 治理

- **Status**: Accepted
- **方案 A**: InterSubMod handoff 記實際 bip7/bip8 paths；LongLineage public docs 只用 aliases。
- **方案 B**: 兩 repo 都公開絕對路徑；rejected，違反 LongLineage public-safety governance。
- **方案 C**: 兩 repo 都只用 alias；rejected，降低本機交接可用性。
- **採用 A 理由**: 符合使用者固定政策並保持 public safety。
- **Revisit if**: machine path 被判為敏感資訊。
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

## 未決問題（Open Questions）

### [2026-08-13 10:38] LongLineage visibility incident（P0）

- GitHub API 在本輪實測回傳 `liaoyoyo/LongLineage` 為 `PUBLIC`，與 private-first
  governance 衝突。
- 已執行 containment 並回讀為 `PRIVATE`；receipt：
  `InterSubMod/research/20260813_complete_research_handoff/receipts/longlineage_visibility_containment_20260813.json`。
- containment 不能證明過去無 clone/download；在 source-origin、license、notices、SBOM、
  history-safety 與 P8 全數通過前，維持 `FAIL_CLOSED`。

### [2026-08-13 10:40] 真正 bip8 session 是否可取得

- **Status**: open
- **Question**: 能否在 hostname=`bip8` 的獨立 session 執行 fresh-clone receipt？
- **Context**: bip7 可見的 `/bip8_disk` NFS 不等於 bip8 主機驗收。
- **Options**: 本輪取得／交由 bip8 執行者補 receipt／保持 blocked。
- **Default if no answer**: `BIP8_INDEPENDENT_RECEIPT_BLOCKED`。
- **Revisit if**: 取得 SSH/session 或外部簽署 receipt。
- **Priority**: critical
- **Evidence tier**: L5 ⚠

### [2026-08-13 10:40] GitHub admin gate 可否在 PR 合併前配置

- **Status**: open
- **Question**: token 是否有 branch-protection 與 visibility 管理權限？
- **Context**: 技術準備可先完成，但不應在 checks 尚未存在時鎖錯 required context。
- **Default if no answer**: 只建立 draft PR 與配置說明，不改 main protection/visibility。
- **Revisit if**: workflow 合併且 admin scope 經 read-only check 確認。
- **Priority**: major
- **Evidence tier**: L5 ⚠

### [2026-08-13 10:40] LongLineage 10 筆 NOT_VERIFIED source mapping 的 legal disposition

- **Status**: open
- **Question**: 每筆來源是否可證明 original／clean-room／compatible derivative？
- **Context**: 無 disposition 不得切 public。
- **Default if no answer**: repo 維持 private、preview PR 可繼續。
- **Revisit if**: mapping reviewer 簽署與 notices/SBOM 完成。
- **Priority**: critical
- **Evidence tier**: L5 ⚠

## Lore（Prior Gotchas / Non-obvious Constraints）

### [2026-08-13] Git worktree checkout 可被短輸出視窗中斷

- **Constraint**: InterSubMod 初次 worktree add 在 41% 中斷，留下空 index＋部分 files。
- **Why it matters**: `HEAD` 正確不代表 checkout clean；必查 tracked、missing、dirty。
- **Evidence**: 修復後 `TRACKED=6953`、`MISSING_TRACKED=0`、`DIRTY=0`。

### [2026-08-13] 08-12 claim inventory 被 ignore

- **Constraint**: 158-row TSV 在來源 workspace 有 checksum，但未被 Git tracking。
- **Why it matters**: 報告引用存在不等於 fresh clone 可取得 evidence matrix。
- **Evidence**: `git check-ignore` 命中；SHA-256=`d5d8794c…642ddbc`。

### [2026-08-13] 目前只驗到 bip7 加 NFS data roots

- **Constraint**: hostname=`bip7`；`/big8_disk`、`/bip8_disk` 為 NFS mounts。
- **Why it matters**: 雙 mount 不是雙主機驗收。
- **Evidence**: red-team hostname/findmnt receipt。

### [2026-08-13] 無副檔名 CLI 必須依 shebang 分類

- **Constraint**: `scripts/site/bootstrap` 與 `scripts/site/doctor` 是 Python，不是 shell。
- **Why it matters**: 依檔名猜語言會讓 `bash -n` 產生假失敗；registry 與驗證命令均應讀 shebang。
- **Evidence**: Python `py_compile` PASS；shell-only 檔案 `bash -n` PASS。

## Provenance Footer

- **Commit**: `ddd8909a838318d8a77969313e9561c8ff9d01c2`
- **Built**: 2026-08-13 10:40 +08:00
- **Skill**: `/implementation-notes` v0.1
- **Pre-decision**: `InterSubMod/research/20260813_complete_research_handoff/pre-decision-audit.md`
- **Publication state**: blocked until acceptance receipt explicitly passes every release gate.
