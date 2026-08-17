<!--
建立時間: 2026-08-13 09:52
目標: 公開文件 P0 修正循環開始前的證據與風險審計
處理範圍: 34 個 InterSubMod P0 claim families + CCU OLD-P0/REGRESSED findings
cycle_id: 20260813_public_docs_p0_correction
topic: public_docs_p0_correction
status: verdict_GO
audit_version: 0.1
關聯檔案:
  - InterSubMod/.claude/skills/pre-decision-audit/SKILL.md
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/00_INDEX.md
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv
-->

# Pre-Decision Audit：公開文件 P0 修正

> **Verdict：GO（70/100）**。既有全面稽核已提供固定母體、逐 claim evidence 與最小修正文案；本輪是可逆的文件與驗證工程，不是重新推導 cellular lineage。

## §0 Cynefin Domain Gate

- **Domain**：Complicated。
- **Test**：相同「claim inventory → bounded rewrite → residual gate → browser QA」流程已有 2026-08-11／12 可重現結果。
- **Rationale**：科學 boundary 複雜，但本輪決策已由 frozen authority 固定；主要風險是跨版本面同步與 generated artifact 漂移，可用專家規則和測試管理。

## §1 Observation Completeness Checklist

| Observation | 狀態 | Tier | 來源 |
|---|---|---|---|
| 158 個 InterSubMod claim families 已完整盤點，34 個為 P0 | ✓ | L2 ⭐⭐⭐⭐ | `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv` |
| CCU 32 個 problem-focused delta findings 已重驗；29 未完全解決 | ✓ | L2 ⭐⭐⭐⭐ | `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/lab_tutorial_delta_reaudit.tsv` |
| Current authority 明禁 confirmed cellular clone／lineage／CN-LOH-corrected CCF | ✓ | L1 ⭐⭐⭐⭐⭐ | `InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json:78` |
| Exact denominator grain 已 machine-readable 鎖定 | ✓ | L1 ⭐⭐⭐⭐⭐ | `InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv:2` |
| CCU current remote 仍是 audited commit `9eb1618` | ✓ | L2 ⭐⭐⭐⭐ | `git -C /tmp/lab_tutorial_p0.GIEknc/repo rev-parse HEAD` → `9eb1618d...` |
| GitHub About／remote merge／deploy 可由本機檔案直接完成 | ✗ | L1 ⭐⭐⭐⭐⭐ | 外部平台 state；本輪未獲 push／publish 授權 |

## §2 Credibility Score

| 維度 | 評分 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | claim boundary 與 denominator contract 已 formalized |
| 觀察支撐 | 20 | 全公開母體、fresh runtime、current CCU commit 均已盤點 |
| 機制清晰度 | 20 | 核心錯誤鏈是 molecule→cell、candidate→truth、association→causality 的層級提升 |
| 反例風險 | 10 | dirty tree、generated mirrors、remote state 與 wording false-positive 仍需防護；依 conflict scan 降一階 |
| 所需資源 | 0 | 全 P0、多 repo、browser QA 預期超過 6 小時 |
| **TOTAL** | **70 / 100** | **GO** |

**Falsifier observable**：若修正假說錯，claim guard 或獨立 reader 仍能在 target surfaces 找到未豁免的 cellular-clone／因果／錯分母／錯 CLI claim，或 HTML／patch validation 失敗。

**Reality-test 三個反例觀察**：

1. 文字改成「candidate」但圖中的 node labels／caption 仍叫 clone/cell，代表 semantic fix 不完整。
2. 本機 README/Wiki/Pages 全綠，但 GitHub About／main／Wiki remote 仍舊，代表只能報 local-ready，不能報 public-fixed。
3. CCU source 改對但 `print-all.html` 或 data JS 仍舊，代表 generated aggregate 漂移。

## §3 Assumption Map

| Assumption | Importance | Known | Quadrant |
|---|---|---|---|
| 2026-08-01 authority 是本輪 scientific ceiling | HIGH | KNOWN | (1) |
| 34 個 P0 inventory rows 無遺漏核心 public claim family | HIGH | KNOWN within audited version | (1) |
| 現在本機 public target files 未被其他 dirty work 改動 | HIGH | KNOWN：target diff 為空 | (1) |
| CCU `9eb1618` 可代表 current remote | HIGH | KNOWN：fresh clone HEAD match | (1) |
| GitHub About／remote default branch 可由 local edit 同步 | HIGH | UNKNOWN／false without external action | (2) ⚠ |
| P1/P2 可延後而不影響 P0 correctness | LOW | KNOWN，只影響完整 remediation，不影響本輪 P0 定義 | (3) |

**MUST validate first**：remote state 必以 `external_action_required` 單列，禁止納入 local fixed numerator。

## §4 Quick Pilot

GO 已達，不需另跑 PROBE。第一個 checkpoint 仍採最小 reversible pilot：先修 README EN/ZH 一個 claim family並跑雙語一致性／禁語檢查；PASS 才擴到 Wiki／Pages／CCU。

## §5 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| GitHub About edit authority | HIGH：live landing metadata 仍會舊 | 0.1 h（需 owner 操作） | P0 external |
| Remote merge／Wiki publish／Pages deploy authority | HIGH：local-ready 不等於 live-fixed | 0.5–1 h + review | P0 external |
| CCU repo push/deploy authority | HIGH：本輪只能交 patch | 0.5–1 h + maintainer review | P0 external |
| portable public fixture | HIGH：quickstart 仍 internal-only | 2–6 h | P1 follow-up |
| cellular truth + allele-specific CN/LOH/purity/CCF | HIGH：限制 biological claim ceiling | >100 h / new data | research dependency |

## §6 Evidence Conflict Scan

- `MEMORY.md`：repo 目前沒有此檔，不能假裝已查；以 append-only `InterSubMod/research/autoresearch/evidence_ledger.jsonl` 與 validated reports 替代。
- `find InterSubMod/docs/reports/validated -iname '*NEGATIVE*'`：找到 `InterSubMod/docs/reports/validated/2026/04/20260401_LOH_weekly_review/06_methylation_hypothesis_negative.md`。

| Prior conclusion | Tier | Relation | Source |
|---|---|---|---|
| Public audit：34 P0、CCU 29/32 unresolved | L2 | support；直接定義修正母體 | `InterSubMod/research/autoresearch/evidence_ledger.jsonl:145` |
| Cellular clone/lineage confirmed = 0 under current authority | L1 contract | support；禁止文字升格 | `InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json:78` |
| Methylation association ladder不支持 clone/lineage | L1 contract | conflict with existing public wording | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md:39` |
| 多條 methyl/filter/TO 方向 NEGATIVE | L2–L3 | conflict with general rescue／validation wording | `InterSubMod/docs/experiments/INDEX.md` |

**Conflict count：3 個實質 conflict families**；反例風險由 20 降為 10。

## §7 Decision Path

- **TOTAL**：70/100。
- **Verdict**：GO。
- **Decision lock**：Y。不得把 current claim ceiling升格；只有新增 C1 cellular truth、C2 identifiable CN/CCF model 或 C3 upstream truth/calibration gate 才能 reopen。
- **Next action**：建立 implementation notes；按 README → Wiki/Pages → CCU → validator → browser QA 的順序執行。

### GO 前獨立紅隊

1. **Failure mode A**：直接編 generated HTML，未修 source-of-truth，下一次 rebuild 回歸。緩解：先找 generator；不存在時明示 HTML 為 tracked source並加 residual validator。
2. **Failure mode B**：dirty worktree 中誤覆寫其他 agent work。緩解：只碰事前 target-diff 為空的檔；每階段檢查 path-scoped diff。
3. **Failure mode C**：把 patch validated 說成 public fixed。緩解：disposition 分 `local_fixed`、`PATCH_VALIDATED_ON_PINNED_CLONE / NOT_APPLIED / NOT_DEPLOYED`、`external_action_required`、`blocked_by_data`。

紅隊不降級：三個 failure mode 均可觀察且有 fail-closed gate。

## Provenance Footer

- Commit：`83741469`
- Branch：`chore/handoff-20260813`
- Build time：2026-08-13 09:52 +08:00
- Cycle：`20260813_public_docs_p0_correction`
- Audit JSON：`InterSubMod/state/cycles/20260813_public_docs_p0_correction/audit.json`
