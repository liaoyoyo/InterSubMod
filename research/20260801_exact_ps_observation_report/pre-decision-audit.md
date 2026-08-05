<!--
建立時間: 2026-08-01 17:01 +0800
目標: exact-PS 全資料完成後自動產生 standalone HTML observation report 的實作前審查
處理範圍: frozen all7/topology/methyl authorities、Python builder、HTML/data/receipt、pipeline final observation hook
cycle_id: cycle_20260801-exact-ps-observation-report
topic: exact_ps_observation_report
status: verdict_GO
audit_version: 0.1
關聯檔案:
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv
  - InterSubMod/.claude/skills/pre-decision-audit/SKILL.md
  - InterSubMod/.claude/skills/implementation-notes/SKILL.md
-->

# Pre-Decision Audit: Exact-PS Observation Report

> **Verdict: GO（80/100）**。這是 authority-driven 的衍生觀察層，不新增生物真值；
> 必須 fail-closed、保留分母與 abstain，且 HTML 不得成為第二份數據 SoT。

## §0 Cynefin Domain Gate

- **Domain**: Complicated。
- **Test**: repo 已有多個 Python→standalone HTML 與 Playwright QA 先例，可預測重複。
- **Rationale**: 資訊設計需要判斷，但 input/output、守恆、hash 與 claim ceiling 可明確驗證。

## §1 Observation Completeness

| Observation | 狀態 | Tier | 來源 |
|---|---|---|---|
| all7 funnel、solver、AF 數字已有 frozen authorities | ✓ | L1 | `authority_manifest.json` |
| exact topology census 已有 71,955 ranked units | ✓ | L1 | `topology_summary` artifact |
| methyl auxiliary 有正式 association-only contract | ✓ | L1 | `methyl_manifest` artifact |
| CN/LOH-aware per-unit result | □ | 尚未完成 | manifest `P0-2/P0-3` |
| 現行 HTML 可作 CN-aware biological tree claim | ✗ | L1 反例 | handoff decision spec |

## §2 Credibility Score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | 純 derived reporting；來源與 denominator contract 已固定 |
| 觀察支撐 | 20 | 7 technical datasets、chr1–22 authorities 已存在 |
| 機制清晰度 | 20 | input manifest → derive data → render → receipt → browser QA |
| 反例風險 | 10 | overclaim、stale authority、ranked-only selection bias 仍需顯式防護 |
| 所需資源 | 10 | 預估 1–6 小時含 browser QA |
| **TOTAL** | **80/100** | **GO** |

**Falsifier observable**：只要任一 authority hash、denominator conservation、
per-sample/cohort conservation、no-JS content、mobile overflow 或 print layout 失敗，
就不能發布 validated HTML。

**Reality-test 三反例**：

1. 若 HTML 只顯示 88.26% 而未顯示 ranked denominator／10,717 abstains，設計失敗。
2. 若換入 hash 不同或 schema 不符的來源仍能正常發布，fail-closed 失敗。
3. 若關閉 JavaScript 後核心數據消失，standalone observation 要求失敗。

## §3 Assumption Map

| Assumption | Importance | Known | Quadrant／處置 |
|---|---|---|---|
| manifest 指向 frozen authorities | HIGH | KNOWN | 啟動時重算 SHA-256 |
| cohort/sample totals 可由現有 JSON 完整取得 | HIGH | KNOWN | schema＋守恆測試 |
| CN/LOH 已有可展示數值 | HIGH | KNOWN=false | 顯示 `NOT_INTEGRATED`，不得推估 |
| 瀏覽器可執行 JavaScript | LOW | UNKNOWN | 核心內容 server-rendered；JS 只增強互動 |
| HTML 可取代 JSON authority | HIGH | KNOWN=false | 明文禁止；HTML 是 derived view |

## §4 Quick Pilot

不需另行 PROBE；直接以 frozen all7 bundle 作 deterministic pilot。Checkpoint：

- source hash 100% match；
- denominator／sample conservation 100% pass；
- external request=0；
- 1440、390 px、no-JS、print 均通過；
- 頁面不得出現 confirmed clone／biological tree claim。

## §5 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| 正式 report-data JSON schema | HIGH | 1h | P0 |
| builder＋standalone CSS/SVG charts | HIGH | 2–4h | P0 |
| post-run integration hook | HIGH | 1h | P0 |
| Playwright/no-JS/print receipt | HIGH | 1–2h | P0 |
| CN/LOH-aware source data | HIGH scientific | 另屬 P0-2/P0-3 | 不在本 report builder 假造 |

## §6 Evidence Conflict Scan

| Prior conclusion | Relation | 處置 |
|---|---|---|
| cellular clone/lineage 尚未驗證 | dependent | 所有卡片標 model-conditional／not integrated |
| methylation 只可 association-only | dependent | 獨立 sidecar 區，不參與 tree selection |
| ranked subset 不可外推全 groups | dependent | funnel 同時展示 98,955、85,941、75,224、71,955 |
| legacy repo 有 report-builder/finalizer split-brain | conflict risk | 新 builder 單一入口、schema/hash fail-closed、無遞迴 finalizer |

沒有與「建立 derived HTML observation」直接衝突的 concluded NEGATIVE；反例集中於
claim boundary 與 release identity。

### Red-team gate

1. **第二 SoT**：HTML 內複製數字後可能與 JSON 漂移。  
   防護：數字只由 report-data 生成，receipt 記 source SHA，HTML 標 derived。
2. **漂亮圖掩蓋 missingness**：ranked-only 88.26% 容易被誤讀。  
   防護：首頁先畫完整 funnel，abstain 與 CN unknown 使用高對比警示。
3. **錯誤整合 pipeline**：builder 失敗卻讓上游結果看似完成。  
   防護：觀察層有獨立 status；`--require-report` 時失敗回傳非零，不能覆寫 canonical data。

Red-team 通過，維持 GO。

## §7 Decision

- **Verdict**: GO
- **Decision lock**: Y
- **Next action**: 實作 authority-driven builder、data contract、receipt 與 browser QA。
- **Claim ceiling**: derived technical observation only。

## Provenance

- Commit: `387a101e6a3292e0d7f230ba8d20271c7434972a`
- Skill: `/pre-decision-audit` v0.1
- Audit JSON: `InterSubMod/state/cycles/cycle_20260801-exact-ps-observation-report/audit.json`

