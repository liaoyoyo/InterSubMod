<!--
建立時間: 2026-07-14 22:16 +08:00
目標: Verification schema v2 與 RegionStratifier 修復進入實作前 pre-decision audit
處理範圍: verification-region-stratifier-schema-repair
cycle_id: cycle_20260714-2216-verification-region-stratifier-schema-repair
topic: verification-region-stratifier-schema-repair
status: verdict_GO
audit_version: 0.1
關聯檔案:
  - /big7_disk/liaoyoyo2001/external_validation/subclone_reconstruction_method_comparison_20260714/20260714_InterSubMod分類與RegionStratifier修復任務規格_01.md
  - InterSubMod/.claude/skills/pre-decision-audit/SKILL.md
  - InterSubMod/.claude/skills/scientific-rigor/SKILL.md
  - InterSubMod/.claude/skills/implementation-notes/SKILL.md
-->

# Pre-Decision Audit: Verification schema v2 與 RegionStratifier 修復

> **Task type**: E Hotfix / Bugfix + B Comprehensive validation
> **服務目標**: G4 reproducibility、G5 externally verifiable engineering
> **Verdict**: `GO_WITH_FAIL_CLOSED`（獨立紅隊降級後，四個 bounded probes 已全數通過）
> **Scope**: 全部指定 source/header/consumer、14-file frozen corpus；不是 subset。

## Frontmatter

- **Topic**: `verification-region-stratifier-schema-repair`
- **Triggered by**: `new-spec`
- **AI session**: Codex 2026-07-14
- **Last updated**: 2026-07-14 22:28 +08:00
- **Cycle ref**: `InterSubMod/state/cycles/cycle_20260714-2216-verification-region-stratifier-schema-repair/`

---

## §0 🟤 Cynefin Domain Gate

- [x] **Domain**: Complicated
- [x] **Test**: deterministic schema mapping、fixed precedence、atomic-output contract 與 golden counts 都有可預測結果。

**Domain decision**: `Complicated`。

**Rationale**: 這是已定位的 semantic/schema defect；修復規格提供 ordered truth table、狀態機、欄位值域與 frozen expected counts。難點是跨 consumer 相容性與 identity fail-closed，不是未知研究機制。

---

## §1 🟢 Observation Completeness Checklist

| Observation | 狀態 | Evidence tier | 來源 |
|---|---:|---:|---|
| Stage④ 目前把 legacy `Strong` 與 `Subclone` 都輸出成 current `Strong` | ✓ | L1 ⭐⭐⭐⭐⭐ | `InterSubMod/src/core/RegionProcessor.cpp:1319` |
| Phase D 目前是固定 threshold stratification，且 precedence 為 LOH→sample ASM→HP ASM→baseline | ✓ | L1 ⭐⭐⭐⭐⭐ | `InterSubMod/src/core/SubcloneAnalyzer.cpp:100` |
| `n_subclones=max_id+1` 對 sparse IDs 會錯計 occupied groups | ✓ | L1 ⭐⭐⭐⭐⭐ | `InterSubMod/src/core/SubcloneAnalyzer.cpp:35` |
| caller 與 analyzer 各自硬編 eligibility，write-back 依 profile position | ✓ | L1 ⭐⭐⭐⭐⭐ | `InterSubMod/src/core/RegionProcessor.cpp:712`; `InterSubMod/src/core/SubcloneAnalyzer.cpp:67` |
| frozen manifest 的 14 檔 SHA-256 與 row counts 全通過 | ✓ | L1 ⭐⭐⭐⭐⭐ | `FROZEN_CORPUS_MANIFEST.tsv` + 本 cycle shell audit (`FROZEN_MANIFEST_PASS files=14`) |
| baseline build、`run_tests`、`ctest` 全綠 | ✓ | L1 ⭐⭐⭐⭐⭐ | build exit 0；234/234 GoogleTest；234/234 CTest |
| target core/consumer 檔沒有與本修復重疊的 content diff | ✓ | L1 ⭐⭐⭐⭐⭐ | 2026-07-14 targeted `git diff` inventory；既有 tests 僅 mode change |
| 所有 executable consumer 的精確用途與 fallback 已逐檔判定 | □ | — | 三路 read-only audit 進行中 |

---

## §2 🔵 Credibility Score

| 維度 | 評分 | 理由 |
|---|---:|---|
| **理論基礎** | 20 | 以 schema normalization、provenance preservation、fail-closed 與 atomic publish 為明確工程原則 |
| **觀察支撐** | 20 | direct code facts + 7-sample/14-file frozen corpus + green baseline |
| **機制清晰度** | 20 | producer→CSV→consumer 與 analyzer→write-back→owned artifact 的 defect chain 明確 |
| **反例風險** | 10 | consumer inventory、ABI/additive compatibility 與 stale-artifact failure path 尚待 targeted tests |
| **所需資源** | 0 | 完整 consumer surface 為 62 支，實際修復與驗證可能超過 6 小時 |
| **TOTAL** | **70 / 100** | 四項紅隊 probes 通過後，允許 bounded engineering GO |

**Falsifier observable (WRAP)**：

> 若這個修復可在不改科學判斷的前提下完成是錯的，將看到至少一項：逐 stable key 的 `Significant` 改變；需要改動既有 predicate/threshold 才能過測；frozen class mapping 不再是 59,910 / 6,931 / 0；或 fail-closed identity test 仍產生部分 write-back。

**Reality-test 三個反例觀察**：

1. `Strong_Bidirectional + ClusterFirstOnly` 不等於原 v1 `Strong=66,841`。
2. sparse IDs `{3}` 或 `{0,3}` 的 occupied count 仍被 `max_id+1` 汙染。
3. paired insufficient、tumor-only、FAILED run 仍留下前次 VALID assignment/summary。

---

## §3 🟡 Assumption Map 2×2

| Assumption | Importance | Known | Quadrant | Resolution |
|---|---|---|---|---|
| 14 個 frozen input 都含有效 legacy provenance | HIGH | KNOWN | (1) | manifest/hash/column audit；migration 仍逐檔 fail closed |
| binary retention 可由同一 predicate、只拆 class 名維持不變 | HIGH | KNOWN | (1) | ordered truth table + per-row Significant golden compare |
| 所有 assignments 可攜帶原始 result index 與 region ID | HIGH | KNOWN | (1) | additive typed assignment struct |
| 四支已知 consumer 以欄名而非 undocumented column offset 讀取 | HIGH | UNKNOWN | (2) ⚠ | consumer agent 逐檔 code audit + fixtures，未確認前不改 |
| 其他 executable consumers 不會被裸字串 fallback 靜默污染 | HIGH | UNKNOWN | (2) ⚠ | bounded `rg` inventory，依 current/legacy/evidence/stratum 分類 |
| 可保留 `SubcloneAnalyzer` 類名作薄 compatibility surface | LOW | UNKNOWN | (4) | correctness 優先；不做無必要 rename |

**右上 (2) quadrant 假設**：consumer 欄位選擇與其他 executable consumer surface。由 read-only agents 與 targeted regression 在實作前/中消除。

---

## §4 🟠 Quick Pilot Guide

本任務不是探索假說，且 frozen full-scope audit 成本低；不另做 subset pilot。最小 safe-to-fail check 已由以下取代：

1. baseline build/test：234/234 pass。
2. frozen 14-file SHA/row audit：14/14 pass。
3. 先寫 isolated unit fixtures，再跑 migration full corpus；任何 falsifier 命中即停止。

---

## §5 🔴 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---:|
| sparse-ID deprecated ABI 是否仍有 ID-index consumer | HIGH | <0.5h | P0 |
| multi-artifact crash 是否能以 status-first/final-last fail closed | HIGH | <0.5h | P0 |
| core C++ API/owned artifacts 精確變更面 | HIGH | <1h | P0 |
| 四支 consumer + 其他 executable consumer inventory | HIGH | <1h | P1 |
| migration preservation/golden harness | HIGH | 1–2h | P2 |

---

## §6 🟣 Evidence Conflict Scan

- `InterSubMod/MEMORY.md` 不存在，已記錄為 repository fact；未假造 concluded entries。
- `InterSubMod/docs/reports/validated/` 沒有檔名含 `NEGATIVE` 的檔案。

| Prior conclusion | Relation | Source |
|---|---|---|
| 甲基化 variant filter 已 DEAD | dependent、非 conflict：本修復明確禁止改 retention/predicate | `InterSubMod/docs/CURRENT_FOCUS.md:150`; `:204` |
| bulk regional state 不得稱 confirmed cellular clone | support：支持把固定規則輸出正名為 RegionStratum | `InterSubMod/docs/CURRENT_FOCUS.md:27`; `InterSubMod/docs/experiments/INDEX.md:32` |
| historical Phase D 曾稱 4-group subclone | affected artifact，需重生敘述；不是阻止 correctness fix 的反證 | `InterSubMod/docs/experiments/INDEX.md:260` |

**Conflict count**: `0` 真衝突；有 2 個 dependent constraints 已鎖入驗收。

---

## §7 ⚫ Decision Threshold + Path

- **TOTAL score**: `70 / 100`
- **Current verdict**: `GO_WITH_FAIL_CLOSED`（22:28 bounded probes 全 PASS）。
- **Scope boundary**: 只驗 schema/serialization/consumer provenance；不把 migration PASS 說成 upstream statistics 或 biological validity 的獨立重證。

### 獨立紅隊最強反方

1. sparse enum ID 與 deprecated `n_subclones=occupied count` 可能破壞把 count 當 ID-index bound 的舊 consumer。
2. per-file atomic rename 不能單獨避免 crash 中斷後 status 與 artifacts 混代。
3. offline golden 若只相信舊 class，不能獨立否證 upstream legacy classification 本身錯誤。
4. executable consumer surface 遠大於規格列出的四支；silent fallback 會把新 structure class 落入 Noise。
5. 若 consumer 把 `ClusterFirstOnly` 或 RegionStratum 當 cellular clone，會違反既有 concluded scope。

### PROBE（≤30 min）

1. `rg` 精確盤點 `n_subclones` / `Subclone_ID` executable readers → checkpoint：沒有 ID-index bound consumer，或可提供 explicit compatibility adapter。
2. 畫出 publish ordering → checkpoint：run 開始先原子寫 `FAILED/RUN_IN_PROGRESS`，data artifacts 與 main CSV 成功後最後原子發布 final status；任一失敗不能呈現 VALID。
3. 對照 handed spec、QA 與 sibling document → checkpoint：`ClusterFirstOnly` 有單一明確 SoT，衝突來源標為 superseded/history。
4. 完成 consumer inventory 分群 → checkpoint：所有 canonical readers 都有 current/legacy/evidence/stratum 明確選擇；無法判斷者 fail closed 或 explicit CLI。

### PROBE 結果

| Probe | Verdict | 可觀察證據 |
|---|---:|---|
| sparse-ID compatibility | PASS | `n_subclones` 無 repo downstream executable consumer；`Subclone_ID` 兩個命中只做歷史空值/文字報告，沒有 ID indexing |
| crash-safe publish | PASS（條件式） | 先 atomic `FAILED/RUN_IN_PROGRESS`，所有 data + main CSV 成功後最後發布 terminal status；single-writer contract |
| enum name authority | PASS | handed spec §207–210、§581 鎖 `ClusterFirstOnly`；`ClusterOnlyAssociation` 僅較早 display wording |
| bounded falsifier | PASS | C++ legacy 2×2 + Stage④ truth table + frozen stable-key/byte preservation 足以否證本次 mapping regression；不涵蓋上游 biological truth |
| consumer surface | PASS（需實作） | `consumer_inventory.md` 完整列 62 VerificationClass、23 LOH、2 RegionStratum-related consumers，分為 C2/L4/E/R1/LOH-L/P/H1 與 B0/B1/B2/H |

### GO hard stops

- target source/test 出現無法安全合併的重疊 user content diff；
- status-first / status-last commit-marker protocol 無法落地；
- 任一 `Significant`、既有 raw cell、TP/FP/FN/F1 或 statistical predicate 改變；
- 仍有未分類 live consumer 或 unknown class 靜默歸 Noise；
- scope 滑向 cellular clone/CCF/tree inference。

### Step → Verify

1. 鎖定 v2 producer schema → 驗證：11 truth-table fixtures + precedence collisions 全通過。
2. 正名並修 RegionStratifier → 驗證：0/49/50、sparse ID、identity mismatch、四種 status 與 atomic artifacts 全通過。
3. 遷移 consumers → 驗證：各 consumer 固定 fixture 明確 current/legacy/evidence/stratum source，unknown 不落 Noise。
4. frozen migration → 驗證：14/14 hash、328,697 rows、59,910/6,931/0、HCC1395 1,516，既有 cells/raw order 不變。
5. 全 repo validation → 驗證：build exit 0、`run_tests`/`ctest` 無新增失敗、`git diff --check` exit 0。

### Decision lock

- [x] **Lock user-specified contracts after red-team pass**: v2 strings/precedence、thresholds、status values、frozen counts、no in-place migration。
- [x] 不改 statistical thresholds、不把 region strata 稱為 cellular clones、不重開 methylation filter。

---

## Provenance Footer

- **Commit hash**: `6067568637088838a9f518955e41d222f057e4f1`
- **Build time**: 2026-07-14 22:16:50 +08:00
- **Skill version**: `/pre-decision-audit` v0.1
- **Audit JSON**: `InterSubMod/state/cycles/cycle_20260714-2216-verification-region-stratifier-schema-repair/audit.json`
- **Linked topic**: `InterSubMod/research/20260714_verification_region_stratifier_schema_repair/`
