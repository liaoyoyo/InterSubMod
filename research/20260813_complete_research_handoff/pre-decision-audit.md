<!--
建立時間: 2026-08-13 10:32
目標: 完整研究交接與 GitHub 公開實作前的 fail-closed pre-decision audit
處理範圍: Task Type B Comprehensive validation + D External handoff；InterSubMod 與 LongLineage
cycle_id: cycle_20260813-1032-complete-research-handoff
topic: complete_research_handoff
status: verdict_PROBE
audit_version: 0.1
build_branch: agent/research-handoff-audit-evidence-20260813
build_commit: ddd8909a838318d8a77969313e9561c8ff9d01c2
worktree: /big7_disk/liaoyoyo2001/worktrees/ism-handoff-20260813
data_sources: InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json,InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv,InterSubMod/docs/reports/validated/2026/08/20260812_InterSubMod_GitHub公開說明與教學完整驗證_01.md
驗證方式: baseline SHA、claim inventory row/count、authority artifacts 與 machine paths 均以命令 receipt 重驗；未重驗項一律不得標 current
證據等級: L2 ⭐⭐⭐⭐（工程方案可執行；公開發布與雙機驗收仍待 gate）
關聯檔案:
  - InterSubMod/.claude/skills/pre-decision-audit/SKILL.md
  - InterSubMod/.claude/skills/scientific-rigor/SKILL.md
  - InterSubMod/research/20260813_complete_research_handoff/implementation-notes.md
-->

# Pre-Decision Audit：完整研究交接與 GitHub 公開實作

> **Verdict: PROBE**：可進行 clean handoff、registry、claim 修正與 private preview PR；InterSubMod tag／GitHub Release、真正雙機完成宣稱與 LongLineage public visibility 均 fail closed。

## Frontmatter

- **Triggered by**: new-spec + comprehensive-validation + external-handoff
- **Task types**: B + D
- **Serves**: G1、G3、G4、G5
- **Cycle ref**: `InterSubMod/state/cycles/cycle_20260813-1032-complete-research-handoff/`
- **InterSubMod baseline**: `ddd8909a838318d8a77969313e9561c8ff9d01c2`
- **LongLineage candidate**: `b9aaa12a11fa00606bd174dabd0f172a5d112359`

## §0 Cynefin Domain Gate

- **Domain**: Complicated
- **Test**: Git freeze、hash replay、schema validation、license audit、build/test、site preflight 都有可重複的工程方法；但資料來源分散、歷史版本多，需專家分析與交叉核對。
- **Boundary**: 若遇到 provenance 無法解釋或 license 來源不明，不以「合理猜測」補洞；該項改標 `MISSING`／`UNVERIFIED`／`BLOCKED`。

## §1 Observation Completeness

| Observation | 狀態 | Tier | 來源與 fresh check |
|---|---:|---:|---|
| InterSubMod release baseline 可建 clean worktree | ✓ | L1 | `git status --porcelain` = 0；HEAD=`ddd8909a…` |
| LongLineage preview candidate 可建 clean worktree | ✓ | L1 | `git status --porcelain` = 0；HEAD=`b9aaa12…` |
| frozen authority 已定義 full scope | ✓ | L1 | authority scope = 7 technical datasets／6 biological IDs／chr1–22 |
| denominator registry 可機械讀取 | ✓ | L1 | TSV data rows = 19 |
| 2026-08-12 public claim audit 有完整母體 | ✓ | L1 | `claim_inventory.tsv` rows = 158；unique IDs = 158 |
| public claim 問題量已可重算 | ✓ | L1 | verdict counts 69/31/28/26/4；P0=34；需處理=58 |
| 兩台機器 fresh-clone 全 gate | □ | — | bip7 尚待跑；bip8 尚待可達性與實機 receipt |
| LongLineage 21 source mappings 與 license gate | □ | — | 待獨立 audit；未通過不得切 public |
| authority 現在仍為 19/19 hash match | ✓ | L1 | 2026-08-13 fresh replay：13 artifacts＋1 binary＋5 source snapshots，全數 MATCH |

## §2 Credibility Score

| 維度 | 暫定評分 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | provenance registry、hash、schema、CI 與 license audit 是成熟做法 |
| 觀察支撐 | 10 | frozen authority 與 public audit 完整，但 18/19/35、雙機與 license disposition 未閉合 |
| 機制清晰度 | 20 | producer→input→artifact→claim→consumer 可建立明確 DAG |
| 反例風險 | 0 | dirty branches、source adjudication、跨機 mount 與 10 筆未驗 source mapping 是高風險反例 |
| 所需資源 | 0 | 全面盤點、雙 repo、雙機與多個 release gate 明顯超過 6 小時 |
| **TOTAL（red-team 後）** | **50 / 100** | `PROBE`；實作可推進，發布不可推進 |

**Falsifier observable**：若任何 final-for-scope artifact 無法以 producer commit、inputs、command receipt、schema、scope、hash 與 machine location 形成 closed evidence chain，或 frozen authority replay 出現 missing/mismatch，則不得宣稱 release-ready；相關 gate 必 fail closed。

**Reality-test 反例**：

1. `authority_manifest.json` 仍存在歷史 source SHA adjudication，可能讓「19/19 artifact match」與「source 可重建」被誤當同一件事。
2. 08-12 的 `claim_inventory.tsv` 雖存在且有 checksum，但在來源 branch 被 ignore，fresh clone 不一定取得。
3. bip8 若無獨立 host 執行 receipt，只能標 `BIP8_DATA_PREFLIGHT_BLOCKED`，不能用 bip7 掛載的 NFS 代替雙機驗收。

## §3 Assumption Map

| Assumption | Importance | Known | Quadrant |
|---|---:|---:|---:|
| 2026-08-01 frozen authority 仍未漂移 | HIGH | UNKNOWN | 2 — 先 hash replay |
| `ddd8909a` 是本次 InterSubMod 唯一 baseline | HIGH | KNOWN | 1 — user locked |
| `b9aaa12` 是 LongLineage 唯一 preview candidate | HIGH | KNOWN | 1 — user locked |
| public claims 可由單一 registry 驅動 | HIGH | PARTIAL | 2 — validator + entrypoint crosswalk |
| bip8 可直接登入並做 fresh clone | HIGH | UNKNOWN | 2 — doctor receipt；否則 blocked |
| 所有來源可再散布 | HIGH | UNKNOWN | 2 — per-asset license audit |
| TiB science 必須本輪重跑 | LOW | KNOWN false | 3 — 本輪明示不重跑 |

## §4 Quick Pilot

本任務不採 subset science pilot；最小 safe-to-fail probe 是四個唯讀 gate：

1. fresh replay authority hashes；任一 mismatch → `AUTHORITY_REPLAY_BLOCKED`。
2. 解析 claim inventory 158 unique IDs；缺列或 checksum mismatch → `CLAIM_REGISTRY_BLOCKED`。
3. 兩個 clean worktree HEAD/dirty check；不符 baseline → `RELEASE_FREEZE_BLOCKED`。
4. LongLineage license/source-origin scan；有未解析來源 → 保持 private。

## §5 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---:|
| authority fresh 19-path replay | HIGH | <1h | P0 |
| 18 manifest rows／19 documented runs／35 directories reconciliation | HIGH | 1–3h | P0 |
| LongLineage license/source-origin/SBOM PASS | HIGH | 2–6h | P0 |
| bip7/bip8 independent receipts | HIGH | 1–6h＋外部可達性 | P0 |
| all public entrypoints from one claim registry | HIGH | 3–6h | P0 |
| full-history secret scan | HIGH | 1–3h | P0 |

## §6 Evidence Conflict Scan

| Prior conclusion | Relation | Source |
|---|---|---|
| Thread B cross-sample whitelist 已撤回 | dependent：交接不能復活舊 filter claim | `InterSubMod/docs/reports/validated/2026/04/20260426_Thread_B_Retraction_Whitelist_Cross_Sample_01.md` |
| confirmed cellular subclone／linear ancestry 皆為 0 | hard ceiling | `InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json` |
| CN/LOH 未整合；methylation association-only | hard ceiling | 同上 |
| public docs 58/158 claims 有問題、P0=34 | release blocker | `InterSubMod/docs/reports/validated/2026/08/20260812_InterSubMod_GitHub公開說明與教學完整驗證_01.md` |
| `73afaeac` 與 08-13 drilldown/CNV 是 active/partial | exclude from release snapshot | 使用者鎖定規格 |
| LongLineage P3/P4/P5/P7/P8 blocked | research-preview ceiling | LongLineage governance，待本輪獨立核對 |

**Conflict count**: 6；都已轉成 release gate 或 claim ceiling，不允許靜默覆蓋。

## §7 Decision Threshold + Path

- **TOTAL**: 50/100
- **Verdict**: `PROBE`
- **Decision lock**: Y
- **Independent red-team**: PASS（已提供有效反例與 fail-closed gate，不代表 release PASS）
- **Scoped GO**: clean worktree、freeze identity、registries、handoff docs、claim corrections、portable interfaces、private LongLineage license/preview PR。
- **NO-GO until receipts**: InterSubMod tag／Release、LongLineage public visibility、LongLineage production release、`BIP8_*_PASS` 宣稱。

### Red-team failure modes

1. 把 `ddd8909a` clean build 誤稱為可重生 2026-08-01 frozen science。
2. 把 bip7 上可讀的 bip8 NFS mount 誤稱為 hostname=`bip8` 的獨立驗收。
3. 把 LongLineage preview 可建置誤稱為 source-origin/license 已可公開或 production ready。

### Mandatory fail-closed order

1. baseline identity → 2. authority 19/19 → 3. 18/19/35 reconciliation → 4. registry validation → 5. 158-claim closure → 6. hygiene → 7. portable execution → 8. dual-host receipts → 9. GitHub consistency → 10. LongLineage public gate。

## Provenance Footer

- **Commit**: `ddd8909a838318d8a77969313e9561c8ff9d01c2`
- **Built**: 2026-08-13 10:32:12 +08:00
- **Audit JSON**: `InterSubMod/state/cycles/cycle_20260813-1032-complete-research-handoff/audit.json`
- **Source workspace observed at**: `95d420f6f397487822ce97dbab8151eb9e56af53`
- **Note**: source workspace 僅供唯讀取證；release worktree 不繼承其 dirty state。
