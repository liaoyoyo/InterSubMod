<!--
建立時間: 2026-07-24 23:50
目標: 將 20260724 exact-PS×HP topology authority 升為 layered_workstation 預設資料前的證據稽核
處理範圍: Task Type B / 7 datasets / chr1-22 / HTML authority promotion
cycle_id: ad-hoc-20260724-layered-workstation-exact-ps
topic: layered-workstation-exact-ps-authority-promotion
status: verdict_GO
audit_version: 0.1
關聯檔案:
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/README.md
  - InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/pre-decision-audit.md
  - InterSubMod/.claude/skills/pre-decision-audit/SKILL.md
-->

# Pre-Decision Audit: layered workstation exact-PS authority promotion

> **Verdict: GO（70/100）**。新資料可升為預設展示，但必須把
> `technical_all_pass`、`ranked_complete` 與 `topology-complete` 分開呈現；
> 不能把 guard-triggered ABSTAIN 隱藏成已解析拓撲。

## §0 Cynefin domain gate

- **Domain**: Complicated。
- **Test**: 已有可重複的 standalone HTML builder、receipt/SHA gate 與
  Chromium regression audit；資料語意複雜，但產生與驗證流程可預測。
- **Rationale**: 這是 authority/schema migration，不是探索未知生物機制。

## §1 Observation completeness

| Observation | 狀態 | Tier | 來源 |
|---|---:|---:|---|
| 7/7 dataset、chr1–22 pipeline technical PASS | ✓ | L1 | `research/20260724_exact_ps_cpp_topology_af_all_samples/20260724_C++_exact_PS拓撲與read_AF全七樣本驗證_01.md` |
| 98,955 final groups 與各 stage 算術守恆 | ✓ | L1 | `research/20260724_exact_ps_cpp_topology_signature_census/data/20260724_exactPS_k_HP與分母重算紀錄_01.md` |
| 71,955 ranked units、680,527 best trees 全部重現 canonical score/tie | ✓ | L1 | `research/20260724_exact_ps_cpp_topology_signature_census/20260724_exactPS最佳樹拓撲Signature精確統計_01.md` |
| 新 exact topology = 88.2579%，舊 legacy 92.18% 不可沿用 | ✓ | L1 | `research/20260724_exact_ps_cpp_topology_signature_census/README.md` |
| 所有 mutation-bearing families topology-complete | ✗ | L1 | cohort receipt 明示 10,717 resource-limit ABSTAIN |
| HTML 逐區域 join 與桌機/手機互動 | □ | — | 本次實作與 Playwright audit |

## §2 Credibility score

| 維度 | 評分 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | exact PS、primary HP 與 strict read-linkage 符合區域分析單位 |
| 觀察支撐 | 20 | 7 datasets 全量 receipt、數量守恆與 signature census |
| 機制清晰度 | 20 | MLHP group → C++ topology row → exact signature row 一對一 key chain |
| 反例風險 | 10 | 有 ABSTAIN、zero denominator 與兩技術版本不可誤當獨立生物樣本 |
| 所需資源 | 0 | 全量 renderer、跨面板 regression 與視覺 audit 超過 6 小時風險 |
| **TOTAL** | **70 / 100** | **GO** |

**Falsifier observable**：若任何 sample 的 MLHP/topology row key 不一對一、
census row 不是 ranked subset、receipt SHA 不符、頁面把 ABSTAIN 當 resolved、
或 Chromium 顯示/互動 regression 失敗，就不得切換預設 authority。

**Reality-test 三反例**：

1. 新資料應不再含 `adjacent_gap<=50000` 作為 primary grouping；若頁面仍出現，
   表示舊資料滲入。
2. resolution 三類加總必須等於 71,955；各 sample 亦須守恆。
3. HCC1395 與 HCC1395_DORADO 必須標成同一 biological sample 的 technical
   concordance，不能計作兩個獨立生物驗證。

## §3 Assumption map

| Assumption | Importance | Known | Quadrant |
|---|---|---|---|
| cohort/census receipt 可重建且 SHA-bound | HIGH | KNOWN | 1 |
| `region_id`/`group_index` 在 MLHP、topology、census 一致 | HIGH | 待本次驗證 | 2 |
| census 只覆蓋 `ranked_complete` 是刻意 subset | HIGH | KNOWN | 1 |
| browser 可承受每頁完整 group index | HIGH | 待本次驗證 | 2 |
| 舊 v5 annotation 可直接搬到新 region keys | LOW | UNKNOWN | 4 |

右上象限先驗：執行 join contract test與 Chromium 大頁面效能/互動測試。

## §4 Quick pilot

本輪 verdict 已達 GO，但先用 HCC1395 做 safe pilot：

1. 驗證 MLHP/topology 11,590 keys 完全相等。
2. 驗證 census 9,130 keys 恰為 `ranked_complete` subset。
3. 產生單頁後以 1440×1000 與 390×844 Chromium screenshot 檢查。

Checkpoint：三項全 PASS 才執行 7 dataset rebuild；任一不符即 fail closed。

## §5 Gap diagnosis

| Missing | Impact | Effort | Priority |
|---|---:|---:|---:|
| 全頁面 authority/join regression | HIGH | 1–3 h | P0 |
| 手機圖表與 region drawer 稽核 | HIGH | 1–2 h | P0 |
| 舊 annotation 到 exact-PS key 的可證 join | LOW | >6 h | P2；本輪不假裝可用 |

## §6 Evidence conflict scan

- Repo root 沒有 `MEMORY.md`，故該強制查詢無檔案可讀。
- `docs/reports/validated/` 的 NEGATIVE 結論集中在 methylation/LOH/FP
  路線，與本次 renderer authority migration 無直接衝突。
- 直接相關的衝突反而來自新 receipt：`all_pass=true` 不等於
  `all_mutation_bearing_families_complete=true`。本次 UI 必須顯示該限制。

## §7 Decision path

- **TOTAL**: 70 / 100
- **Verdict**: **GO**
- **Decision lock**: Y；使用者已明確要求 `layered_workstation/` 用新資料。
- **Next action**:
  1. 以 exact-PS sources 重建 index 與 7 sample pages。
  2. 舊 50 kb 內容只保留為可摺疊 legacy baseline，不進預設統計。
  3. 完成資料 contract、跨面板與 Chromium 桌機/手機驗證。

## Provenance

- Source commit at audit: `0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- Task type: B comprehensive validation
- Skill: `/pre-decision-audit` v0.1
