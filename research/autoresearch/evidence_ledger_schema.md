<!--
建立時間: 2026-05-18
parent: P2 audit M13 fix (P3 紀錄 7 欄補完)
status: validated
purpose: evidence_ledger.jsonl schema 文件化，含 v2 升級 (next_action + identified_issues)
-->

# evidence_ledger.jsonl Schema (v2)

> **2026-05-18 v1 → v2 升級**：依 P2 audit M13 補 2 個 action-oriented 欄位（next_action + identified_issues），對齊業界紀錄 7 欄標準（日期/目標/做法/結果/問題/原因/下一步）。
> **向後相容**：v1 entries 不必 migrate；新 entries 應含新 2 欄。

---

## Schema Fields (v2)

### Core identification (v1, required)

| 欄 | 型別 | 範例 | 說明 |
|---|------|------|------|
| `cycle_id` | string | `"2026-05-18_v6_marker"` | 短 ID，含日期 + 假說簡稱 |
| `hypothesis_id` | string | `"H_THREAD_D_MAIN"` | 對齊 hypothesis_queue.json |
| `hypothesis` | string | `"V6 binary marker coverage > V3F"` | 一句話 |
| `timestamp` | ISO 8601 | `"2026-05-18T14:00:00+08:00"` | 寫入時間 |

### Pipeline & scope (v1, required)

| 欄 | 型別 | 範例 |
|---|------|------|
| `pipeline_track` | enum | `"paired_full" \| "paired_pileup" \| "TO" \| "cross-platform"` |
| `datasets_tested` | array | `["HCC1395_5kHz", "COLO829", "H2009"]` |
| `scale` | enum | `"single_sample" \| "cross_sample" \| "full_genome"` |

### Quantitative results (v1, required)

| 欄 | 型別 | 範例 |
|---|------|------|
| `baseline_f1` | float | `0.7167` |
| `result_f1` | float | `0.7185` |
| `delta_f1` | float | `+0.0018` |
| `delta_tp` | int | `+12` |
| `delta_fp` | int | `-8` |
| `delta_fn` | int | `-4` |

### Verdict & metadata (v1, required)

| 欄 | 型別 | 範例 |
|---|------|------|
| `verdict` | enum | `"POSITIVE" \| "NEGATIVE" \| "NEUTRAL" \| "PARTIAL"` |
| `tier_used` | int | `1-5` （⭐1 最低 / ⭐5 最高 — 對齊 `/scientific-rigor §2`）|
| `human_decision` | enum | `"keep" \| "revert" \| "extend"` |
| `artifacts_path` | string | `"InterSubMod/research/autoresearch/cycles/<id>/"` |
| `methodology_doc` | string | `"InterSubMod/docs/methodology/20260518_xxx.md"` |

### Qualitative observations (v1, required)

| 欄 | 型別 | 範例 |
|---|------|------|
| `key_observations` | string | `"V6 在 5/7 樣本 marker coverage +9%；HCC1395 chr8 例外"` |
| `feature_correlations_noted` | array | `["HP×LOH 強相關", "AF×CN 弱相關"]` |
| `human_notes` | string | `"Pass 2 BAM 待對比"` |

### Forward-looking (v1, optional)

| 欄 | 型別 | 範例 |
|---|------|------|
| `research_potential` | enum | `"high" \| "medium" \| "low"` |
| `orthogonality` | float | `0.0-1.0`（與既有 feature 正交度）|
| `mechanism_clarity` | enum | `"clear" \| "partial" \| "unclear"` |
| `combination_ready` | bool | `true / false` |

### **v2 新增（2026-05-18 M13 fix）— Action-oriented 欄位**

| 欄 | 型別 | 範例 | 為何加 |
|---|------|------|--------|
| **`next_action`** | string | `"擴展到 H1437 + H2009 cross-sample 驗證"` | 對齊「下一步」業界紀錄標準；避免 ledger 變純歷史紀錄無 actionability |
| **`identified_issues`** | array | `["COLO829 truth set 0600 權限", "Pass 2 BAM 缺失"]` | 結構化 issue list；對應 `/scientific-rigor §9.2` Postmortem action items 來源 |

**為何只加這 2 欄不加更多**：
- v1 已有 `human_notes` / `key_observations`（observation 類）+ verdict（result 類）
- 缺的是「下一步行動」+ 「結構化問題」— 對應 PDCA Act 階段 + Postmortem action items
- 對齊 `templates/postmortem.md` 8 段中的 timeline + action items 區段

---

## v2 完整範例 entry

```json
{
  "cycle_id": "2026-05-18_v6_marker_coverage",
  "hypothesis_id": "H_V6_PROD",
  "hypothesis": "V6 binary marker coverage 跨 6 樣本 ≥ V3F baseline",
  "pipeline_track": "paired_full",
  "datasets_tested": ["HCC1395", "H1437", "H2009", "HCC1954", "HCC1937", "COLO829"],
  "scale": "cross_sample",
  "baseline_f1": 0.7167,
  "result_f1": 0.7185,
  "delta_f1": 0.0018,
  "delta_tp": 12,
  "delta_fp": -8,
  "delta_fn": -4,
  "key_observations": "V6 在 5/6 樣本 marker coverage +9.4%；COLO829 deferred",
  "feature_correlations_noted": ["HP×LOH 強相關"],
  "verdict": "POSITIVE",
  "human_decision": "keep",
  "human_notes": "符合 SEQC2 baseline",
  "timestamp": "2026-05-22T16:00:00+08:00",
  "tier_used": 3,
  "artifacts_path": "InterSubMod/research/v6_bam_tpfp_hp_loh_cn/",
  "research_potential": "high",
  "orthogonality": 0.45,
  "mechanism_clarity": "partial",
  "combination_ready": false,
  "methodology_doc": "InterSubMod/docs/methodology/20260518_v6_marker.md",
  "next_action": "W4 啟動 paired Tier 2 Archive TO 7-sample rerun + COLO829 truth set 權限解封流程",
  "identified_issues": [
    "COLO829 truth set 0600 權限阻塞，需 owner 改權限或 accept limitation",
    "Pass 2 BAM 對比待補（5 commits chain 中 ploidy bug 4-30 修了但 Pass 2 仍未跑）",
    "HCC1395 chr8 marker coverage 異常拉高 mean，需用 median 報告"
  ]
}
```

---

## 升級指引

### 新 entry（cycle 結束時寫入）

✅ **包含 v2 新 2 欄**：
- `next_action`: 一句具體 actionable 描述（≥1 個動詞 + ≥1 個 owner / deadline）
- `identified_issues`: 結構化 array（每個 issue ≤ 80 chars）

### v1 既有 entry

- ❌ **不要 migrate**（保留歷史 snapshot 完整性）
- ✅ JSON schema validator 應允許 v2 新欄為 optional

### `/conclude-research` skill 寫入時

依本 schema 強制：
- 6 個 required core fields
- 4 個 quantitative results
- 3 個 qualitative observations
- **v2: next_action + identified_issues 強制**（若 empty array 也要顯式 `[]`）

---

## 相關

- **Skill**: `/conclude-research` P6 寫入 + `/run-evaluator` P5 標 tier
- **Hook**: `evidence_ledger_sync.sh` PostToolUse 同步檢查
- **Audit**: `/provenance-tier-audit` 跨 artifact 一致性
- **Postmortem 對齊**: `templates/postmortem.md` 8 段中 timeline / action items 由 v2 新 2 欄填充
- **Plan**: P2 audit M13 (`InterSubMod/docs/reports/validated/2026/05/20260518_agent_harness_audit_p1_skills_01/p4_to_p7_final.standalone.html`)
