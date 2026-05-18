<!--
建立時間: 2026-05-18
parent: p4_to_p7_final.standalone.html H4-2 fix item
purpose: P3 audit H4-2 follow-up — evaluate 8 PostToolUse Edit|Write hooks 重疊度
status: validated, no deprecation needed
-->

# H4-2 Hook Overlap Evaluation — 8 PostToolUse Edit|Write Hooks

> **2026-05-18 P3 audit H4-2 結論**: **無重疊，無 deprecate 候選**。8 hooks 各偵測不同 file pattern，職責邊界清晰。

---

## 8 Hooks 職責矩陣

| Hook | 觸發 file pattern | 目的 | Output 性質 |
|------|------------------|------|------------|
| `cpp_edit_guard.sh` | `*.cpp`, `*.hpp`, `*.h` | 寫「待編譯」標記到 `/tmp/ism_cpp_pending_compile.txt` | Advisory (echo) |
| `trigger_routing.sh` | 多 path（依 file type） | Context-aware routing guidance | Advisory (echo) |
| `evidence_ledger_sync.sh` | `evidence_ledger.jsonl` | 提醒同步 MEMORY.md / INDEX.md / CURRENT_FOCUS | Advisory (echo) |
| `terminology_guard.sh` | `autoresearch/*.json(l)`, `hypothesis_queue.json` | Canonical 字眼檢核（verdict/track/phase） | Warn non-blocking |
| `md_link_check.sh` | `*.md` | Broken relative refs 偵測（呼叫 check_md_links.py） | Warn list broken |
| `kb_sot_guard.sh` | 任何 file 含 F1 SoT 數字 (+0.0112 / 0.9762 / -0.0206 / 0.9650) | 防 SoT 數字非 KB 來源被改 | Warn (block 候選) |
| `standalone_trigger.sh` | `docs/reports/.../*.md` (validated reports) | 提醒寫 standalone HTML companion | Advisory (echo) |
| `skill_change_audit.sh` | `.claude/skills/*/SKILL.md` 等 | Append audit entry 到 monthly log | Advisory + log |

---

## 重疊度分析

### Pairwise File Pattern Overlap

| Pair | Overlap? | 說明 |
|------|---------|------|
| cpp_edit_guard × trigger_routing | 部分 | trigger_routing 對 .cpp 也 fire，但只提供 routing；cpp_edit_guard 寫具體 marker |
| evidence_ledger_sync × terminology_guard | 部分 | 都對 autoresearch JSON fire，但前者提醒同步，後者檢字眼 |
| md_link_check × standalone_trigger | 部分 | 都對 .md fire，但前者檢 links，後者提醒 HTML companion |
| skill_change_audit × trigger_routing | 部分 | 都可對 skills 目錄 fire，但 audit 寫 log，routing 給 advisory |

→ **「部分重疊」皆為相同 file 觸發但作用不同**（不是功能重複）。

### Output 方向比較

| 維度 | 7 個 hooks | 1 個 hook (kb_sot_guard) |
|------|-----------|------------------------|
| Advisory (echo) | ✅ | — |
| Side effect (log/marker) | 2 個 (cpp / skill_audit) | — |
| Warn non-blocking | 2 個 (terminology / md_link) | — |
| Potentially block (Hard Gate) | — | ✅ SoT 改動候選 |

---

## Deprecate 候選評估

**無 hook 該 deprecate**：
- 每 hook 都有不可由其他替代的具體功能
- 即使部分 file pattern 重疊，output 性質不同（log vs advisory vs warn）
- 性能影響低（each hook ~10-50ms，不影響 turn latency 感知）

**可能整合候選（未來 P4+）**：
- ❌ `cpp_edit_guard` + `trigger_routing`: 整合會把 advisory 邏輯混入 marker 邏輯，降可讀性 — 不建議
- ❌ `md_link_check` + `standalone_trigger`: 不同 output 性質，整合複雜度 > 收益
- ❌ `evidence_ledger_sync` + `terminology_guard`: 同 file 但邏輯獨立 — 整合需大改

**結論**: 維持 8 hooks 現狀。

---

## 與 P3 H4-1 對比

| ID | 原計劃 | 實際發現 | 結論 |
|----|--------|---------|------|
| **H4-1** PreToolUse Bash 5 hooks | 整合 router_pre_bash.sh | 待評估（下一個 fix）| TBD |
| **H4-2** PostToolUse Edit\|Write 8 hooks | deprecate 2-3 個 | **無重疊** | **No-op** |

→ P3 audit estimate over-aggressive，H4-2 不必動。Hook architecture 設計健康。

---

## 後續監控

`InterSubMod/scripts/hooks/skill_change_audit.sh` 已對 skill 改動建立 monthly audit log（`InterSubMod/docs/postmortems/skill_audit_YYYYMM.log`）。
類似的 hook health monitoring 可加：
- `hook_latency_log.sh` (P4 candidate): 每 hook 跑完寫 latency 到 log
- 月度 review log 看是否有 hook latency drift > 100ms

---

## 元層說明

本 evaluation 源自 2026-05-18 P3 audit H4-2 (`InterSubMod/docs/reports/validated/2026/05/20260518_agent_harness_audit_p1_skills_01/p4_to_p7_final.standalone.html`)。原始 fix proposal 為 "deprecate 候選 2-3 個 hook"，實際 audit 後發現無重疊，**結論不採納 deprecate 動作**。
