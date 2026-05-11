---
id: ism-kb-09-conclusions-index
name: "Conclusions 索引"
description: "研究結論狀態對照表：positive findings / characterization / NEGATIVE / ongoing；連結到 docs/reports/research_landscape/。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "conclusions status against MEMORY 20+ concluded entries"
related_ids:
  - ism-kb-09-conclusions-positive-findings
  - ism-kb-09-conclusions-characterization-only
  - ism-kb-09-conclusions-concluded-negative
  - ism-kb-09-conclusions-research-landscape-index
  - ism-kb-09-conclusions-hypothesis-queue-snapshot
tags: [conclusions, index, status, research]
canonical_paths: [09_conclusions/00_index.md]
alias_paths: []
---

# Conclusions 索引

- 一句結論：ISM 研究結論 20+ 條，分為 positive / characterization / NEGATIVE / ongoing；**查結論前必讀**
- 適用對象：AI agent（避免重新調查已結案方向）、新進研究者
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  cat /bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/MEMORY.md | head -60
  ```

---

## 狀態對照表

| 狀態 | 含義 | 文件 |
|------|------|------|
| 🟢 **Positive** | 已驗證為真實信號，可作分析特徵 | [01_positive_findings.md](01_positive_findings.md) |
| 🟡 **Characterization** | 真實現象，**不可做 variant filter** | [02_characterization_only.md](02_characterization_only.md) |
| 🔴 **NEGATIVE / NO-GO** | 已證偽，別重新調查 | [03_concluded_negative.md](03_concluded_negative.md) |
| 🟣 **Ongoing** | 進行中 | [../10_research_status/](../10_research_status/) |

---

## 文件列表

| 檔案 | 主題 |
|------|------|
| [01_positive_findings.md](01_positive_findings.md) | Positive findings 索引（HPFineNGroups / LOH×AF / Self-Phasing / Zone-Aware H1 等） |
| [02_characterization_only.md](02_characterization_only.md) | 能描述但不可做 filter（HPFineNGroups filter / Zone QS 等） |
| [03_concluded_negative.md](03_concluded_negative.md) | NEGATIVE 目錄（O9-O13 / Option C / Fine-Pairwise / TO Germline FP / ...） |
| [04_research_landscape_index.md](04_research_landscape_index.md) | 11 份 research_landscape 報告的導航 |
| [05_hypothesis_queue_snapshot.md](05_hypothesis_queue_snapshot.md) | hypothesis_queue.json 快照 |

---

## 閱讀順序建議

1. **先查 status**：查某主題狀態 → 00_index（本文件）
2. **positive 想深入**：01_positive_findings → 07_derived_features 特徵細節 → docs/reports/research_landscape/
3. **避免重做**：03_concluded_negative → 看失敗原因
4. **整體推論鏈**：04_research_landscape_index

---

## 核心原則

- **本目錄只做索引**，詳細內容在 `docs/reports/research_landscape/`
- 狀態變更時同步更新本索引與對應 docs 文件
- AI agent 回答結論類問題時，**優先查此索引**，再跳轉 landscape

---

## 相關

- Research landscape 權威：[../../docs/reports/research_landscape/00_INDEX.md](../../docs/reports/research_landscape/00_INDEX.md)
- Derived features：[../07_derived_features/00_index.md](../07_derived_features/00_index.md)
- Research status：[../10_research_status/00_index.md](../10_research_status/00_index.md)
