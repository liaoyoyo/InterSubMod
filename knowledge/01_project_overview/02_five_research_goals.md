---
id: ism-kb-01-project-overview-five-research-goals
name: "Five Research Goals"
description: "ISM 五大研究目標定義、主軸優先順序、成功標準；對應 MEMORY project_research_vision_five_goals。"
status: active
last_verified: 2026-04-22
content_nature: frozen-decision
doc_type: explanation
verified_scope: "five goals against MEMORY project_research_vision_five_goals"
related_ids:
  - ism-kb-01-project-overview-index
  - ism-kb-01-project-overview-breakthrough-strategy
  - ism-kb-01-project-overview-project-summary
tags: [overview, goals, research, vision]
canonical_paths: [01_project_overview/02_five_research_goals.md]
alias_paths: []
---

# Five Research Goals

- 一句結論：5 大目標定義 ISM 學術貢獻；目標 1 為基礎（subclone detection），其餘 4 目標依賴目標 1
- 適用對象：論文 writing、研究優先級決策
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  grep -A 30 "五個研究目標定錨" /bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/MEMORY.md 2>/dev/null
  ```

---

## 五大目標（優先順序）

### 🎯 目標 1：Subclonal Methylation Detection（基礎）
**定義**：從 long-read 甲基化數據中準確偵測 tumor subclonal structure

**成功標準**：
- HPFineNGroups ≥4 region TP rate 顯著高於 baseline（已驗證 +3.7pp）
- 跨 7 樣本一致性驗證

**狀態**：✅ POSITIVE characterization（但非 variant filter）

---

### 🎯 目標 2：Epigenetic Landscape Characterization
**定義**：繪製 tumor 甲基化異質性全景

**成功標準**：
- Zone-Aware framework 可描述 5 個 zone
- LOH × AF × Methylation 三維關係（已驗證：Inter AF→NGroups +0.705）

**狀態**：✅ POSITIVE（依賴目標 1）

---

### 🎯 目標 3：Cross-Region Subclonal Order
**定義**：跨 region 的 two-hit order 分析（哪個 subclonal event 先發生？）

**成功標準**：
- Region-level temporal ordering 可行
- P4 pilot 驗證

**狀態**：⏸ **暫緩**（P4 pilot: region-level 無 two-hit order 訊號；待目標 1 更成熟；MEMORY: `project_p4_target3_depends_on_target1`）

---

### 🎯 目標 4：Variant Filter Support
**定義**：ISM 作為 variant filter 提升 F1

**成功標準**：
- ΔF1 > +0.01 跨樣本一致

**狀態**：⚠️ **部分完成**（paired_full +0.0112 locked；TO -0.0206；Phase 1A 凍結，不再追求）

---

### 🎯 目標 5：Normal Methylation Reference Framework
**定義**：整合 normal BAM + CN/purity correction 的 subclone baseline

**成功標準**：
- Normal baseline 讓 subclone detection 更穩
- Purity-aware correction 可應用

**狀態**：🚀 **Phase 2 進行中**（HCC1395 pilot POSITIVE 97.3% sig；待 7 樣本全量驗證）

---

## 目標關係圖

```
       ┌─────────────────────────────────────┐
       │ 目標 1 Subclonal Detection（基礎）   │
       └────────────┬────────────────────────┘
                    │
        ┌───────────┼───────────┬──────────────┐
        ▼           ▼           ▼              ▼
    目標 2        目標 3      目標 4         目標 5
    Landscape     Order       Filter         Normal Ref
    ✅ POSITIVE  ⏸ 暫緩      ⚠ Locked      🚀 進行中
```

---

## 主軸優先順序

1. **目標 1** — 基礎（已有 HPFineNGroups + LOH×AF positive findings）
2. **目標 5** — 目前最高優先（Phase 2 A+D 方向）
3. **目標 2** — 與目標 5 並行（Zone-Aware characterization）
4. **目標 3** — 暫緩，待目標 1 更成熟
5. **目標 4** — 已 locked，不再追求（characterization-only）

---

## MEMORY 對應

- `project_research_vision_five_goals`
- `project_p4_target3_depends_on_target1`（目標 3 暫緩）

---

## 相關

- 專案總覽：[01_project_summary.md](01_project_summary.md)
- 突破策略：[04_breakthrough_strategy.md](04_breakthrough_strategy.md)
- 當前 phase：[03_current_phase.md](03_current_phase.md)
