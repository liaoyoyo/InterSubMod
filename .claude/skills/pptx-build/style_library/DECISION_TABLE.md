# 模板取用決策表 — 關鍵字 → A/B/C 變體選擇

> 2026-05-12 起。Self-Phasing PPT 系列用。當收到 PPT slide 需求時，依關鍵字查此表取對應變體 + 物件數量配置。

## 查表規則

```
收到關鍵字 → 查決策表 → 取 A/B/C 變體 yaml → 填 slide 數據 → ≤ 6 element 紅旗檢查 → 套入 build_html.py
```

## 主決策表

| 關鍵字 / slide 目的 | A 變體 | B 變體 | C 變體 | 總物件（邏輯元素）| 建議佈局 |
|---|---|---|---|---|---|
| **priority bug 機制完整講解（slide 06）** | A1 Multi-read | B2 Compact | C2 Single | ~45（3 邏輯） | `upper_lower_main_dual.yaml` |
| **baseline ↔ V5 雙向深入對照（backup）** | A2 Single-read | B1 Full | C1 Dual | ~50 | 左 50% A2+B1 / 右 50% C1 |
| **per-read scope 反駁 aggregate 誤解** | A2 Single-read | B3 Inline | C3 Linear | ~18 | A2 主軸 + B3/C3 註腳 |
| **HP tag 概念介紹（教學早期 slide）** | A3 Minimal | — | C3 Linear | ~10 | A3 視覺 + C3 流程線 |
| **752 victims 數據展示（slide 08）** | A3 Minimal | — | — | ~6 | A3 引用 + 表格主體 |
| **5-commit fix timeline（slide 07/10）** | — | B2 per commit | C3 時序 | ~15 | B2 列 + C3 timeline |
| **球員兼裁判 phasing 層（slide 05）** | A3 Minimal（引用 somatic edge）| — | C3 Linear（3-step）| ~12 | 概念 step 圖 + A3 補充 |
| **backup / appendix 深入解釋** | A1 Multi-read | B1 Full | C1 Dual | ~60 | 單張 full 或拆 2 張 |

## ≤ 6 element 紅旗檢查

playbook §20.E：每張 slide ≤ 6 視覺元素。違規時：

1. **合併子組件計為 1 element**（如 4 個 position headers + 1 read + 4 markers + countMap 合算 1 個 "IGV view"）
2. **拆分到下一張 slide**（如 backup B1/C1）
3. **降級變體**（A1→A2 / B1→B2 / C1→C2）

## 物件數量速查

| 變體 | 物件數 | 適用區域 |
|---|---|---|
| A1 Multi-read | ~30 | 單張上半 |
| A2 Single-read | ~12 | 半張 |
| A3 Minimal | ~6 | 內嵌 |
| B1 Full | (大，~30 行 diff) | 單張左半 |
| B2 Compact | (中，~10 行 diff) | 半張 |
| B3 Inline | (小，~3 行 highlight) | inline |
| C1 Dual mermaid | (大，~12 nodes × 2) | 單張右半 |
| C2 Single mermaid | (中，~5 nodes) | 半張 |
| C3 Linear HTML | (小，~5 div boxes) | inline |

## 顏色語義約定（避免混用）

- **baseline** → 紅系（`#FEE2E2` bg / `#DC2626` stroke / `#7F1D1D` text）
- **V5/V6 修法** → 綠系（`#DCFCE7` bg / `#16A34A` stroke / `#166534` text）
- **Layer 1 germline** → 綠系
- **Layer 2 somatic encoding** → 藍系（`#DBEAFE` / `#1E3A8A`）
- **PoN ✓ germline het** → 藍系
- **PoN ✗ somatic** → 黃 + 棕邊（`#FBBF24` / `#B45309`）
- **caveat / warning** → 黃系（`#FEF9C3` / `#CA8A04`）

## 取用範例（slide 06 完整應用）

```yaml
slide: 06
target: priority bug 機制完整講解
layout: upper_lower_main_dual.yaml
components:
  upper:
    object: igv_read_block.yaml
    variant: A1_multi_read
    data:
      positions: [pos1_germ_het, pos2_germ_het, pos3_somatic, pos4_germ_het]
      reads: [read1, read2, read3]  # 不同覆蓋 subset
      countMap_per_read: [...]
      tag_calc: [baseline_hp_per_read, V5_hp_per_read]
  lower_left:
    object: code_diff_pane.yaml
    variant: B2_compact
    data:
      path_baseline: longphase-to/HaplotagProcess.cpp:506-530
      path_V5: longphase-to-mod/HaplotagProcess.cpp:512-560
      diff_lines: [...]  # key -/+ only
  lower_right:
    object: flow_diagram_box.yaml
    variant: C2_single
    data:
      mermaid: "flowchart TD\n  Start --> L1 --> L2 --> R"
emphasis: V5 修法 + Layer 1.5 caveat 跨 slide 16
```

## 未來擴充

- 新關鍵字 / 新 slide 主題 → 加入主決策表新 row
- 新變體（如 A4 cross-region / B4 multi-commit）→ 對應 yaml 新增 + 在此表加 row
- 衝突案例 → 加入「失敗模式」section
