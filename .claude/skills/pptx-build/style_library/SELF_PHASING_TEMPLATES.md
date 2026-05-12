# Self-Phasing PPT 模板索引（A/B/C 三模板 × 3 變體）

> 2026-05-12 起新增。專門用於 longphase priority bug / per-read getVote / phasing 球員兼裁判等概念視覺化。
> 既有 README.md 列 15 個通用物件；本檔列 self-phasing 專用模板（5 物件 + 1 layout + 決策表）。

## 模板清單

| 檔案 | 類別 | 用途 |
|------|------|------|
| `objects/igv_read_block.yaml` | A 模板 | IGV-style read 視覺化（pentagon + position marker） |
| `objects/code_diff_pane.yaml` | B 模板 | GitHub +/- diff 風格程式碼差異 |
| `objects/code_source_header.yaml` | B 配件 | 檔案路徑/行號標題框 |
| `objects/flow_diagram_box.yaml` | C 模板 | mermaid / linear 流程圖 |
| `objects/position_header.yaml` | A 配件 | 位點 header (PoN ✓ / PoN ✗) |
| `layouts/upper_lower_main_dual.yaml` | layout | slide 06 用：上半 A1 / 下半左 B2 / 下半右 C2 |
| `DECISION_TABLE.md` | 決策表 | 關鍵字 → A/B/C 變體選擇查表 |
| `examples/slide_06_priority_bug.md` | 案例 | slide 06 套用 A1+B2+C2 的完整紀錄 |

## 變體規格速查

每個模板有 3 變體（依使用情境選擇）：

| 模板 | 變體 1 (Full) | 變體 2 (Compact) | 變體 3 (Inline) |
|------|---------------|-----------------|----------------|
| **A igv_read_block** | A1 Multi-read（3+ reads 各自 tag panel）| A2 Single-read（1 read scope 證明）| A3 Minimal（1 read + 2 pos 示意）|
| **B code_diff_pane** | B1 Full（雙 header + 完整 +/- ）| B2 Compact（inline path + 關鍵 hunk）| B3 Inline（1 行 path + 高亮 2-5 行）|
| **C flow_diagram_box** | C1 Dual mermaid（baseline + V5 並排）| C2 Single mermaid（只強調修法）| C3 Linear（純 HTML div+arrow）|

## 取用流程

```
收到關鍵字 → 查 DECISION_TABLE.md → 取 A/B/C 變體 yaml → 填數據 → ≤ 6 element 紅旗檢查 → 套入 build_html.py
```

## 覆蓋優先級（與既有 playbook 一致）

1. inline (build_html.py style attr) ≈ inline style
2. layouts/*.yaml 特殊規則 ≈ #id selector
3. objects/*.yaml 預設值 ≈ .class selector
4. colors/palette.yaml token ≈ :root vars

## 後續擴充項

- A 模板可慢慢補：A4 cross-region paired/TO 對照、A5 LOH BED 區段標記、A6 methylation track 等
- B 模板可補：B4 git log -p 多 commit 連續變動
- C 模板可補：C4 sequence diagram（call stack）、C5 state machine

## 驗證來源

- C++ 對證：`/big8_disk/liaoyoyo2001/longphase-to/HaplotagProcess.cpp:506-530`（baseline）+ `/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp:512-560`（V5/V6）
- per-read concept memory：`memory/project_getvote_per_read_concept.md`
- source 行號 memory：`memory/reference_longphase_getvote_source.md`
- demo HTML：`InterSubMod/docs/presentations/validated/2026/05/self_phasing_synthesis_PI/preview/igv_read_template_demo_v2.html`
