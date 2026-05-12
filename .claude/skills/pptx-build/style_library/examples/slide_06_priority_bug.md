# Slide 06 套用案例 — tagging 層 priority bug (per-read getVote)

> 2026-05-12 完成。使用 `upper_lower_main_dual.yaml` 佈局 + A1 + B2 + C2 變體組合。

## 套用組合

| 區域 | 物件 | 變體 | 來源 yaml |
|---|---|---|---|
| 上半 | igv_read_block | A1 Multi-read | `InterSubMod/.claude/skills/pptx-build/style_library/objects/igv_read_block.yaml` |
| 下半左 | code_diff_pane | B2 Compact | `InterSubMod/.claude/skills/pptx-build/style_library/objects/code_diff_pane.yaml` |
| 下半左前 | code_source_header | — | `InterSubMod/.claude/skills/pptx-build/style_library/objects/code_source_header.yaml` |
| 下半右 | flow_diagram_box | C2 Single mermaid | `InterSubMod/.claude/skills/pptx-build/style_library/objects/flow_diagram_box.yaml` |
| 上半子件 | position_header × 4 | germline_het + somatic_mut | `InterSubMod/.claude/skills/pptx-build/style_library/objects/position_header.yaml` |

## 資料填入（per-read 框架）

**位點配置**（4 個 variant 位點）：
- pos 1: chr19 某 germline het（PoN ✓）— marker = HP2 alt blue
- pos 2: chr19 某 germline het（PoN ✓）— marker = HP2 alt blue
- pos 3: chr19 某 somatic mut（PoN ✗）— marker = HP1_1 yellow w/ brown border
- pos 4: chr19 某 germline het（PoN ✓）— marker = HP2 alt blue

**3 條 read 各自覆蓋 subset**（不同 read 經過位點不同）：
| read | 經過位點 | countMap | baseline hp | V5/V6 hp |
|---|---|---|---|---|
| read 1 | pos 1, 2, 3 | HP2=2, HP1_1=1 | **11 ❌** | **21 ✅** |
| read 2 | pos 1, 2, 4 | HP2=3 (no somatic) | 2 ✓ | 2 ✓（同） |
| read 3 | pos 3, 4 | HP2=1, HP1_1=1 | **11 ❌** | **21 ✅** |

**程式碼 diff（B2 Compact）**：
- path: `longphase-to/HaplotagProcess.cpp:506-530` → `longphase-to-mod/HaplotagProcess.cpp:512-560`
- 內容：baseline iterate variantKeys vector + break early；V5/V6 重寫為 two-layer

**流程圖（C2 Single mermaid）**：
- 起點：`this read's countMap`
- Layer 1: germline 決方向
- Layer 2: somatic 加 sub-tag
- 終點：`hp = 21 ✅`

## 為什麼這組合

1. **主視覺 A1 Multi-read** — 凸顯 per-read scope（每條 read 獨立 tag）+ 為何 752 victims 不是 aggregate
2. **下半左 B2 Compact** — 程式碼依據（與圖文並列不喧賓奪主）
3. **下半右 C2 Single** — 只強調 V5 修法（baseline 錯誤已在 A1 右側 panel 顯示 hp=11 ❌）
4. **不選 B1 Full / C1 Dual** — slide 06 已有 A1 + 4-grid 基本概念，再放雙完整對照會超過 ≤6 element 紅旗

## 邏輯元素計算

| 邏輯元素 | 含子件 |
|---|---|
| 1. IGV view (A1) | 4 pos + 3 read + 3 tag panel + observation footer = ~30 子件 |
| 2. Code diff (B2) | path inline + 7 -/+ lines = ~10 子件 |
| 3. Flow diagram (C2) | 5 mermaid nodes = ~8 子件 |

**總計 3 邏輯元素 ≤ 6 紅旗 ✓**

## Speaker note 重點（per-read 框架）

- getVote 是 per-read 操作（C++ HaplotagProcess.cpp:533 每條 read reset countMap）
- 每條 read 各自的 countMap → 各自的 getVote → 各自的 HP tag
- baseline iterate variantKeys vector：① somatic pair → break / ② mixed / ③ germline 永遠輪不到
- V5/V6 重寫為 two-layer：Layer 1 germline 決方向 / Layer 2 somatic 加 sub-tag encoding
- read 1/read 3 受 priority bug 影響（有 somatic + germline mix）；read 2 純 germline 不受影響
- 752 victims = 752 條符合「somatic + germline mix」條件的獨立 read，非 aggregate

## Cross-link 到 slide 16

A1 右側 panel 顯示 V5/V6 hp=21 ✅，**但** germline-absent 區域 V5 Layer 1.5 仍有 caveat（priority bug 4.19:1 偏 HP1）— 詳 slide 16

## 驗證 checklist

- [x] per-read 框架明示（A1 紅框警告 + speaker note）
- [x] 4 個位點 / 3 條 read 物件齊全
- [x] 每 read tag panel 含 baseline + V5 對照
- [x] B2 path 標明 baseline + V5/V6 來源
- [x] C2 只顯示修法側（綠系）
- [x] 顏色語義一致（baseline 紅 / V5 綠 / germline 藍 / somatic 黃）
- [x] ≤ 6 element 紅旗 ✓（3 邏輯元素）

## 後續可改進

- 加 hover tooltip 看每條 read 完整 cigar / variant detail
- 連接到 IGV 真實截圖（slide 04a/04b 已有 SP1/2/3）
- backup slide 用 A1 + B1 + C1 完整版（appendix）
