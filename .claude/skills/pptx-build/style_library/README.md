# myPPT Style Library

樣式調用資料夾索引。所有 PPTX 物件、版型、顏色 token 集中於此，由 `ppt_toolkit.style_library_loader` 反序列化為 Python dataclass 後，於 `build_pptx.py` 中以 `load_object()` / `load_layout()` API 取用。

對應 playbook 章節：§12（Style Library 結構）、§13（物件 YAML schema）、§14（顏色 palette）、§15（通用版型分類）、§16（對齊規則）、§17（toolkit 調用 API）、§19（v3 24-slide 物件抽取來源）、§23（tier_recommendation + focal_point_zone）。

---

## 目錄結構

```
style_library/
├── README.md                # 本檔（索引與調用約定）
├── colors/
│   └── palette.yaml         # 6 色 token + use_for/do_not_use_for + WCAG + colorblind
├── objects/                 # 15 個標準物件 YAML
├── layouts/                 # 12 個通用版型 YAML
└── examples/
    ├── wireframes/          # 由 Stream D 產生 PNG 樣板（本 stream 不負責）
    └── case_studies.md      # 何時用何種版型 + v1/v2/v3 案例對照
```

---

## 物件清單（15 個）

| 檔名 | category | tier_recommendation | 主要用途 |
|------|----------|---------------------|---------|
| `caveat_red_strip.yaml` | warning | A | 標 reviewer 校準 / 來源錯誤 / 過度推論警示 |
| `insight_green_card.yaml` | conclusion | S | 標核心發現 / 改進結論 |
| `method_blue_box.yaml` | method | S | 標方法論 / 演算法定義 / pipeline 階段 |
| `motivation_amber_strip.yaml` | motivation | A | 標研究動機 / 為何在意 |
| `neutral_grey_caveat.yaml` | caveat_neutral | B | 中性語氣的限制聲明（非 reviewer 警示） |
| `tldr_big_number.yaml` | summary | S | 大字 delta 數字（17.3:1、F1+0.0112、Δ AMB 17.5→8.0%） |
| `thesis_title_bar.yaml` | title | S | Assertion-evidence 標題列（每張 slide 強制） |
| `citation_footnote.yaml` | reference | B | 右下角引文 / 數據來源 |
| `code_diff_box.yaml` | data | A | 雙欄程式碼 diff（baseline vs fixed） |
| `kpi_card.yaml` | summary | S | 單一 KPI 數字 + label + delta 三層卡片 |
| `decision_tree_node.yaml` | structure | S | 決策樹節點（root / branch / leaf） |
| `pipeline_stage_box.yaml` | structure | S | Pipeline 階段方塊（含階段名 + 輸入 / 輸出） |
| `side_card_5cell.yaml` | structure | S | 5-cell 並列卡片（用於 5 大目標、5 階段 audit） |
| `focal_point_marker.yaml` | meta | A | 視覺焦點標記（Storyboard 階段標 ≤ 20 字 focal point） |
| `data_main_panel.yaml` | data | S | 主數據面板（大圖 + 軸標 + 圖例） |

---

## 版型清單（12 個）

| 檔名 | template_type | focal_point_zone 位置 |
|------|---------------|----------------------|
| `before_after_split.yaml` | improvement, comparison | delta_strip (bottom center) |
| `4quadrant_matrix.yaml` | comparison, data_showcase | center cross-axis |
| `3card_tldr.yaml` | executive_summary | center 3-card row |
| `pipeline_flowchart.yaml` | improvement, concept | flow center horizontal axis |
| `side_by_side_compare.yaml` | comparison | metric delta strip (bottom) |
| `timeline_horizontal.yaml` | improvement, concept | end-of-timeline state |
| `kpi_dashboard.yaml` | executive_summary, data_showcase | top-left primary KPI |
| `decision_tree_layout.yaml` | concept, academic_defense | leaf node decision |
| `data_main_with_caveat.yaml` | data_showcase, academic_defense | main_panel 中央 |
| `root_cause_tree_with_trigger.yaml` | improvement, academic_defense | trigger root box |
| `family_tree_with_2x2.yaml` | concept, academic_defense | 2×2 quadrant intersection |
| `5card_grid.yaml` | concept, executive_summary | center card (target #1) |

---

## 顏色 token

詳見 `colors/palette.yaml`。6 色：red / green / blue / amber / grey / black。每色含：
- `hex` — 16 進位色碼
- `use_for` — 推薦語意用途清單
- `do_not_use_for` — 反向清單
- WCAG AA 對比規則
- 色盲安全（紅綠不相鄰、形狀 + 顏色雙編碼）

---

## 調用約定（API 範例）

```python
from ppt_toolkit.style_library_loader import load_object, load_layout

# 載入物件規格
caveat = load_object("caveat_red_strip")
caveat.render(slide, x=1.0, y=5.5,
              content="⚠ 來源報告寫 38%, 實為 14/85=16.5%")

# 載入版型骨架
layout = load_layout("before_after_split")
layout.populate(slide, {
    "title": "Baseline tag voting vs V5 三層投票",
    "left.heading": "Baseline (somatic-first)",
    "left.main_visual": "figures/baseline_getvote.png",
    "right.heading": "V5 (germline-first + Layer 1.5)",
    "right.main_visual": "figures/v5_getvote.png",
    "delta_strip.content": "Δ AMB 17.5% → 8.0%",
})
```

實作上：YAML 規格在 build 時 deserialize 成 Python dataclass，配合 `python-pptx` 渲染。所有對齊與顏色都從 YAML 讀取，**不再 inline hardcode**。

---

## 命名約定

- **物件 / 版型**：小寫底線（`caveat_red_strip.yaml`）
- **YAML 內 `name` 欄位**：與檔名一致，可用連字號（`caveat-red-strip`）
- **顏色 token**：英文常用名（red / green / blue / amber / grey / black）
- **欄位名**：snake_case（`tier_recommendation`、`focal_point_zone`）

---

## 新增物件 / 版型 SOP（playbook §18）

1. 確認需求：哪個 template_type、解決什麼版面問題
2. 找 example：v1/v2/v3 舊 slide 是否已有原型
3. 寫 YAML：依 §13 schema 填欄
4. 寫例子：good + bad 各 2-3 個
5. 跑 wireframe 渲染：產 PNG 進 `examples/wireframes/`
6. 進 `examples/case_studies.md`：補何時用、何時不用
7. PR review：是否與既有物件重疊、命名是否一致、是否符合對齊規則

---

## 與其他 stream 的邊界

- 本 stream（Stream B）只負責：YAML schema + 物件 + 版型 + colors + case_studies.md
- Stream A：playbook、SKILL.md、templates、prompts
- Stream C：ppt_toolkit Python 模組 + audit
- Stream D：`examples/wireframes/` PNG 圖片渲染（依賴本 stream 完成）
