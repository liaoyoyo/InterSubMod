# Style Library Case Studies — v1/v2/v3 物件對照

本檔記錄歷史簡報的物件 / 版型套用案例，作為新 slide 設計的參考。
對應 playbook §19（v3 24-slide 物件抽取來源）+ §22（v3 主軸聚焦 audit 示範）。

---

## 何時用何種版型（Decision Tree）

```
是否量化 Before/After 對照？
├── 是 → before_after_split
└── 否
    ├── 是否多方法 / 多樣本對比？
    │   ├── 2 方/樣 → side_by_side_compare
    │   ├── 3 方 → 3card_tldr 或 5card_grid
    │   ├── 4 方 → 4quadrant_matrix（若交叉維度）
    │   └── 5 方 → 5card_grid
    │
    ├── 是否流程展示？
    │   ├── 線性 3-5 階段 → pipeline_flowchart
    │   ├── 時序 / 版本歷程 → timeline_horizontal
    │   ├── 分支邏輯 → decision_tree_layout
    │   └── 根因分析 → root_cause_tree_with_trigger
    │
    ├── 是否資料展示？
    │   ├── 單主圖 + caveat → data_main_with_caveat
    │   ├── 多 KPI 並列 → kpi_dashboard
    │   └── 矩陣 / 混淆矩陣 → 4quadrant_matrix
    │
    ├── 是否 executive summary？ → 3card_tldr
    │
    └── 是否生態 / 領域定位？ → family_tree_with_2x2
```

---

## v3 24-slide 物件對照（playbook §19 七案例）

### 案例 1: S1 amber 動機 strip

**原 slide 描述**：
- v3 第 1 張（cover 後第 1 張內容）
- 主訊息：「Self-Phasing 17.3:1 artifact 是當前所有下游 ISM 統計失效的根因」
- 視覺：頂部琥珀色 strip + 大字 17.3:1

**套用的物件 + 版型**：
- 版型：`data_main_with_caveat`（簡化版，主視覺為大字數字）
- 主物件：`motivation_amber_strip`（頂部，動機句）
- 輔物件：`tldr_big_number`（中央 17.3:1）
- 標題：`thesis_title_bar`（「Self-Phasing 是當前最大 artifact」）

**為何選此**：
- 動機 slide 需要琥珀色強調「為何在意」
- 大字數字 17.3:1 立刻吸引注意，比文字更具衝擊力
- 不需 caveat（cover 後第一張不適合放 limitation）

**變體建議**：
- 若加入 sample-level 17.3:1 分布圖 → 升級為 `data_main_with_caveat` 完整版
- 若有多個動機並列（罕見）→ 改用 `3card_tldr`（amber tone 變體）

---

### 案例 2: S11 紅色根框 + 3 葉樹

**原 slide 描述**：
- v3 第 11 張（root cause 章節）
- 主訊息：「Self-Phasing 17.3:1 是 LOH region somatic-first phasing 觸發」
- 視覺：紅色 root box（17.3:1）→ 3 個藍色 branch（PON-only / V3-Fixed / V5）→ 3 個 leaf（驗證結果）

**套用的物件 + 版型**：
- 版型：`root_cause_tree_with_trigger`
- Root 物件：`decision_tree_node`（type: root, RED tone）+ 內嵌 `tldr_big_number`（17.3:1）
- Branch 物件：`decision_tree_node` × 3（type: branch, blue）
- Leaf 物件：`decision_tree_node` × 3（type: leaf, green for V5、grey for V3-Fixed superseded）
- 底部：`caveat_red_strip`（「若不解決 → 下游 ISM 失效」）

**為何選此**：
- Root cause 分析的標準版型，紅色 root 強調觸發點
- 3 分支 / 3 leaf 是根因 + 解決方案 fit-to-screen 的最佳數量
- 底部紅 caveat 與 root 紅 box 形成上下「問題語氣」呼應

**變體建議**：
- 若 leaf 為演算法決策（非根因驗證）→ 改用 `decision_tree_layout`（root 中性）
- 若分支 >3 → 拆 slide（一次太多分支違反 Rule 3 ≤6 elements）

---

### 案例 3: S13 雙欄程式碼 diff

**原 slide 描述**：
- v3 第 13 張（method 章節）
- 主訊息：「V5 三層投票將 AMB% 從 17.5 降至 8.0」
- 視覺：左欄 baseline getvote 程式碼 / 右欄 V5 getvote 程式碼 / 下方 Δ AMB strip

**套用的物件 + 版型**：
- 版型：`before_after_split`
- 左欄：`code_diff_box`（baseline，紅色 - 行強調）
- 右欄：`code_diff_box`（V5，綠色 + 行強調）
- 底部 delta_strip：`tldr_big_number`（「Δ AMB 17.5% → 8.0%」）
- 標題：`thesis_title_bar`（「V5 三層投票將 AMB% 從 17.5 降至 8.0」）

**為何選此**：
- 程式碼 diff 是 before/after 的具體實作 evidence
- before_after_split 的 delta_strip 是 focal point 焦點區
- 雙欄等寬讓兩段程式碼對齊比較

**變體建議**：
- 若程式碼 >15 行 → 拆 slide 或改為 `data_main_with_caveat`（單欄主圖 + caveat）
- 若不是程式碼而是演算法 pseudocode → 仍用 `before_after_split` + `method_blue_box` 替代 `code_diff_box`

---

### 案例 4: S15 主圖 + 數據 caveat

**原 slide 描述**：
- v3 第 15 張（result 章節）
- 主訊息：「reviewer 校準口徑：來源報告 38%，實為 14/85=16.5%」
- 視覺：中央 KDE 分布圖 / 右側 reviewer 校準說明 / 下方紅色 caveat strip

**套用的物件 + 版型**：
- 版型：`data_main_with_caveat`
- 主視覺：`data_main_panel`（KDE 分布）
- 右側 annotation：`insight_green_card` + `method_blue_box` + `neutral_grey_caveat`
- 底部：`caveat_red_strip`（「⚠ 來源報告寫 38%, 實為 14/85=16.5%」）
- 右下角：`citation_footnote`

**為何選此**：
- 資料 slide 標準版型
- reviewer 校準口徑必須用紅色 caveat 顯著標註
- 右側 annotation 區可分層（insight / method / scope）

**變體建議**：
- 若無 reviewer 校準需求 → 改為 `neutral_grey_caveat`（中性語氣）
- 若多圖並列 → 改用 `4quadrant_matrix`

---

### 案例 5: S20 業界家族樹

**原 slide 描述**：
- v3 第 20 張（discussion / related work 章節）
- 主訊息：「InterSubMod 在 long-read 表觀遺傳工具家族中佔據獨特位置（高甲基化深度 + 高變異整合）」
- 視覺：上半部家族樹（按工具類別分支）+ 下半部 2×2 矩陣（本研究位於右上象限）

**套用的物件 + 版型**：
- 版型：`family_tree_with_2x2`
- 家族樹：`decision_tree_node` × 多個（root + 3 branch + 6-8 leaf）
- 2×2 矩陣：mini `4quadrant_matrix` 變體
- 右側 callout：`insight_green_card`（本研究定位 highlight）+ `method_blue_box`
- 底部：`citation_footnote`（多工具 paper references）

**為何選此**：
- Related work slide 需同時呈現「生態」+「定位」
- 家族樹給 context、矩陣給 differentiation
- focal point 在 2×2 右上象限（本研究位置）

**變體建議**：
- 若僅需家族樹（無定位）→ 改用 `decision_tree_layout`
- 若僅需定位（無生態）→ 改用 `4quadrant_matrix`
- 若工具 >10 → 拆 slide 或進 backup

---

### 案例 6: S22 五大目標卡片

**原 slide 描述**：
- v3 第 22 張（conclusion / future work 章節）
- 主訊息：「InterSubMod 五大研究目標：特徵化 / 亞克隆偵測 / 甲基化-變異關聯 / LOH 整合 / 跨樣本」
- 視覺：5 cell 並列卡片，每 cell 含 icon + 目標名 + 簡短說明 + 狀態 badge

**套用的物件 + 版型**：
- 版型：`5card_grid`
- 5 cell：`side_card_5cell`（單一物件含 5 cell array）
- 底部：`tldr_big_number` mini（「主軸：目標 1 + 目標 4」）
- 右下：`citation_footnote`

**為何選此**：
- 5 cell 是 5 大目標的最佳數量（≤6 elements 規則邊界）
- focus_card_index 可指定當前主軸（目標 1 + 4）highlight
- 底部 summary strip 強調主軸選擇

**變體建議**：
- 若 4 目標 → 改用 `4quadrant_matrix`（若交叉維度）或 `3card_tldr`（並列）
- 若每目標需詳細說明 → 拆多 slide，每 slide 一目標 + `data_main_with_caveat`

---

### 案例 7: S24 Take-home + Next + Q&A

**原 slide 描述**：
- v3 第 24 張（最後 1 張，conclusion）
- 主訊息：「Take-home: V5 完成；Next: Phase 2 Normal BAM；Q&A」
- 視覺：3 張卡片並列 / 底部 contact 資訊

**套用的物件 + 版型**：
- 版型：`3card_tldr`
- 左卡：`insight_green_card`（Take-home 結論）
- 中卡：`motivation_amber_strip`（Next step / call-to-action）變體（升級為 card 格式）
- 右卡：`method_blue_box`（Q&A / contact info）
- 底部：`citation_footnote`（commit / data version）

**為何選此**：
- 結論 slide 的標準三段式：What we did / What's next / Open questions
- 三色語氣（綠 / 琥珀 / 藍）對應三種訊息類型
- 視覺權重平均，不偏向任何單一卡片

**變體建議**：
- 若無 Q&A 需求 → 改用 `before_after_split`（左 take-home / 右 next）
- 若 take-home 含多重結論 → 拆 conclusion slide + 獨立 next slide

---

## 物件未對應 v3 的補充說明

以下物件在 playbook §19 v3 抽取列表中未顯式對應，但仍納入樣式庫：

| 物件 | 補充原因 |
|------|---------|
| `kpi_card` | v3 未使用，但 executive summary template 必需 |
| `pipeline_stage_box` | v3 S2 用於 4 階段 pipeline，與 method_blue_box 區別需獨立 |
| `focal_point_marker` | meta 元素，storyboard 階段強制標註，非最終 slide 視覺 |
| `data_main_panel` | v3 S15 主視覺底層元件，常被嵌入其他物件，獨立列出方便調用 |

---

## 物件 / 版型互斥規則速查

| 場景 | 應用 | 不應用 |
|------|------|--------|
| Before/After 同對象量化對照 | `before_after_split` | `side_by_side_compare` |
| 兩方法 / 兩樣本對比 | `side_by_side_compare` | `before_after_split` |
| 線性 3-5 階段流程 | `pipeline_flowchart` | `timeline_horizontal` |
| 時序 / 版本歷程 | `timeline_horizontal` | `pipeline_flowchart` |
| 演算法決策樹（focal: leaf） | `decision_tree_layout` | `root_cause_tree_with_trigger` |
| 根因分析（focal: root） | `root_cause_tree_with_trigger` | `decision_tree_layout` |
| 強警示 / reviewer 校準 | `caveat_red_strip` | `neutral_grey_caveat` |
| 中性限制聲明 | `neutral_grey_caveat` | `caveat_red_strip` |
| 單一 hero 數字 | `tldr_big_number` | `kpi_card` |
| 多 KPI dashboard | `kpi_card` × N + `kpi_dashboard` layout | `tldr_big_number` |

---

## 新案例新增 SOP（playbook §18）

新增案例至本檔時：
1. 確認新 slide 的物件 / 版型已存在於樣式庫
2. 若不存在 → 先寫 YAML（依 §13 schema）
3. 寫案例段落（原描述 / 套用 / 為何 / 變體）
4. 加入 v1/v2/v3 編號（v?-S?）
5. 更新「物件未對應」表格（如有新增 v? 對應）
