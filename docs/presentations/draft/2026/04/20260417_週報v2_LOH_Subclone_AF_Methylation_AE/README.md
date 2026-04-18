<!--
建立時間: 2026-04-17 05:10
目標: 週報 PPT v2 — Assertion-Evidence 精簡版的設計理念與使用說明
處理範圍: 16 投影片 + 5 重製圖 + Okabe-Ito 色盲友善配色
關聯檔案:
  - theme_okabe_ito.py
  - figure_specs.py
  - build_pptx_v2.py
  - ../../../validated/2026/04/20260414_研究週報_LOH_Subclone_AF_Methylation/build_pptx.py (v1 參考)
-->

# 週報 PPT v2 — Assertion-Evidence 精簡版

## 定位

此為 `20260414_研究週報_LOH_Subclone_AF_Methylation/` 的**完全重新設計版本**，並存於 `draft/`，讓用戶對照兩版後決定採用。

## 為何重做（與 v1 差異）

| 面向 | v1 (validated) | v2 (draft) |
|------|----------------|------------|
| 投影片數 | 41 張 | **16 張** |
| 標題風格 | 名詞片語（「LOH 定義與 Subclone 分層」） | **完整 claim 句**（「LOH 區域 NGroups 每樣本皆顯著較高」） |
| 配色 | 暖色系 (#C76B50 赭紅 + #6E9D8A 薄荷) | **Okabe-Ito 色盲友善 8 色** |
| 圖表 | 部分手繪 rect/text、部分引用 png | **5 張關鍵圖全部 matplotlib 重製 300 DPI** |
| 程式碼 | 3460 行 / 41 函數 | **~1200 行 / 16 函數 + 8 helpers** |
| 排版 | 元素密度高，手動調整 | **統一模板 15/70/15 垂直切分** |

## 設計參考

1. **Michael Alley, Penn State** — *The Craft of Scientific Presentations* — Assertion-Evidence 方法
2. **Masataka Okabe & Kei Ito** — 色盲安全 8 色調色盤（<https://jfly.uni-koeln.de/color/>）
3. **Edward Tufte** — *The Visual Display of Quantitative Information* — 資料墨水比最大化
4. **Nature Methods Points of View** — 科學視覺化專欄（2010-2013）

## 16 投影片大綱

| # | Claim Title (zh) | Evidence |
|---|------------------|----------|
| 1 | LOH 區域發現 subclonal 表觀遺傳結構 | 封面 + #2 縮圖 |
| 2 | 三個核心問題定義本週進度 | 3-col Q1/Q2/Q3 |
| 3 | 中間 AF 區間富含 LOH | fig 01 |
| 4 | LOH 區域 NGroups 每樣本顯著較高 | fig 02 |
| 5 | AlleleDelta 在 LOH 呈大效應分離 | fig 03 |
| 6 | Segment-level spatial correlation 中等 | fig 04 |
| 7 | 四分層策略分離 high-TP 亞群 | fig 05 |
| 8 | HPFineNGroups 作為 somatic heterogeneity 標記 | schematic |
| 9 | 雙證據鏈 AF × methylation | flow diagram |
| 10 | 研究定位：read-level epigenetic context | 對比表 |
| 11 | 已關閉方向 prevent re-investigation | 2-col |
| 12 | Phase 2 方向 A+D 為最高優先 | roadmap |
| 13 | 限制：7 樣本 × 有限 coverage | caveats |
| 14 | 下週 3 個可驗證交付物 | checklist |
| 15 | 問題與討論 | — |
| 16 | Appendix 方法與資料連結 | refs |

## 快速使用

```bash
cd docs/presentations/draft/2026/04/20260417_週報v2_LOH_Subclone_AF_Methylation_AE/

# 步驟 1：重新產出 5 張關鍵圖 (300 DPI)
python3 regenerate_figures.py

# 步驟 2：生成 PPT
python3 build_pptx_v2.py

# 輸出：outputs/20260417_週報v2_AE.pptx + _snapshot.json
```

## 驗證清單

- [ ] `figures/` 有 5 張 PNG，每張 ≥ 300 DPI
- [ ] `outputs/20260417_週報v2_AE.pptx` 含恰 16 張投影片
- [ ] 每張標題為完整 claim 句（含主詞 + 動詞 + 結論）
- [ ] 配色只使用 Okabe-Ito 8 色 + 中性灰
- [ ] 雙語：中文 28pt + 英文 17pt (60% 灰)
- [ ] Color-blind simulator (Coblis) 三種色盲模式皆可辨識

## 廢棄條件

若用戶視覺審查後仍傾向 v1，此資料夾可安全刪除；v1 完整保留於 `validated/` 下不受影響。
