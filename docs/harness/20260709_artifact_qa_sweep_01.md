<!--
類型: 全面 artifact QA sweep 紀錄
建立: 2026-07-09
工具: tools/batch_render_qa.py (render→自動偵測) + Claude Read PNG 目視
build_branch: research/subclonal-reconstruction-202606
data_sources: scratchpad qa_batch/qa_summary.json
-->

# 交付物全面 render QA sweep（2026-07-09）

> **範圍**：docs/ 下 48 個活躍 HTML 交付物（layered_workstation 7 樣本 + topology_workstation 7 樣本 + standalones + phylo/subcluster pilots + methyl_lineage）。排除 66MB multisample 巨檔（已被 per-sample 版取代）。
> **方法**：`tools/batch_render_qa.py`（playwright render → pageerror/broken-image/blank 自動偵測）→ Claude Read 代表性 PNG 目視（tofu/溢出/layout）。= verify-workstation W7 閉環的批次版。

## 裁決

| 結果 | 數 | 說明 |
|---|---|---|
| ✅ OK（auto-check 過）| 47 | 無 pageerror、無 broken-image、非空白 |
| 🔴 BROKEN_IMG | 1 | `20260623_geometry_divergence_observation_dashboard.standalone`（見下）|
| 目視確認 | 4 | Q5 / H2009 methyl（7+10 圖 OK）+ 分層樹 flagship（07-07，完美）+ geometry(破) |

**flagship 分層樹工作站（07-07，14.3MB）目視 = 完美渲染**：U1-U6 層級表、HP-multiplicity 彩條、determinacy 徽章、比例字典、region 瀏覽器、CJK 全正常；數字 U1=23,810 somatic 與 §13.0 修正一致。

## 🔴 唯一缺陷（H-008）
`20260623_geometry_divergence_observation_dashboard.standalone.html`（07-02 舊 pilot）引用 `figs_geomdiv/geom_TP_*.png`（400 張），但 **`figs_geomdiv/` 目錄完全不存在** → **400/400 破圖**。「standalone」名不副實（圖沒內嵌也沒 bundle）。非活躍 flagship → 未修，待決定是否重生 400 圖。

## 工具自我修正（本 sweep 揪出 batch_render_qa 2 bug，已修）
- **H-009 檔名碰撞**：stem-only → layered/topology 同樣本名 PNG 互相覆蓋 → 改 parent-dir namespace。
- **H-010 lazy-load 漏抓**：首輪只報 28/400 破圖（fold 下 lazy img `complete==false`）→ 加逐屏 scroll 觸發 lazy-load → 正確 400/400。

## 結論
**整體工作流交付物健康**：47/48 活躍 HTML 渲染乾淨；唯一破的是一週前的舊 pilot（figs 目錄缺）。閉環 QA 工具本身也被此 sweep 磨出 2 個 bug 並修正 + 記 log（H-008~010），示範 self-harness 失敗記錄的閉環價值。
