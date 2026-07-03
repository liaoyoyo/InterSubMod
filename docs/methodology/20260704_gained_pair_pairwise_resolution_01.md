---
title: gained-pair pairwise 定序 — 「缺中間群等機率」的修正 + 區域內所有子read provenance
date: 2026-07-04
sample: HCC1395 (frozen, SEQC2 truth)
status: in_progress
tier: L2（單樣本 read-level 實證 + VAF 獨立佐證；⭐3 characterization cap）
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/resolution_subreads.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/candidate_trees.json, docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json
build_branch: research/subclonal-reconstruction-202606
---

# gained-pair pairwise 定序 + 區域內所有子read provenance

## 1. Context（為何做）
使用者 catch 到 `enumerate_candidate_trees.py` 的一個盲點：它把「缺中間群」(A_ambiguous_order) 區標成 **等機率候選**（37 區全 `equiprobable`），並宣稱 read 無法排序（我先前照 production 現況說「0/37 可排」）。但**枚舉器只用「全跨 read 的 population(節點)」判斷，沒查那兩個缺失中間群突變的 pairwise 2×2** —— 而該 2×2 其實已含定序資訊。

## 2. 演算法（gained-pair pairwise 定序，3 層判定）
對 tree 每條 gain ≥2 突變的邊，查 gained 突變間的 pairwise 2×2（cells RR/RA/AR/AA, coread）：

| gained pair 關係 | 判定 | 意義 |
|---|---|---|
| co_linked（RA=AR=0, AA≥2）| **BLOCK** | 兩突變永遠一起 → 無中間群，記「一起獲得」一個 event，不排序 |
| nested（小邊 <ε=2%）| **ORDERED** | 中間群被觀測到（大邊）→ 順序可定；方向由大邊決定 |
| 4-gamete（雙邊 ≥ε）| **CONFLICT** | 真衝突，不成單樹 |
| 無 coread（>50kb/未共讀）| **AMBIG** | 才是真等機率/欠定 |

**ε=2%** 對齊 ONT 錯誤率（小邊 <2% coread = 噪聲可去）；nested 方向再用 **pos_vaf 交叉核對**（較早突變應較高 VAF=較 clonal）。

## 3. 完整結果（HCC1395 3885 有向量區）

| 區級分類 | 數 |
|---|---:|
| clean 單突變邊（本就清楚）| 3845 |
| **BLOCK 共event** | **15** |
| **ORDERED 可定序**（solid 5 + tentative 19）| **24** |
| CONFLICT 真衝突 | 1 |
| AMBIG_NOCOREAD（真等機率）| **0** |

**修正「0/37 可排」**：production 原標 37 缺中間群等機率 → 實為 **15 block + 24 可定序 + ~1 真衝突、0 真等機率**。等機率標籤完全是「沒查 gained-pair pairwise」的產物（L2）。

**驗證證據**：
- ✅ **VAF 交叉核對 24/24 全一致**（0 不一致）— 定序方向獨立佐證
- ✅ **transitivity 0 違反**（nested 無成環）
- ✅ **chr17:48357368-48365161 重現 canonical**（S2/S4 co_linked → BLOCK）
- ✅ ε-tiering 對齊 production（raw CONFLICT 8 → ε 後 1）
- ✅ truncated 42 區反查：全區衝突已由 determinacy 捕捉（15 incompatible + 34 collapse-single）→ per-edge 層無漏判

## 4. 盲點與無法判斷的邊緣案例（4 類，誠實標記）
1. **block@CN-gain（12）** — co_linked 是真共event 或 CN-multiplicity 假共現？需 allelic CN（SAVANA）才分得清。
2. **ε-tentative ordered（20）** — 中間群證據 <2% coread → 可能 ONT 噪聲；需「≥2 read 全面觀察統計層」+ 噪聲 caveat 才浮出。
3. **conflict-drop@CN-gain（4）** — 4-gamete 小邊 <ε 被降級 ordered，但 CN-gain 區小邊可能是真 multiplicity → 🔴 保留 conflict-or-order，不直接判 ordered。
4. **truncated/dense（42）** — 非 resolver 盲點（determinacy 已捕捉）。

→ 3 個盲點裡 2 個根因是缺 allelic CN → 再次指向 **SAVANA 為最高槓桿**。

## 5. 區域內所有子read（provenance 輸出）
每區產出所有子read 的 genotype 向量（含 X=未覆蓋）+ 覆蓋分布 + full-span 絕對群 + pairwise 2×2 + 邊解析 + **絕對/組合 provenance**（絕對＝單分子跨全部；組合＝pairwise 推得）。42 個多突變區完整落 `resolution_subreads.json`。範例 chr17：絕對比 31%（3 群 32 read 單分子跨全 4 點），18 種子read 向量。

## 6. 工作站顯示
`build_topology_workstation.py` 加 `resolveBlock()` 面板：點多突變區 → detail 底部顯示解析分類 + pairwise 2×2 + 區域內所有子read + provenance + 盲點旗標。HTML: `20260629_multisample_topology_workstation.standalone.html`（HCC1395 分頁，node --check + resolveBlock 實跑驗證通過）。

## 7. 邊界（守牆）
- ⭐3 characterization cap（單樣本、single-pipeline）；影響 ~40 區（3885 的 ~1%），headline determinacy 44.8% 幾乎不動，但「缺中間群」敘述被實質修正。
- 定序純用 read 共現 + VAF（CN-clean）；**絕不用甲基**（循環）。
- assembled/composed 永不升 A_determined；block@CN-gain / conflict-drop@CN-gain 需 CN 才定案。

## 8. 產物
- 資料：`data/resolution_subreads.json`（42 區）、`data/resolve_v2_complete.json`（全 3885 分類，scratchpad→待入庫）
- 腳本：`scripts/resolve_v2_complete.py`（resolver）、`scripts/gen_resolution_subreads.py`（子read+provenance 生成）
- 下一步：擴全 3885 區子read（per-chrom 平行）／落地 `enumerate_candidate_trees.py`（定序邏輯正式接入）／roll-out 7 樣本。
