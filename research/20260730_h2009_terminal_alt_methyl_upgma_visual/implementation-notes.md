<!--
建立時間: 2026-07-30
目標: 紀錄 H2009 chr18:563687 ALT-read 甲基距離與 UPGMA HTML 的實作決策、偏離、折衷與未決事項
處理範圍: 單一位點 Task F 視覺化；不修改 production caller 或分析方法
關聯檔案: InterSubMod/research/20260730_h2009_terminal_alt_methyl_upgma_visual/pre-decision-audit.md
-->

# Implementation notes：H2009 terminal-ALT methyl UPGMA visual

## 設計決定

- 統計與主圖固定使用全部 86 個 after-peel focal-ALT reads。
- UPGMA 採與 formal M1 一致的 Bernoulli distance、`max(D,Dᵀ)` 對稱化與 average linkage。
- 樹與 86×86 熱圖共用同一個 canonical leaf order。
- HTML 以 canonical portable artifact 打包；原生 chart 顯示 10-read deterministic display subset，完整 86-read 樹＋矩陣以可重現 PNG 嵌入。

## 偏離

- Portable artifact 沒有 native dendrogram，因此 UPGMA 不能用原生 chart；改用 sandboxed data-URI PNG。
- 現成 region `tree.nwk` 是 124 ALT+REF reads，與 formal M1 範圍不同，明確不重用。

## 折衷

- 完整 86×86 tidy matrix 有 7,396 rows，超過單 dataset 2,000-row contract；主圖仍顯示全部 86，native heatmap 改顯示 5 個大群 medoid-nearest reads＋全部 5 個小群 reads，使 desktop/mobile 不需捲動即可看完整 10×10。
- read 名稱以 aliases 與 read-name SHA8 呈現，避免在報告中暴露完整 raw UUID。

## 未決

- 若後續要將 methyl group 升級為 clone/subclone，仍需 matched-normal、CN/purity/multiplicity/CCF 與獨立細胞層級驗證。
- 本 Task F 不估計跨樣本 prevalence，也不寫 evidence ledger。
