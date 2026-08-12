<!--
建立時間: 2026-08-12
目標: 用原始論文與官方頁面限定 InterSubMod 公開文字中關於 short read、single-bulk、lineage 與 methylation 的外部主張
處理範圍: 僅針對公開文字中可被外部文獻反駁或限縮的一般性命題；不以外部 paper-scope 效能替代本專案驗證
關聯檔案:
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv
  - InterSubMod/docs/reports/validated/2026/08/20260812_InterSubMod_GitHub公開說明與教學完整驗證_01.md
-->

# 外部原始來源與可支持邊界

| 公開主張 | 原始來源的直接證據 | 對 InterSubMod 文字的影響 | 不能反向推論 |
|---|---|---|---|
| short reads 只有 marginal VAF | [PairClone](https://academic.oup.com/jrsssc/article/68/3/705/7058355) 明確建模同一 short read 跨過的 proximal mutation pairs，用於 subclone inference | 這個絕對句被反駁；應改為「長 reads 提供更遠距離、更密集的 joint linkage」 | 不代表 short reads 與 long reads 有相同可連結距離或性能 |
| single bulk 原理上只能 characterize、不能作任何 subclonal reconstruction | [Nature Methods practical guide](https://pmc.ncbi.nlm.nih.gov/articles/PMC7867630/) 指出 single 或 multiple bulk samples 皆可進行模型化重建，但依賴 purity、copy number、multiplicity 與模型假設 | 通用性「不可能」應改為本專案的 evidence ceiling：「現有 frozen analysis 未確認 cellular clones」 | 不代表 single-sample reconstruction 無偏、唯一或等同真值 |
| single-sample 重建只是不可驗證概念 | [2024 DREAM benchmark](https://www.nature.com/articles/s41587-024-02250-y) 以 31 個演算法、51 個有真值的模擬腫瘤與 12,061 runs 量化 single-sample reconstruction | 說明這是可 benchmark 的推論問題，不是無法進行的問題；同時強烈支持「演算法、purity-adjusted depth、CN 與 mappability 影響結果」 | 模擬 truth 效能不能直接證明 InterSubMod 的真實樣本 accuracy |
| methylation 不能用於腫瘤演化重建 | [Mazor et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC4573399/) 在多空間／多時間 glioma 樣本中獨立建立 phyloepigenetic 與 mutation phylogenetic trees，並報告兩者高度關聯 | 全領域絕對句被反駁；應限縮為「本專案的 single-bulk read-level methylation 因 ASM/LOH/CN 混雜，目前不作 backbone」 | 多區域樣本的群體平均甲基樹，不等於單一 bulk read cluster 可直接定義 cellular clone |
| 同分子共現等於直接看到 cellular lineage | PairClone/TreeClone 將 same-read mutation pairs 作為輸入，但 clone number、genotype、frequency 與 tree 仍是模型推論 | 可以說「直接觀測單分子共現」，不得跳成「直接觀測細胞譜系」 | 不否定 same-read linkage 比 marginal-only 多一層可辨識資訊 |

## 整合裁決

- 外部文獻反駁的是「全領域的絕對說法」，不是 InterSubMod 對自己 frozen evidence 採保守上限的正當性。
- 現行可防守文案應同時保留兩件事：同分子 joint evidence 比 marginal VAF 多資訊；但 InterSubMod 現有資料仍只支持 local molecular candidate structure，不確認 cellular lineage。
