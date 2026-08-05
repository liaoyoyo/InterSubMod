<!--
建立時間: 2026-07-12
目標: 紀錄 HCC1395 跨來源 VAF／clone／區域 read／拓撲精確整合驗證的設計決定與未決限制
處理範圍: chr1-22；HCC1395 與 HCC1395_DORADO；固定 5,720 個 historical layered-v2 complete-both regions
關聯檔案:
  - InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/pre-decision-audit.md
  - InterSubMod/research/20260712_hcc1395_pair_site_topology_containment_validation/data/
  - InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/report/artifact.json
狀態: validated_with_claim_ceiling
-->

# Implementation Notes：HCC1395 位點、區域、clone 與 topology 精確整合驗證

## 設計決定

- [決策] 任務類型為 B（comprehensive validation），固定 chr1-22 與完整 5,720-region population；所有 conditional rate 同時保留 fixed denominator。
- [決策] 證據依解析度分成 raw exact-locus VAF、PyClone-VI separate fits、shared-site relation、full-region candidate forest、selected shape；禁止跨層把相似度直接當 accuracy。
- [決策] C 表示 read-supported mutation-state groups／Steiner-tree nodes，不是 PyClone clone count；T 表示樹結構；Topo 表示 rooted-unlabeled shape。
- [決策] 外部 clone 主要 endpoint 只使用兩來源獨立 separate fits；joint fit 只作條件模型描述。
- [決策] 子結構一致只接受 genomic shared-site induced projection 或已凍結 candidate-set containment；不接受任意 existential subtree 配對。

## 偏離

- [偏離] clean layered-v3 全量重建仍在 RUNNING；本輪 regional topology 只能使用明確標示的 historical layered-v2 frozen population，待 v3 final receipt 後以相同腳本重算。

## 折衷

- [折衷] site-level binomial-noise 模型可量化 read-count sampling，但不能取代 BAM split-half／depth-matched topology rerun；因此「抽樣造成區域拓撲差異」只能列為部分支持，不作單一根因。
- [折衷] DORADO 暫共用 HCC1395 SAVANA allele-specific CN 與 purity=0.96；PyClone 結果屬 shared-CN sensitivity，不是 source-independent CCF validation。

## 未決

- [未決] DORADO source-specific allele-CN／purity、兩來源 BAM union-site recount、within-source split-half/downsample 與 single-cell/colony/synthetic truth 尚缺。
- [已解決] PyClone mutation key 與 5,720 regions 的 join coverage 已完成：雙側 joined 14,369/15,713；cluster/CP compatibility、uninformative 與 conflict 已分欄輸出。
- [未決] clean layered-v3 完成後，需更新區域分母、k strata、topology ladder 與 HTML evidence release。

## 執行結果

### Regional precision

- [結果] 固定 5,720 regions 全量 join 與 32/32 守恆檢查 PASS；k=1/2/3/>=4 為 0/3,121/1,506/1,093。
- [結果] Site-set equal 5,534/5,720（96.75%）；shared allele identity 15,713/15,713。
- [結果] Coarse agreement 3,969/5,720（69.39%），但 k>=4 的 83.81% 主要由 unresolved↔unresolved 888/1,093 堆高。
- [結果] Read strict+true-induced 1,599/5,720；VAF strict+true-induced 1,790/5,720；VAF conflict 1,234/3,860 evaluable。

### Clone-to-region bridge

- [結果] 14,369/15,713 regional alleles 可接到兩來源 separate PyClone fits；regional subclonal Jaccard=0.239、kappa=0.379，低於 global 0.381/0.544。
- [結果] 5,028/5,189 表面 cluster partition exact 中，5,007 為 both-single-cluster；真正 both-multicluster 為 21/34 exact（61.76%）。
- [結果] Assignment>=0.8 subset 的 perfect concordance 只有 1 個 two-source subclonal-union mutation，判為 selection-induced degeneracy。
- [結果] 598 個兩側都有 singleton read direction 的 pairs 中 596 同向、2 反向；但 high-confidence CP edges 全為 same-cluster，無不同 cluster ordering 資訊。

### Cause decomposition

- [結果] 16/16 守恆檢查 PASS；baseline raw VAF exact join 15,713/15,713，latest cross-release 為 15,636/15,713。
- [結果] k>=5 vs k=2 的 read/VAF conditional exact-or-induced 差 -38.62/-57.73 pp；read candidate >20 vs 1 差 -73.57 pp，VAF candidate >=2 vs 1 差 -46.71 pp。
- [結果] Global 93.07% binomial-compatible 與 1.089x noise 只支持 locus-level sampling contribution；depth imbalance 接近零，未建立 topology sampling causality。
- [結果] 703/5,720 HP count mismatch 令細粒度 endpoint fail closed；最終原因分類為 identifiability-dominant、pipeline-contributed、sampling-plausible、biological-divergence-unresolved。

### 最終交付與 QA

- [結果] Portable HTML 已輸出至 `InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/report/20260712_HCC1395跨來源VAF_CCF與subclone驗證_01.html`；SHA256=`8b407e9395e2f59af730bf3ac232c0dadb9b716d70dded87ba2732d5d1cc1c77`。
- [結果] Portable package verification PASS；獨立 Playwright desktop/mobile 26/26 PASS，console/page error=0。
- [結果] Evidence ledger precision revision 已寫入；正式 claim ceiling 為 partial technical reproducibility，true-tree accuracy NO-GO。
