---
id: ism-kb-11-layered-reconstruction-external-bridge-20260710
name: "2026-07-10 Layered reconstruction × external evidence bridge"
description: "將 paired ClairS、SEQC2、LongPhase-S/TO、perfect-phylogeny/CN/CCF 與 methylation regime boundaries 對齊 InterSubMod 七月區域 mutation-state tree 主線。"
status: active
last_verified: 2026-07-11
content_nature: reference-summary
doc_type: explanation
verified_scope: "InterSubMod 20260711 ClairS/LongPhase-S contract audit + external_validation CONTEXT/manifest；7/10 downstream rates 改為 engineering baseline。"
decision_refs: []
related_ids: [ism-kb-11-external-literature-index, ism-kb-11-paper-claim-search-memory-router-20260711, ism-kb-11-latest-research-theory-contract-delta-20260711]
tags: [external-literature, paired-clairs, seqc2, longphase, regional-tree, perfect-phylogeny, cn, ccf, methylation-boundary]
canonical_paths: [11_external_literature/11_20260710_layered_reconstruction_external_bridge.md]
alias_paths: []
---

# 2026-07-10 layered reconstruction × 外部證據橋接

> 🔴 **2026-07-11 07:00 OVERRIDE**：先讀 [12_20260711_paper_claim_search_and_memory_router.md](12_20260711_paper_claim_search_and_memory_router.md)。完整 LongPhase-S recalibration input 已改為 normalized paired ClairS raw-all；ClairS PASS 只作 sensitivity/audit，canonical tree input 是同一 run 的 `_sc.vcf` PASS。PASS-only producer 已在 0/7 完成時中止。兩個 raw-all bounded probes 通過，但 7/7 closeout 未完成；歷史 tags 6/7 另受 truth BED 限制，舊 downstream rates繼續暫停。

## 一句結論

現行 InterSubMod 以 **normalized paired ClairS raw-all → LongPhase-S `_sc` PASS tree input** 的同分子共現為骨幹，ClairS PASS 另作 backbone sensitivity；在每個 germline-HP family 內枚舉/拒絕區域 mutation-state compatible trees。CN 是後驗解讀，methylation 是 bounded auxiliary。region tree 能支持局部突變狀態推論，但在沒有 single-cell/multi-region orthogonal truth 前，不等於 biological clone tree。

## 知識對照

| 問題 | 現行答案 | 外部錨點 |
|---|---|---|
| paired caller 是什麼？ | ClairS paired v0.4.0 production；不可寫成 ClairS-TO 或最新 v0.5.0 | `clairs-paired-2026` |
| PASS 是真值嗎？ | 不是；raw-all 是 recalibration universe、ClairS PASS 是 sensitivity subset、LongPhase-S `_sc` PASS 才是 tree input | `clairs-paired-2026` + `seqc2-hcc1395-somatic-truth-2021` |
| 7/7 downstream 結果可寫嗎？ | 尚不可；舊 tagging 只有 1/7 genome-wide，clean rerun pending | 20260711 contract audit + `ISM-R03` |
| truth 能驗 tree 嗎？ | SEQC2 只驗 variant-level，不能驗 clone/tree | `seqc2-hcc1395-somatic-truth-2021` |
| HP family 是 clone 嗎？ | 不是；它是 germline allelic partition/discriminator | `longphase-s-somatic-haplotagging-2025` |
| tumor-only 名詞怎麼寫？ | LongPhase-TO 應寫 tumor DNA fraction；legacy `tumor_purity` 不等於 cellular purity | `longphase-to-2025` |
| 多解如何報？ | 枚舉相容最小樹；無唯一性就報 ambiguity，不任選 | Gusfield / Pe'er / SPRUCE / minimum-uncovering cards |
| read AF 是 CCF 嗎？ | 不是；缺 purity/CN/multiplicity 時不可替代 | ASCAT / DeCiFer / Tarabichi cards |
| 缺 CN 怎麼辦？ | unavailable，不可預設 neutral | ASCAT / DeCiFer；InterSubMod 2026-07-10 P0 修正 |
| methylation 能 rank tree 嗎？ | 目前不能；只作 bounded annotation/characterization | Sgootr/Gaiti/Epiclomal/MethylTree/PCDH regime boundary |

## 投稿紅線

- 可寫：`regional mutation-state compatible tree inferred from long-read sSNV co-occurrence`。
- 不可寫：`biologically validated subclone tree`、`HP family equals clone`、`read AF equals CCF`。
- 可寫：`methylation was evaluated as a bounded auxiliary layer`。
- 不可寫：`methylation improves caller F1 / ranks lineage trees`，除非另有新正交證據。
- paired ClairS、SEQC2、LongPhase-S 的版本、reference、chemistry、model、BED/VCF checksum 必個別記錄。

## 最新入口

- 研究 SoT：`/big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md`
- 分層單位/分母：`/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md`
- 外部完整橋接：`/big7_disk/liaoyoyo2001/external_validation/_landscape/09_intersubmod_research_bridge_20260710.md`
- 機械 catalog：`/big7_disk/liaoyoyo2001/external_validation/_landscape/library_manifest.tsv`
- 全庫稽核：`/big7_disk/liaoyoyo2001/external_validation/_provenance/LIBRARY_AUDIT_20260710.md`
- Claim/search/memory：`/big7_disk/liaoyoyo2001/external_validation/_landscape/10_paper_claim_evidence_search_and_memory_20260711.md`

> Root gate：current research tree 是 `/big7_disk/liaoyoyo2001/InterSubMod`；`/big8_disk/...` 是三月快照，不可用來覆蓋七月數字與主張。
