---
id: ism-kb-11-latest-research-theory-contract-delta-20260711
name: "2026-07-11 latest literature, theory terminology, and raw-all contract delta"
description: "濃縮 Rosenski、McAuley、PCDH、Blood、Foltz 最新查證，rooted three-gamete 術語更正與 ClairS raw-all→LongPhase-S 契約，供研究理解與 anti-regression。"
status: active
last_verified: 2026-07-12
content_nature: reference-summary
doc_type: reference
verified_scope: "primary VOR/preprint/repo facts plus InterSubMod 20260712 seven-dataset raw-all producer/receipt closeout; fresh layered downstream remains pending"
decision_refs: [ISM-M02, ISM-M03, ISM-R01, ISM-R03, ISM-M06, ISM-R05, ISM-X04]
related_ids: [ism-kb-11-layered-reconstruction-external-bridge-20260710, ism-kb-11-paper-claim-search-memory-router-20260711]
tags: [external-literature, latest-research, rooted-three-gamete, raw-all, longphase-s, methylation-boundary, anti-regression]
canonical_paths: [11_external_literature/13_20260711_latest_research_theory_contract_delta.md]
alias_paths: []
---

# 最新文獻、理論術語與 raw-all contract 增量

> **2026-07-12 override**：本檔是7/11歷史增量。raw-all producer/receipt現已7/7 PASS；fresh layered downstream仍未7/7。Blood review已取得並親讀VOR全文；ASMS已升peer-reviewed 2026 VOR。現況與完整方法裁決請以[14 full method/process external validation](14_20260712_full_method_process_external_validation.md)為準。

## 一句結論

最新外部研究讓 methylation 的 lineage potential 與 confounding 都更清楚，但沒有替 InterSubMod 補上 biological-clone truth；最新內部 contract 已由 ClairS-PASS-only 更正為 **normalized paired ClairS raw-all → LongPhase-S `_sc.vcf` PASS**。7/7 producer已完成，clean downstream仍未完成。

## 本輪五個知識更新

1. **Rosenski et al. 2026**：官方 publication date 2026-04-15；39 種 purified cell types 的 WGBS 顯示 33,574 個 sequence-dependent ASM regions；是 germline/cell-type confound 的最新正式錨點，不是 cancer somatic lineage evidence。Code/Source Data/GEO 可公開定位，但 human EGA 個體資料是 controlled access。
2. **McAuley et al.** 已於 2026-05-26 正式發表於 Scientific Reports；仍是單病人 pre/post-treatment 跨平台 pilot，不是 regional-tree benchmark。
3. **PCDH preprint** 全文證實特殊 methylation barcode 可在 PD6646 ONT WGS 回復 clade fractions並對照既有 single-cell-derived-colony WGS tree；但無 PERMANOVA、branch lengths 不穩、不可外推任意 locus。
4. **Blood 2026 review** 是 CC BY-NC-ND OA；2026-07-12已透過ScienceDirect VOR親讀全文，卡片升`fulltext_verified`。全文將lineage epimutations與cell identity/transcription/environment難以區分列為outstanding question。
5. **Foltz SomaticHaplotype** 有正式公開 repo；linked-read somatic co-occurrence 是 genetic-spine prior art，但 shared barcode/haplotype 不等於完整 clone tree。

## 理論術語

Gusfield 1991 的 directed all-zero-root gain-only model 檢查 `10/01/11` 三種 non-root patterns，應稱 **rooted three-gamete condition**；若明示 virtual `00` root，也可稱 **root-augmented four-gamete condition**。InterSubMod solver 原行為正確，本輪是名稱與 citation correction，不改演算法結果。

Complexity 也要分 regime：AncesTree 的 NP-complete 對象是 multi-sample VAF Factorization Problem，ancestry arborescence 另須 sum condition；Foulds–Graham 證的是 unrooted/undirected minimum-Hamming Steiner Problem in Phylogeny。兩者都不能被簡寫成「任一 rooted incompatibility 即 NP-hard」。

## Input contract 與 anti-regression

- `ISM-M02`：normalized paired ClairS raw-all biallelic sSNVs 是完整 recalibration universe；同一 run 的 `_sc PASS` 是 canonical tree input；ClairS PASS 只作 sensitivity/audit。
- `ISM-X04`：PASS-only producer 已淘汰；它在 0/7 完成時中止，因會讓 non-PASS→PASS rescue branch 永遠不可達。
- bounded evidence：HCC1954 chr22 426→426；patched HCC1395 chrX 35,823→35,823；兩者 record/payload/transition/sidecar gates 通過。
- 已完成：7/7 raw-all fail-closed producer/receipt closeout。尚未完成：fresh layered consumer、U3–U7 final verifier、sensitivity/comparison與immutable closeout；`ISM-M03/R03`保持pending。

## 寫作邊界

- 可寫 `regional mutation-state tree within parental haplotype families`。
- 不可寫 `biologically validated clone tree`。
- 可寫 methylation 在 single-cell/special-barcode regime 可攜帶 lineage information。
- 不可用該成功結果推論 bulk ONT methylation clusters 就是 subclones，或用 methylation rank trees/filter variants。
- 外部 paper 對內部 production/result claims 通常是 `calibrates/complements`，不是 `supports`。

## 權威入口

- 完整增量與數字：[external landscape 11](/big7_disk/liaoyoyo2001/external_validation/_landscape/11_latest_research_theory_and_contract_delta_20260711.md)
- Claim registry：[paper_claims.tsv](/big7_disk/liaoyoyo2001/external_validation/_schema/paper_claims.tsv)
- Current contract audit：[20260711 contract audit](/big7_disk/liaoyoyo2001/InterSubMod/research/20260710_layered_reconstruction_v2/20260711_ClairS_LongPhaseS_sSNV位點與ReadTag契約稽核_01.md)
- Search/memory router：[12_20260711_paper_claim_search_and_memory_router.md](12_20260711_paper_claim_search_and_memory_router.md)
