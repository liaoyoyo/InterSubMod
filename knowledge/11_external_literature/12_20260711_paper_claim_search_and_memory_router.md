---
id: ism-kb-11-paper-claim-search-memory-router-20260711
name: "Paper claim × external evidence search and memory router"
description: "以 stable claim IDs、paper sections、evidence roles 與 claim status 路由107張逐卡curated外部卡；納入2026-07-12 producer closeout與downstream pending邊界。"
status: active
last_verified: 2026-07-12
content_nature: reference-summary
doc_type: reference
verified_scope: "external_validation 107-card manifest + 234 semantic links + 134-artifact strict audit + InterSubMod 20260712 raw-all producer closeout；clean 7/7 downstream rates pending"
decision_refs: [ISM-A01, ISM-M02, ISM-M03, ISM-R01, ISM-R03, ISM-D01]
related_ids: [ism-kb-11-external-literature-index, ism-kb-11-layered-reconstruction-external-bridge-20260710, ism-kb-11-latest-research-theory-contract-delta-20260711]
tags: [external-literature, claim-registry, semantic-search, memory, invalidation, longphase-s, regional-tree]
canonical_paths: [11_external_literature/12_20260711_paper_claim_search_and_memory_router.md]
alias_paths: []
---

# Paper claim × 外部證據搜尋與記憶路由

> **2026-07-12 live state**：107 cards / 106 entity keys / 234 source→claim links / 107 curated / 0 axis-default；readiness=`direct84/qualified19/background4`，strict audit 0 errors/0 warnings。完整方法裁決先讀[14 full method/process external validation](14_20260712_full_method_process_external_validation.md)。

## 一句結論

論文外部知識不再靠自由文字 `maps_to_thesis` 或舊 R#/V# 搜尋，而是先以 stable `ISM-*` claim ID 判斷主張狀態，再從 97-card semantic catalog 找 support/calibration/tension/limitation 來源，最後回讀 CONTEXT 與一手 paper/repo。catalog 是路由層，不是可引用來源。

2026-07-11 07:00 contract audit 是目前最高優先 override：LongPhase-S 會 rescue input non-PASS，因此完整 recalibration input 已改為 normalized paired ClairS raw-all；ClairS PASS 只作 sensitivity/audit，canonical tree input 是同一 run 的 `_sc.vcf` PASS。PASS-only producer 在 0/7 完成時中止並登錄 `ISM-X04`。HCC1954 chr22 與 patched HCC1395 whole-chrX bounded probes 通過，但 7/7 raw-all closeout 仍 pending。歷史 tagged BAM另只有 COLO829 genome-wide、6/7 受 `--truth-bed` 限制；舊 downstream results 只屬 engineering baseline，定量率由 `ISM-R03=pending_rerun` 阻擋。

外部證據庫目前有 **107 cards、20 claims、234 source-to-claim links**；**107/107 cards均為curated，axis default=0**。citation readiness是檢索優先級：direct84、qualified19、background-only4。strict live/manifest/catalog/link/coverage/artifact audit為0 errors/0 warnings。任何來源即使是direct，正文引用前仍需回讀該卡§3/§4/provenance與原始來源。

## 路由層

| 層 | 入口 | 功能 |
|---|---|---|
| Claim state | `/big7_disk/liaoyoyo2001/external_validation/_schema/paper_claims.tsv` | active、pending、retired、gate、禁止說法 |
| Semantic catalog | `/big7_disk/liaoyoyo2001/external_validation/_landscape/paper_evidence_catalog.tsv` | card 級 section/role/readiness |
| Normalized links | `/big7_disk/liaoyoyo2001/external_validation/_landscape/paper_evidence_links.tsv` | source × claim × relation |
| Source truth | `/big7_disk/liaoyoyo2001/external_validation/axis*/<slug>/CONTEXT.md` | 外部來源內容與 verification depth |
| Project truth | audit/DAG/v2 gate 文件 | 內部結果、失效與 validation scope |

## 術語與角色

| 術語 | 中文意義 | 研究角色 |
|---|---|---|
| `active_method` | 方法/邊界可寫 | Methods、Discussion 承重 |
| `active_result` | 目前內部 audit 結果可寫 | 僅限 canonical wording 與 scope |
| `pending_rerun` | 語意已定但 production 結果未過 gate | 不可變成 Results 百分比 |
| `retired` | 已被新 evidence 推翻 | anti-regression memory |
| `supports` | 直接支持一般方法定義/原理 | 不自動支持我方數值 |
| `calibrates` | 校準術語、estimand、分母、regime | 常用於 Related Works/Discussion |
| `tensions` | 不同 regime 成功，限制我方 wording | 不是同口徑 conflict |
| `complements` | 可作其他 modality/正交驗證 | future validation target |
| `citation readiness` | 文獻驗證深度的檢索分層 | 不取代逐句 source verification |

## 現行 claim 狀態速查

- 可寫方法：`ISM-A01`、`ISM-M01/M02/M04/M05/M06/M07`。
- 可寫 audit：`ISM-R01` record continuity + 1/7 historical genome-wide；`ISM-R02` solver consistency，必標非 biological truth。
- 可寫邊界：`ISM-R04/R05/R06`、`ISM-D01/D02`。
- 暫停：`ISM-M03` clean tag/sidecar contract、`ISM-R03` downstream rates。
- 淘汰：`ISM-X01` 23,810/35,332 backbone、`ISM-X02` 舊 7/7 comprehensive PASS、`ISM-X03` methylation filter/rescue/tree-ranker、`ISM-X04` ClairS-PASS-only LongPhase-S input contract。

## 內部搜尋

```bash
cd /big7_disk/liaoyoyo2001/external_validation
python3 _schema/query_library.py --list-claims --claim-status pending_rerun
python3 _schema/query_library.py --claim ISM-M02 --relation calibrates --tier A
python3 _schema/query_library.py --claim ISM-D01 --relation tensions,complements --mapping-confidence curated
python3 _schema/query_library.py --section S2.4 --readiness direct --format paths
```

寫作流程固定為：claim ID → status/gate → support/calibration → limitation/tension → CONTEXT/source → sentence with boundary。找不到 claim ID 的句子先不寫，避免舊 narrative 重新產生新 claim。

## 長期記憶邊界

應記：穩定定義、決策、失效來源、claim ID、gate、DOI/PMID/commit/hash、查詢入口。

不應記：running job 進度/ETA、尚未 atomic-complete 的 outputs、pending downstream rates、search snippet、整段 PDF 摘錄、由方法相似性推測 biological validity。若查到 axis-default mapping，應視為 schema regression 而非可保存知識。

更新時不刪舊 claim：先標 `pending_rerun` 或 `retired` 並記 invalidating audit；新 gate 真正完成後再更新同一 claim ID 的 status 與 `last_reviewed`。這讓後續 session 能知道「為何不能再用」而不是只看到新數字。

## 已知前提與未驗證界線

| 已驗證 | 尚未驗證 |
|---|---|
| 7/7 historical PASS-subset continuity；raw-all chr22/whole-chrX probes record/payload守恆 | normalized raw-all 7/7 fail-closed producer/closeout與 sidecar consumer U0–U7 完整通過 |
| 6/7 historical truth-BED scope mismatch | corrected funnel/determinacy/multi-HP final rates |
| HP 是 allelic family；H3?/H4/unphased 非 primary lineage | region-tree node/edge 對應 biological clone |
| solver engineering consistency | patient generalisation、second-caller robustness |
| CN missing=unavailable；read AF≠CCF；methylation 不 rank/confirm | purity/CN/multiplicity-corrected CCF 與 orthogonal clone truth |

## 建議閱讀順序

1. [7/11 contract audit](/big7_disk/liaoyoyo2001/InterSubMod/research/20260710_layered_reconstruction_v2/20260711_ClairS_LongPhaseS_sSNV位點與ReadTag契約稽核_01.md)
2. [Layered v2 remediation gates](/big7_disk/liaoyoyo2001/InterSubMod/research/20260710_layered_reconstruction_v2/00_INDEX.md)
3. [Claim-boundary DAG](/big7_disk/liaoyoyo2001/InterSubMod/docs/concepts/DAG/20260710_layered_reconstruction_claim_boundaries_01.md)
4. [Full search/memory guide](/big7_disk/liaoyoyo2001/external_validation/_landscape/10_paper_claim_evidence_search_and_memory_20260711.md)
5. [Paper claim registry](/big7_disk/liaoyoyo2001/external_validation/_schema/paper_claims.tsv)
6. [97-card semantic catalog](/big7_disk/liaoyoyo2001/external_validation/_landscape/paper_evidence_catalog.tsv)
7. [Strict semantic/provenance audit](/big7_disk/liaoyoyo2001/external_validation/_provenance/LIBRARY_SEMANTIC_AUDIT_20260711.md)

> 任何與本頁衝突的 2026-07-10 downstream percentage，先降為 engineering baseline；以 `ISM-R03` gate 判斷能否恢復。
