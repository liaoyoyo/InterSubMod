<!--
建立時間: 2026-07-12 10:00
目標: HCC1395 與 HCC1395_DORADO 區域粗拓撲、癌症基因與藥物註記整合驗證前審計
處理範圍: chr1-22、全 7 dataset rows、同一 HCC1395 生物樣本的兩套技術處理
cycle_id: 20260712-1000-hcc1395-pair-coarse-topology-validation
topic: hcc1395_pair_coarse_topology_gene_drug_validation
status: verdict_GO_engineering_report_only
audit_version: 0.2
關聯檔案:
  - InterSubMod/docs/CURRENT_FOCUS.md
  - InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/regional_topology_composition.json
  - InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json
-->

# Pre-Decision Audit：HCC1395 兩流程粗拓撲與基因／藥物一致性

> **初始 verdict：PROBE（50/100）**。描述性統計可直接做，但「一致＝方法有效」涉及 region boundary、basecaller、CN、共同上游與外部註記等混淆；先以 HCC1395 區域配對與 null model 作 safe-to-fail probe。

## §0 Cynefin Domain Gate

- **Domain**：Complex。
- **Test**：相同分析流程尚未對同一 biological sample 的兩套 raw-all／basecalling 結果反覆產生可預測的 region-level topology agreement。
- **Rationale**：全域粗比例的相近不保證同一基因組區域相同；癌症基因／藥物富集也可能由 gene length、mutation density、CN 與註記密度造成。

## §1 Observation Completeness Checklist

| Observation | 狀態 | Tier | 來源 |
|---|---|---|---|
| HCC1395 與 HCC1395_DORADO 是同一 biological sample 的技術處理比較 | ✓ | L2 | `docs/research/methylation_methodology/2026/03/20260307_5kHz主實驗與方法學驗證藍圖_01.md` |
| Historical layered-v2 已有 7 dataset rows、exact T 與 regional topology census | ✓ | L2 | `research/20260711_read_group_C_tree_T_topology_report/data/` |
| 7/11 23:58 coarse-topology census 65/65 engineering checks PASS | ✓ | L2 | `docs/CURRENT_FOCUS.md` |
| Clean layered-v3 aggregate `_SUCCESS` 與 verification 在 7/11 23:58 尚未完成 | ✗ | L1 | `docs/CURRENT_FOCUS.md` |
| Exact/overlap-matched HCC1395 pair 的五類 topology agreement | □ | — | 本輪 probe |
| 癌症基因與藥物 annotation join coverage／重複鍵／版本 | □ | — | 本輪盤點 |
| 區域樹有 single-cell 或 multi-region biological ground truth | ✗ | L1 | `docs/CURRENT_FOCUS.md`；目前明示缺正交 ground truth |

## §2 Credibility Score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | 同一生物樣本技術重處理應保留穩定的大尺度訊號；樹形分類可明確定義。 |
| 觀察支撐 | 10 | 有全資料 engineering census，但尚無 region-matched pair concordance。 |
| 機制清晰度 | 10 | 有 read-state→Steiner T→coarse class 邏輯，但 boundary/basecaller/CN 會改變觀測。 |
| 反例風險 | 0 | 共同上游、非獨立樣本、region boundary、dominant class、CN/coverage 皆可能製造一致。 |
| 所需資源 | 10 | Python 全量 join／null／HTML 約 1–6 小時。 |
| **TOTAL** | **50 / 100** | **PROBE** |

**Falsifier observable**：若技術穩健性假設錯，exact 或 reciprocal-overlap matched regions 的五類 agreement 不會高於 chromosome-preserving permutation null；重要癌症基因區也不會比可比背景更一致，且差異會集中於低 read support／Topo>1／CN-altered 區。

**Reality-test 三個反例**：

1. 全域比例相似但 region-level Cohen's κ 接近 0。
2. 高 agreement 完全由「無 HP 內關係」dominant class 造成；排除該類後消失。
3. 癌症基因／可藥物基因的一致性不高於 gene-length／region-count matched background。

## §3 Assumption Map

| Assumption | Importance | Known | Quadrant |
|---|---|---|---|
| 五類在 mutation-bearing observed states 上可互斥且守恆 | HIGH | UNKNOWN | 2 — 必先驗 |
| 兩資料列代表同一 biological sample、不同技術處理 | HIGH | KNOWN | 1 |
| region boundary 可用 exact 與 reciprocal-overlap 兩口徑比較 | HIGH | UNKNOWN | 2 — 必先驗 |
| 癌症基因／藥物註記可作生物合理性支持但非 topology truth | HIGH | KNOWN | 1 |
| clean layered-v3 已完整可用 | HIGH | UNKNOWN | 2 — 必先查 |
| 外部註記版本不影響主要結論 | LOW | UNKNOWN | 4 — 記錄版本即可 |

## §4 Quick Pilot Guide

1. 讀 historical regional topology JSON 與最新 producer lifecycle，固定可用 snapshot。
2. 對 HCC1395 pair 建 exact-coordinate 與 reciprocal-overlap matched region 表，套五類 coarse topology。
3. 檢查 raw agreement、balanced accuracy、Cohen's κ、category Jaccard，並跑 chromosome-preserving permutation null。

**Checkpoint**：

- GO：matched n≥1,000，κ 高於 null 95% 上界，且排除最大類後仍高於 null；join 守恆 100%。
- PROBE：matched n<1,000 或只在一個 matching definition 顯著。
- NO-GO：κ 不高於 null，或五類定義不能守恆。

## §5 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| clean layered-v3 terminal status／aggregate verification | HIGH | 0.5h | P0 |
| 五類 classifier 與 golden cases | HIGH | 1h | P0 |
| pair interval matching＋null | HIGH | 1–2h | P0 |
| cancer gene／drug source inventory 與 join QA | HIGH | 1–2h | P0 |
| single-cell／multi-region topology truth | HIGH | 外部依賴 | P1 |
| CN/purity/read-depth sensitivity | MED | 1–2h | P1 |

## §6 Evidence Conflict Scan

| Prior conclusion | Tier | Relation | Source |
|---|---|---|---|
| Region tree ≠ confirmed biological clone；需 single-cell／multi-region truth | L1/L2 | conflict with “proof” wording | `docs/CURRENT_FOCUS.md` |
| Cross-region correlation 曾受 shared-read-count confound | concluded NEGATIVE | dependent; 必須避免把空間相關當獨立驗證 | `docs/experiments/INDEX.md` O13 |
| Window aggregation AUC 可由 genomic autocorrelation 製造 | concluded NEGATIVE | dependent; 癌症基因 enrichment 必須用 matched null | `docs/experiments/INDEX.md` P3 |
| VAF ranking 受 CN/purity/read depth 混淆 | L2 | dependent; 只能稱 VAF-supported candidate | 2026-07-11 topology report |
| clean layered-v3 尚未在 7/11 cutoff 完成 | L1 | conflict with “latest final” | `docs/CURRENT_FOCUS.md` |

**Conflict count**：5；反例風險維度維持 0。

## §7 Decision Threshold + Path

- **TOTAL**：50/100。
- **Verdict**：PROBE。
- **Next action**：先完成 HCC1395 pair 的五類守恆、區域配對與 null；通過 checkpoint 後將觀察支撐升 20、反例風險升 10，重評 70/100 GO，再擴全 dataset 與 annotation。
- **Decision lock**：N；最後 claim ceiling 依 clean-v3 與正交 truth 狀態決定。

## Provenance Footer

- **Commit hash**：`6067568637088838a9f518955e41d222f057e4f1`
- **Build time**：2026-07-12 10:00:00 +08:00
- **Skill version**：`/pre-decision-audit` v0.1
- **Audit JSON**：`InterSubMod/state/cycles/20260712-1000-hcc1395-pair-coarse-topology-validation/audit.json`

## §8 Post-Probe Reassessment（2026-07-12 00:54 +08:00）

> **更新 verdict：GO（70/100），但只授權 engineering reproducibility 報告；biological validity／efficacy proof 仍 NO-GO。**

| Checkpoint evidence | 結果 | 判定 |
|---|---:|---|
| Exact matched regions | 6,252；complete-both 5,720 | PASS，遠高於 1,000 |
| 五類 raw agreement / κ | 69.39% / 0.4973 | 中度重現，不是可互換 |
| Chromosome-preserving null | 39.51%；κ 95% null upper=0.0231；p=1/5001 | PASS |
| Binary resolved/unresolved | agreement=73.15%；κ=0.4587 | PASS；主要差異是能否解出唯一 topology |
| 排除 dominant 同格影響 | 移除 Topo>1/Topo>1 同格後 agreement=45.50% | qualified PASS；仍高於 headline null mean，但未另建 restricted null |
| Exact/overlap sensitivity | agreement 69.22%–69.43%；κ 0.492–0.496 | PASS，matching rule 不改方向 |
| 分類與配對守恆 | 64/64 checks PASS；9 TSV 獨立重跑 byte-identical | PASS |
| Annotation matched null | CGC +3.15 pp，p=0.585；approved-antineoplastic +3.11 pp，p=0.430 | 不支持 annotation 額外驗證 topology |
| Clean-v3 release gate | producer 6/7；canonical/sensitivity roots 與 verification absent | FAIL；historical-v2 PARTIAL 維持 |

### 更新 Credibility Score

| 維度 | 初始 | 更新 | 理由 |
|---|---:|---:|---|
| 理論基礎 | 20 | 20 | classifier 與 forest semantics 不變 |
| 觀察支撐 | 10 | 20 | 5,720 exact complete pairs、五種 matching sensitivity、5,000 permutations |
| 機制清晰度 | 10 | 10 | coarse→shape-set→exact-tree 三層重現率揭露解析度差異 |
| 反例風險 | 0 | 10 | dominant-class、phase swap、annotation strata、shared upstream 均已量化；仍缺 biological truth |
| 所需資源 | 10 | 10 | 全量分析已在可重現腳本內完成 |
| **TOTAL** | **50** | **70 / 100** | **GO for qualified engineering report** |

**Decision lock**：本輪可產出 HTML 並回答 technical reproducibility；不得把 GO 延伸成 accuracy、biological clone truth、clinical actionability 或 method efficacy proof。clean-v3 canonical+sensitivity `_SUCCESS` 與正交 topology truth 到位後才重開 scientific gate。
