---
title: layered_workstation 癌基因/用藥 + 多層證據摺疊 build 操作紀錄與檢核報告
date: 2026-07-19
build_branch: research/subclonal-reconstruction-202606
scope: docs/methodology/_assets/layered_workstation（current-v5, 7 datasets / 6 biological samples, chr1-22）
status: built + verified (pilot HCC1395 → finalized 7/7)
scientific_release: NO-GO（annotation-stratified context；非拓撲真值/非臨床可用性）
---

# layered_workstation：癌基因/用藥 + 多層證據摺疊 — 操作紀錄與檢核報告

> **TL;DR** — 在 canonical layered_workstation 的每個 region 加入：①癌基因/用藥摺疊面板（COSMIC CGC v104 + DGIdb，放 Observed site evidence 之前）+ verdict 癌基因 badge；②Observed 證據多層摺疊（HP1/HP2 各自 + 唯一第一/候選空間 network/locus 突變排列 各自展開）；③勾選狀況 checkbox（只看癌基因/有藥）；④GRCh38 全基因組分布標癌基因（358 紅點）。全改在**產生 HTML 的 Python SoT**，7/7 樣本重生、canonical 覆蓋、playwright+VLM 逐項驗證、**0 pageerror**。

## 1. 需求（使用者）
1. 每個 region 的 **Observed site evidence 之前**加摺疊區放癌基因位置 + 用藥資訊。
2. Observed 底下整理成**多層摺疊**：HP1/HP2 各自展開、唯一第一順位、候選空間 network 各自展開。
3. **主要資訊 vs 輔助資訊**分清 + 文字顯示與大小依重點程度調整。
4. **勾選要觀察的狀況** checkbox。
5. 每個 region 的 **locus 突變排列**觀察表示。
6. **GRCh38 全基因組分布**清楚標記有交集的癌基因區域。
7. 整體美觀調整。
8. 確認後整理操作紀錄與檢核報告（本檔）。

## 2. 根因與關鍵修正（v2/v3 → current-v5 backbone 漂移）
「之前找的癌基因/用藥」資料在 `research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/hcc1395_pair_cancer_drug_region_hits.tsv`（COSMIC CGC v104 + DGIdb + GENCODE v46）。
- 🔴 **問題**：該 TSV 的 `region` key 來自**舊 layered-v2 backbone**，與 current-v5 workstation 的 region key **不匹配**（例：PTEN 舊 `chr10:87840023-87928739` vs v5 `chr10:87818272-87928739`）。exact-key join 幾乎全落空。
- ✅ **修正**：改用 **GENCODE 基因座標重疊 join**（`gene_start`/`gene_end` 為 backbone-independent），對 current-v5 region 做座標 overlap → 正確命中。**8,222 regions 中 358 命中 COSMIC CGC 癌基因**。

## 3. 操作紀錄（改的都是產生 HTML 的 Python SoT）
### 3.1 `docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py`（單樣本 builder）
| # | 改動 | 用途 |
|---|------|------|
| 1 | `load_gene_drug()` + `genes_overlapping()` + `_parse_region_span()` + `_safe_int()` | 讀 gene-drug TSV → per-chrom 基因區間；座標重疊 join |
| 2 | `main()`：`SM_GENE_DRUG_TSV`（optional overlay，不破壞既有嚴格 SHA 契約）→ 每 region 附 `gene_drug`、index 附 `has_cgc`/`cancer_genes`/`n_antineo_drugs`、`page_data.gene_drug_meta` | 資料注入（缺 env 則優雅空白） |
| 3 | JS `geneDrugPanel()` + `cancerGenesOf()` | 癌基因/用藥摺疊面板（放 Observed **之前**）+ claim ceiling |
| 4 | `renderDetail()`：verdict 加 `cgcBadge`（主要）；`out+=geneDrugPanel(region)` before Observed | 主要 badge + 面板注入 |
| 5 | `laneHtml()` 重構：`<section>`→`<details>`（HP1/HP2 各自摺疊）；內部 `唯一第一` / `候選空間 network` / `locus 突變排列` 各自 `<details class="sub-drawer">` | 多層摺疊 + 主要/輔助 |
| 6 | JS `locusArrangement()` | read-state × 位點 A/R/· 網格（locus 突變排列） |
| 7 | filter UI + `applyFilters()` 加 `fcgc`/`fdrug` checkbox + 事件 | 勾選狀況（只看癌基因/有抗腫瘤藥） |
| 8 | ideogram region mark 加 `.ideogram-cgc` 紅點 + 圖例 | GRCh38 全基因組分布標癌基因 |
| 9 | CSS：`.badge.cgc` / `.gene-drug` / `.sub-drawer` / `.locus-table` / `.filter-checks` / `.ideogram-cgc` + 文字階層 | 主要大/粗、輔助小/muted；紅=癌、藍=結構 |

### 3.2 `docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py`（7 樣本 driver）
- 加 `GENE_DRUG_TSV_BY_SAMPLE`（HCC1395 → cancer_drug_region_hits.tsv）+ `collect_rows()` 依樣本帶 `SM_GENE_DRUG_TSV`。→ 之後全量重建自動帶 gene-drug；其他 6 樣本無資料則優雅空白。

## 4. 檢核報告（playwright + VLM，§13.7 fresh 驗證）
| 項 | 驗證方法 | 結果 |
|----|---------|------|
| 癌基因面板在 Observed **之前** | DOM `geneBeforeObserved` | True ✓ |
| verdict 癌基因 badge（主要） | DOM `.badge.cgc` | 🧬 PTEN ✓ |
| PTEN CGC T1·TSG + 60 抗腫瘤藥 | VLM 讀截圖 | 正確 ✓ |
| HP1/HP2 各自摺疊 + 3 sub-drawer（唯一第一/network/locus） | DOM `details.lane` + `.sub-drawer` | 2 lanes + 3 sub ✓ |
| locus 突變排列網格 | DOM `.locus-table` + VLM | 渲染（XRXX→·R·· 75 reads）✓ |
| 勾選「只看癌基因命中」 | 勾選後 count | 8,222 → 358 ✓ |
| GRCh38 分布標癌基因 | DOM `.ideogram-cgc` + VLM | 358 紅點跨 22 染色體 ✓ |
| canonical HCC1395.html | playwright | gene=True, badge=🧬PTEN, subDrawers=3, locus=1, **errs=0** |
| canonical COLO829.html（無 gene 資料優雅） | playwright | gene=False, subDrawers=6, locus=2, **errs=0** |
| canonical index.html | playwright | 渲染, **errs=0** |

**Finalize**：driver 重生 **7/7 hash-bound pages + canonical index**（22:20-22:21 覆蓋），SUMMARY_SHA256=`71da78b6…`，CANONICAL_RUN=`20260713_layered_reconstruction_v3_raw_all_lps_pass_v5`。

## 5. 誠實界線（claim ceiling，全程維持）
- 癌基因命中＝「region 座標與 GENCODE 基因交集」；用藥＝「該基因 DGIdb 關聯 approved antineoplastic 藥物」。**非拓撲為真、非候選排序依據、非臨床可用性**（沿用 `gene_drug_inventory` 定案）。
- HCC1395 兩 dataset 為同一 biological sample 的技術處理（非獨立生物重現）。
- locus 突變排列＝read 直接觀測，**非 phasing 推斷、非 CCF**。
- 資料版本＝current-v5（`20260713_…_v5`）；scientific release 仍 **NO-GO**。

## 6. 待辦
- 其他 6 樣本目前無 CGC/DGIdb region-hit 資料（優雅空白）；若日後產出，加進 `GENE_DRUG_TSV_BY_SAMPLE` 即自動生效。
- 全案 469,849-sSNV 新 run render 進 workstation 後，gene-drug overlay 因採**基因座標** join，可直接沿用（無需改 join 邏輯）。
