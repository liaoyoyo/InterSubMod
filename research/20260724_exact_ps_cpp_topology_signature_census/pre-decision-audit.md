<!--
建立時間: 2026-07-25
目標: 在 exact-PS layered workstation 加入全 7 組癌症基因與 approved-antineoplastic 藥物同基因交集前，固定資料 grain、反例與 claim ceiling
處理範圍: 7 technical datasets / 6 biological samples / GRCh38 chr1-22 / 98,955 exact-PS final groups
關聯檔案:
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/gene_drug_inventory.md
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/implementation-notes.md
-->

# Exact-PS 癌症基因 × 藥物顯示 pre-decision audit

## §0 Cynefin domain

**Complicated / GO with explicit invariants**。GENCODE interval join、HGNC bridge、CGC 與 DGIdb boolean intersection 都有可重複的 deterministic 規則；不把資料庫交集升格成臨床結論。

## §1 Observation completeness

| 觀察 | 狀態 | Tier | 證據 |
|---|---|---|---|
| GENCODE v46、CGC v104、DGIdb 本地來源可讀且有 schema | ✓ | L2 | `gene_drug_inventory.md` §2–3 |
| 舊 HCC region TSV 不是 all7 authority | ✓ | L2 | 只含 HCC technical pair；新版 exact-boundary match 低 |
| actual locus join 可覆蓋 all7 98,955 groups | ✓ | L2 | 只讀全量稽核：CGC=3,554 groups、同基因 CGC∩drug=1,252 groups |
| HCC1395 舊 PTEN envelope 在新版不構成 PTEN locus hit | ✓ | L1 | exact loci 87,818,272 / 87,840,023；PTEN body 87,862,638–87,971,930 |
| gene-drug context 可驗證拓撲或療效 | ✗ | L4 | 既有 inventory scientific NO-GO；只准 annotation context |

## §2 Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | actual locus → GENCODE gene → HGNC → CGC/DGIdb 同基因鏈清楚 |
| 觀察支撐 | 20 | 已全 7 組、chr1-22 只讀重算 |
| 機制清晰度 | 20 | deterministic 1-based inclusive join；AND 規則可直接測試 |
| 反例風險 | 10 | DGIdb release/version 混合；需強制 provenance 與 claim ceiling |
| 資源 | 10 | producer、7 HTML 重建與 browser audit 約 1–6 小時 |
| **TOTAL** | **80 / 100** | **GO** |

**Falsifier observable**：若 producer 的 all7 CGC/同基因交集不等於 3,554/1,252、任何樣本不等於預先重算值、PTEN target 被標成命中、或 browser 的 AND filter 回傳 non-CGC / 無 approved-antineoplastic gene，則 fail closed，不發布 HTML。

## §3 Assumption map

| 重要性 | 已知 | 未知 |
|---|---|---|
| 高 | GRCh38、1-based inclusive、actual sSNV position、HGNC same-gene AND | DGIdb snapshot 的完整 release identity |
| 低 | CGC/GENCODE boundary 差異只作 cross-check | drug claim 對特定 mutation、癌別與劑量的臨床方向 |

高重要未知的處理：DGIdb 僅稱 `local snapshot`，保留 source/version；不得寫 actionability。

## §4 Quick pilot checkpoint

1. 由四個 raw sources建立 annotation sidecar → 驗證：schema/hash/721 autosomal CGC GENCODE genes 守恆。
2. HCC1395 exact loci join → 驗證：CGC=480、CGC∩drug=174；PTEN target 為 `NO_HIT_EVALUATED`。
3. HCC1395 Chromium 桌機/手機 → 驗證：AND filter、locus badge、detail table、0 console error、0 page overflow。

只有三項全過才重建 all7。

## §5 Gap diagnosis

| 缺口 | 影響 | 工時 | 優先級 |
|---|---|---:|---|
| DGIdb 完整 release receipt | 限制對外版本聲明 | 外部協調 | P1 |
| mutation-specific resistance/actionability source | 不能推療效 | 另案 | P2 |
| COSMIC redistribution terms receipt | 對外散布需另確認 | 外部協調 | P1 |

## §6 Conflict scan

- `gene_drug_inventory.md` 的 scientific NO-GO 與「把 CGC/DGIdb 當 topology truth / clinical validation」衝突；本實作明確不做此 claim。
- `InterSubMod/MEMORY.md` 在 repo 中不存在，無法完成該路徑 grep；改由 `docs/CURRENT_FOCUS.md` 與 validated NEGATIVE 文件掃描，未找到禁止 annotation-context 顯示的結論。
- `docs/reports/validated/2026/04/20260401_LOH_weekly_review/06_methylation_hypothesis_negative.md` 與本 UI annotation 無直接方法衝突。

## §7 Decision

**GO，decision lock=Y**。只允許：

`actual sSNV locus ∈ GENCODE gene body ∩ COSMIC CGC ∩ DGIdb(approved=true ∧ anti_neoplastic=true)`。

任何 span-only、舊 region-key、不同 gene 的 OR join、topology truth、mutation sensitivity 或臨床 actionability 解讀均不進本版。

## §8 Post-implementation outcome

- Producer all7_v2：721 CGC genes、279 drug-linked；224 non-HGNC DGIdb
  rows 中 0 筆通過 GENCODE symbol gate、0 keys overlap production CGC。
- Falsifiers：all7 `3,554/1,252` 與七組逐樣本計數全部符合；HCC1395
  PTEN target 為 `NO_HIT_EVALUATED`；browser AND regions 都含同基因
  approved-antineoplastic row。
- Python regression 34/34 PASS。
- Chromium final：16 page states、44 PNG、84 six-mode interaction
  records、0 console error、0 external request、0 post-detail overflow。
- Claude Code Opus 5：四項 P1 全部關閉，P0/P1=0，結論
  `AGREE_WITH_NOTES` 且「合理，可提交」。

因此 pre-decision `GO` 已由 post evidence 支持；claim ceiling 不變。
