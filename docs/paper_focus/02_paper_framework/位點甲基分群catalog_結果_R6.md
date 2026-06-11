<!--
建立時間: 2026-06-09
狀態: results (位點甲基分群標籤 catalog — 骨架結果 + 7-TAG 分類 + 觀察分佈圖)
報告類型: paper_focus_catalog_results
受眾: 廖子游 · PI · 論文 session（Results R6 素材）
task_type: B_validation（全基因組 + 全 6 樣本）
partial_flag: "subset — cis/copy 欄僅 HCC1395 (816 HP-axis loci) 有真值；mod_type/gene 多數 NA；single-pipeline ClairS-TO 封頂 ⭐3"
data_sources: research/tsg_promoter_asm_reviewer/genome_survey_v2/catalog/catalog_skeleton.tsv,research/tsg_promoter_asm_reviewer/genome_survey_v2/catalog/catalog_tag_counts.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/catalog/catalog_audit.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/catalog/discovery_venn.json
scripts: research/tsg_promoter_asm_reviewer/scripts/98_build_catalog.py,research/tsg_promoter_asm_reviewer/scripts/99_dbeta_venn.py,research/tsg_promoter_asm_reviewer/scripts/100_catalog_figures.py
verification: workflow wf_5644ed77-082 (6/6 catalog 數字獨立 adversarial 重驗 PASS)
-->
<!-- provenance-verified: 全數字直接由 catalog_skeleton.tsv / catalog_tag_counts.json / discovery_venn.json grep；scripts 98-100 可重跑；wf_5644ed77-082 獨立重驗。 -->

# 位點甲基分群標籤 catalog — 骨架結果（Results R6）

> **L0 一眼結論**：把 6 樣本 **332,705 個 ISM 掃描位點** 各貼一個 TAG，回答「哪些位點甲基分群明顯、是什麼性質」。結果是**誠實的 characterization**：**reliable 甲基分群有 12,868 個（TAG-C），但 cis-test 證它們是 germline-allelic 背景，不是 somatic-ASM；全基因組只有 1 個位點（chr17/TBC1D16, TAG-A）是乾淨 somatic cis-ASM**。BRCA2 是 copy/subclone-confounded（TAG-B）。**甲基分群真實但壓倒性地不是 somatic 判別器 → 強化「characterization 非 filter」主軸。**
>
> **L1 重點邏輯**：
> ① **核心欄 `clustering_reliability` 三態**（reliable 11,689 / latent 28,254 / none）—— reliable 用 ISM 原生 CramersV Cochran gate（=`Significant`）；latent 用 PERMANOVA 救回被 gate 歸零的真分群。
> ② **每位點收斂 1 個 TAG（7 類）** → P5 per-TAG 長條 = R6 主圖。
> ③ **A/B/F 三個 cis-based TAG 只在 HCC1395**（816 HP-axis loci 有 cis 真值）；C/D/E/G 全基因組可判。
> ④ **守 characterization≠filter**：TAG-D（FP-prone）+ TAG-E（latent minN<15 FP-leaning）**只標註不當 filter**（latent enrich≈base 47:1，無判別力）。

---

## L2 — 1. 七 TAG 分類結果（P5 = Results R6 主圖）

| TAG | 意義 | 數量 | 佔比 | 代表位點 |
|-----|------|-----:|-----:|---------|
| **A** | 乾淨 cis somatic-ASM | **1** | 0.0003% | **chr17:79991120 / TBC1D16**（唯一）|
| **B** | copy/subclone-confounded | 5 | 0.0015% | **BRCA2/ZAR1L**（d_within −0.023 ≪ d_copy −0.11）+ chr5/chr14/chr20×2 |
| **C** | germline-allelic 背景 | 12,868 | 3.87% | reliable 分群但 cis-test=germline 背景（B2）|
| **D** | high-Δβ FP-prone artifact | 5 | 0.0015% | regression-to-extreme（**只標註不當 filter**）|
| **E** | latent（PERMANOVA 救回）| 28,254 | 8.49% | CramersV Cochran-gated 但 PERMANOVA 顯著 |
| **F** | untestable | 54 | 0.016% | 53 pure-ALT clonal + chr18:11741161 CpG-creation caveat |
| **G** | 無訊號 | 291,518 | 87.62% | 背景 |
| | **合計** | **332,705** | 100% | 6 樣本 × TP(294,723)/FP(3,177)/FN(34,805) |

> 🔑 **論文敘述**：「332,705 個位點分 7 類；reliable 甲基分群存在於 ~3.9% 位點（TAG-C），但 cis-test 證其為 germline-allelic 背景而非 somatic-ASM；全基因組僅 **1 個乾淨 somatic cis-ASM exemplar（chr17/TBC1D16）**。BRCA2 經 within-subclone copy-partition 證為 copy/subclone-confounded。」
> （圖：`InterSubMod/docs/paper_focus/04_figures/P5_per_tag_counts_R6.png` + `catalog_overview.png`）

---

## L2 — 2. 核心欄 `clustering_reliability`（T-CODE-2 audit）

**定義（ISM 原生 + script 88 canonical）**：

| 態 | 定義 | 全基因組計數 |
|----|------|-----------|
| **reliable** | `Significant` = PassedGating & GlobalP≤0.05 & **CramersV≥0.1** & NumReads≥20（CramersV 經 Cochran 最小期望格≥5 gate）| 11,689（TP 11,655 / FP 34）|
| **latent** | NOT Significant & PassedGating & NumReads≥20 & **ClusterPermanovaValid&P≤0.05** & LabelHPPermanovaValid&P≤0.05 & min(HP1FamN,HP2FamN)≥5 | 28,254 |
| **none** | 皆非 | 292,762 |

**過嚴審查結論（ISM-2）**：CramersV Cochran gate **確實偏嚴**——把稀疏表（小 HP group，如 3 vs 51）的 CramersV 歸零，PERMANOVA（distance-based, robust）能救回真分群。但 **latent 納回集無判別力**：

| minN | 救回 TP | 救回 FP | rec TP:FP | vs base 47:1 |
|-----:|-------:|-------:|----------:|-----:|
| 0 | 25,980 | 310 | ~84 | enrich ~1（characterization-only）|
| 5 | 25,380 | 279 | ~91 | ~base |
| 15 | 21,033 | 82 | ~257 | minN 越高越偏 TP 但 minN<15 FP-leaning |

> **決策（不直接改碼）**：建議 ISM 把 `clustering_reliability` 從二元 0/1 改 **三態 reliable/latent/unreliable**（latent 標但不當判別）；**🔴 若真要改 C++（RegionProcessor.cpp:1592）→ 走 /cpp-change + 編譯 Hard Gate**。納回 latent **只提升 characterization 分辨，絕不可當 filter**（enrich≈base = 無 TP/FP 判別力，與甲基-filter NEGATIVE 一致）。

### §2.5 — clustering_reliability 非 coverage / CN-dosage 假象（confound 驗證 ⭐2026-06-10）

> 答 schema P6「reliability 是否被覆蓋驅動」：**否**。core 欄 `clustering_reliability` 的可信度由此 backstop。完整：`InterSubMod/docs/experiments/in_progress/2026/06/20260610_methylation_clusterability_vs_coverage_CN_characterization_01.md`（ledger 97/98）。

| 驗證 | 結果 |
|------|------|
| **coverage（read 數）**| 微弱 ρ=0.109（HCC1395 TP）；18 sample×class cells 全 |ρ|<0.25；decile = power floor(<50×)後 plateau；**多變量控制後 coverage OR=0.91 轉負** → 非 coverage driver |
| **CN（拷貝數）**| 強相關但反 dosage 方向：clustered-rate neutral 0.222 > gain 0.179 ≫ **loh 0.022**；**LOH logistic OR=0.07**；gain 固定覆蓋亦不增 → 非 dosage 假象 |
| **真 driver = n_CpG**| logistic OR=1.20（唯一正向 power）；cov+nCpG+CN 僅解釋 8%（pseudo-R²=0.083）|
| **LOH-抑制機制跨樣本**| **6/6 樣本 OR<1**（0.01–0.107，3 癌種）；機制 = read-clustering 需 ≥2 haplotype 等位甲基差，LOH 移除第二 allele |
| **FP 例外**| FP 不顯示 LOH 抑制（OR=0.98）→ 暗示 FP 分群=noise（⚠ loh n=70 低 power；characterization-only 不可 filter）|

> 🔑 **意義**：分群是**真等位甲基生物學（需兩個 haplotype），非覆蓋/拷貝數技術假象** —— 與 fast_cnv（dosage 非量級 driver）+ HD-4（NGroups=phasing 非甲基）一致。⚠ 區分：**NGroups（HPFineNGroups）與 coverage 強相關 ρ=0.51 但那是 HP-tag occupancy 取樣非甲基分群**。圖 `04_figures/clusterability_vs_cov_cn.png` + `clusterability_cov_cn_extended.png`。

---

## L2 — 3. Δβ 互補發現通道（T-CODE-3，ISM-3 Venn）

ISM gate 在 discrete clustering（CramersV/PERMANOVA）；Δβ 是 per-position allele 平均差。兩者正交。**union 發現通道**的 Venn（全基因組 TP/FP，base TP:FP = 92.8:1）：

| cell | TP | FP | TP:FP | vs base | 解讀 |
|------|---:|---:|------:|--------:|------|
| **dbeta_only**（大Δβ·無分群）| 1,633 | 137 | **11.9** | **0.13** | 🔴 **最 FP-enriched** = D5 regression-to-extreme |
| both（分群+大Δβ）| 3,379 | 85 | 39.8 | 0.43 | |
| ISM_only（分群·無大Δβ）| 33,656 | 228 | **147.6** | **1.59** | 🟢 最乾淨（BRCA2 即此類：clustering-only 低 Δβ）|
| neither | 256,055 | 2,727 | 93.9 | 1.01 | ≈ base |

> **設計結論**：**Δβ 不可 naive 加**（dbeta_only cell FP-enriched 7×）。正確 union = **ISM-clustered OR（大-Δβ 且有 clustering）**——排除 dbeta_only 的 FP cell。Δβ 是**互補通道**（抓 ISM 漏的高平均差位點），**不取代 CramersV，必 clustering-gated**。
> 🔴 **加 per-position `dbeta_max` 欄進 ISM 輸出 = C++/Python 改動** → 走 /cpp-change。本 catalog 已用既有 `AlleleDelta`/`HPMergedDelta` 算 `dbeta_max`（證明值可得，差在做成一級發現通道）。

---

## L2 — 4. cis_status 欄（T-VAL-1 部分；HCC1395 816 HP-axis loci）

| cis_status | 條件 | 數量 | TAG |
|-----------|------|-----:|-----|
| **cis** | 有 within-control 且 d_within > 1.5×d_copy 且 perm p≤0.05 | 1（chr17, d_within 0.142 ≫ d_copy 0.029, p=0.001）| A |
| **copy** | within-control 顯示 d_copy 主導或邊際 | 5（BRCA2 −0.11 主導；chr5 邊際 0.059≈0.055）| B |
| **cis-mechanical** | within-dominant 但 CpG-creation caveat | 1（chr18:11741161 CREATES CpG）| F |
| **untestable** | T3 cis-candidate 但無 within-control（pure-ALT clonal）| 53 | F |
| NA | 非 cis-tested（全基因組 / 非 HCC1395 HP-axis）| ~332,645 | — |

> **誠實口徑**（對齊 2026-06-07 capstone）：chr17/TBC1D16 是**唯一**乾淨 cis exemplar；BRCA2 是 copy/subclone（HP1-1 是 longphase-S somatic subclone tag 非 copy；**% split 不 robust**）；6 個最強 untestable 需 COLO829/新法。**大規模協議仍卡「最強位點測不了」（pure-ALT/CGI-desert）非算力**——須 anchor SEQC2 CN + COLO829。

---

## L2 — 5. T-OBS-2 觀察例 + chr17 第二樣本

- **per-TAG canonical 例**：見 `catalog_tag_counts.json` `examples`（每 TAG ≤5 例，HCC1395 優先）。read×CpG heatmap 重生用既有 ISM display_v2（需 figs/ 重生）。
- 🔴 **chr17 第二樣本 = 無法直接重現**：chr17:79991120 的 somatic SNV **不在 HCC1937/HCC1954 的 variant set**（HCC1395 private somatic）→ 同位點 cis 無法在第二乳腺樣本直接驗。**升信心唯一路徑 = COLO829（不同癌種 private somatic）或 matched-normal diffASE**。**這正是 single-locus single-sample nominal p 的天花板。**

---

## L1 — 缺漏 / 前置依賴（誠實未補）

1. `mod_type`（5mC/5hmC/mixed）：依賴 MSA Level1 dup-bug 修復，骨架階段 = NA。
2. `gene`/`is_TSG_promoter`：無 genome-wide gene model bundled → 僅 curated（chr17/TBC1D16, BRCA2/ZAR1L）有值，其餘 NA。
3. `clustering_origin` ARI 欄：per-locus blind-ARI 僅 high-ARI 子集有值（q2 map）；其餘用 reliability+axis 推 germline-allelic 預設（B2）。
4. `cis_status` 全基因組未鋪：只 HCC1395 816 HP-axis loci 有真值；anchor SEQC2 CN 待補（T-VAL-1 完整協議）。
5. **single-pipeline**：所有 somatic_status 來自 ClairS-TO → tier 封頂 **⭐3**。

---

## L3 — Provenance + 可重跑

- **建構**：`scripts/98_build_catalog.py`（catalog + TAG + audit）/ `99_dbeta_venn.py`（Venn）/ `100_catalog_figures.py`（P1-P6）。
- **資料源**：`ism_existence_scan/<sample>_<cls>/significance_summary.csv`（117 欄/位點）+ `cis_scan_full.json`（816）+ `copy_partition_confirm.json` + `survivor_permutation.json`（chr17 p=0.001）。
- **獨立驗證**：workflow `wf_5644ed77-082` 6/6 catalog 數字 adversarial 重驗 PASS（332,705 列 / A=chr17 唯一 / BRCA2=B / sum / Venn FP-enrichment / reliable-latent）。
- **修正紀錄**：(a) `(p_cis or 1)` falsy-zero bug 修（chr17 p_cis=0.0 曾被誤判 copy）；(b) cis margin 1.5× + chr18 CpG-creation caveat → A=chr17 唯一（對齊 capstone）；(c) CRLF→LF 行尾。
