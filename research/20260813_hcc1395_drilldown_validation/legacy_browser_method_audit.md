<!--
建立時間: 2026-08-13
目標: 稽核 2026-07-26 legacy methyl allele-axis browser 的資料來源、分母、定義、校正、缺失、選案偏差與可重現性，並與 HCC1395_v1/current dashboard source 做 schema crosswalk
處理範圍: HCC1395、chr1-22、legacy 30,077 TP loci、current topology 19,849 sSNV；read-only，不修改 legacy/current dashboard 或外部 bundle
關聯檔案: InterSubMod/docs/reports/in_progress/2026/07/20260726_methyl_allele_axis_backbone_coenrichment/20260726_ssnv_branch_x_methyl_browser.standalone.html; /bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/index.html; InterSubMod/scripts/drilldown/
任務類型: B — Comprehensive validation（指定 HCC1395 與 22 條常染色體全量；非多樣本生物學 validation）
研究目標: G3 / G4 / G5
狀態: COMPLETE_READ_ONLY_AUDIT
-->

# Legacy methyl browser 方法與資料稽核

敘述框架：Claim–Evidence–Caveat + claim-to-source mapping（技術稽核）。

**TL;DR：legacy 數字可由落檔資料重算，但科學語意與重現鏈未達 validation 等級；不得把舊 A/B 類別映射成 current axis，也不得沿用「HP 陰性對照」或「density 已證實否決」的推論（影響：高，信心：高）。**

> Claim ceiling：`OBSERVATION_ONLY_NO_TRUTH_SET`。本稽核確認的是資料與方法 contract，不是生物學真值，也不是多樣本重現。

## 1. 結論與優先級

### P0 — 在任何科學比較或 current dashboard 重用前必修

1. **legacy A/B 不是文件所暗示的互斥 allele-axis vs HP-axis 比較。** 實際分類有未揭露的 B-first priority：`unstable → E`、`not aligned → D`、`V_hp≥0.7 → B`、`V_allele≥0.7 → A`、else C。全 30,077 列重建 mismatch = 0；但 4,025 個 B 中有 **3,544（88.0%）也同時滿足 `V_allele≥0.7`**，只有 481 是 HP-only。B 因此不能稱為純 HP「陰性對照」。
2. **headline 20,535 linkage hits 不是 `COREAD_MIN=6` 的 powered/confirmed linkage。** TSV 的 `n_linkage_partners` 與所有已記錄 pair（source code 只要求 `coread≥2`）逐列完全一致；powered hit 是 17,749，powered 且兩端 somatic-confirmed 是 17,141。A/B headline 由 `90.6% vs 68.3%` 改用 powered+both-somatic 後為 `82.8% vs 55.8%`。方向未反轉，但 estimand、分母與語意都改變。
3. **region graph 混合 TP 與 FP，且 branch gate 未要求 powered。** 7,415 regions 是在 TP∪FP raw-pair graph 中「至少含一個 legacy TP annotation」的 connected components；其中 1,341 regions 含 FP。472 個 branch+methyl regions 中 68 個含 FP、18 個沒有 powered branch edge、27 個沒有 powered+both-somatic branch edge。不可直接稱為「TP somatic backbone confirmed regions」。
4. **legacy 產物不是端到端可重現。** `build_page.py` 與 `build_browser.py` 只把 `/tmp/.../scratchpad/*.json` render 成 HTML/JSON；產生 `methyl_backbone.json` 與 `branch_methyl_cases.json` 的 analytic script 未在目標目錄或 repo 搜得。`n_branch_methyl_regions=472` 又是 build script hard-code。scratchpad 目前仍存在且可雜湊，但 `/tmp` 不是可持久 SoT。
5. **current HCC1395_v1 methyl provenance 也未封口。** receipt 對 topology/MLHP 有檔案 SHA256，但 ISM 只記 directory `/bip7_disk/liaoyoyo2001/ism_lineage_v1`，`sha256=""`；沒有 build command、git revision 或逐染色體 summary/run_params hash。built bundle 也尚未含 working-tree source 已提供的 `observationScope`、axis definitions 與 multiplicity metadata。

### P1 — 下一次 rebuild 應完成

1. 將 funnel 明確落檔：`30,077 TP loci → 7,415 annotated components → 518 with A → 472 branch+A+n=2..12 → 18 with A≥3 → 14 displayed`，每一步都附排除原因。
2. 將所有 region 統計分拆成 `raw pair`、`powered`、`powered+both-somatic`，另列 component 是否含 FP/unannotated TP。
3. 把局部候選數視為 **linkage opportunity/risk set**，用 pair-level powered rate，或 site-level count model + candidate-pair offset；不要再以「有任一 partner」的飽和率直接下 confounder/root-cause verdict。
4. 對 current v1 重新 build，讓 receipt 封存 ISM 22 份 `significance_summary.csv`、22 份 `run_params.json`、build argv、source revision、output hashes；並消除 receipt 中 MLHP `1,873 broken` 與 strict endpoint `0 / 11,590 broken` 的並存敘述。
5. 把 legacy/current coordinate crosswalk 另存成 versioned table，至少補 REF/ALT；目前只能做 `chrom:pos` crosswalk，不能證明 multiallelic site 的 allele identity。

### P2 — 可改善瀏覽與多樣本擴充

1. UI 同時顯示 sample、cohort、window、edge tier、axis contract 與 claim ceiling；legacy badge 明示 `w5000 / Cramér V / Tier-R raw pair`。
2. 對 A/HP 採 multi-label 呈現（A-only、HP-only、both、neither），移除 B-first 類別造成的語意遮蔽。
3. 熱圖選案只用於案例瀏覽；population metric 必須用完整母體。匯出目前 filter 結果與 exclusion ledger。
4. 多樣本擴充時先逐樣本輸出相同 schema、missingness 與 opportunity-adjusted effect，再做 random-effects/meta-analysis；不可直接 pool 不同 window/cohort/axis contract。

## 2. 稽核範圍、假設與驗證方式

### 2.1 關鍵假設

- legacy 與 current 都是 HCC1395/hg38 chr1-22；sample identity 有路徑與 receipt 支持。
- legacy `phylo.pos` 是 window start；`+5000` 後才是 SNV position。raw join = 1，校正後全 34,736 records 都可對上 linkage census；legacy report 自身也記錄交集由 1 升到 24,220 的當時診斷片段。
- current L1 只內嵌 `chrom:pos`，legacy TSV 也不含 REF/ALT，因此 crosswalk 僅為座標交集，**不是 allele-level identity**。
- working-tree dashboard source 正在與 bundle 分離演進；本文把 `/bip7_disk/.../HCC1395_v1` 視為已建 artifact，把 `InterSubMod/scripts/drilldown/` 視為 audit 時工作區 source，不假設二者由同一 revision build。

### 2.2 Step → Verify

1. 追溯 legacy raw source → annotation → locked JSON/HTML。  
   → 驗證：source hashes、row counts、join offset、class rule與 headline counts 可獨立重算。
2. 重建 raw-pair graph 與 browser funnel。  
   → 驗證：`regions.total=7,415`、`with_A=518`、`with_A2=104`、`A_and_B=239`、eligible `472` 均精確命中。
3. 解析 current bootData/L4/receipt。  
   → 驗證：current L1 `n=19,849`、ISM join `16,195/19,849`、window `±1,000 bp`、axis bit contract 與 source code一致。
4. legacy/current 全量座標 crosswalk。  
   → 驗證：兩側 unique count 無 duplicate；列出 overlap/only/Jaccard 與 class×current-axis matrix。
5. 對每項方法判定 reusable/non-reusable。  
   → 驗證：每項都有 evidence path/embedded key，以及不可推論理由。

## 3. Input 與 provenance

### 3.1 Legacy stable inputs

| 角色 | 輸入 | audit SHA256 / embedded key |
|---|---|---|
| standalone browser | `InterSubMod/docs/reports/in_progress/2026/07/20260726_methyl_allele_axis_backbone_coenrichment/20260726_ssnv_branch_x_methyl_browser.standalone.html` | `ecdc1890...c78ed` |
| browser locked summary | `InterSubMod/docs/reports/in_progress/2026/07/20260726_methyl_allele_axis_backbone_coenrichment/20260726_ssnv_branch_x_methyl_browser.data.json` | `0a1fde91...93ac7`; keys `n_cases=14`, `n_branch_methyl_regions=472` |
| main locked data | `InterSubMod/docs/reports/in_progress/2026/07/20260726_methyl_allele_axis_backbone_coenrichment/20260726_methyl_allele_axis_backbone_coenrichment.data.json` | `6cbe9fed...5185947`; `n_tp=30077`, `verdict=REFUTED_BY_DENSITY` |
| locus annotation | `InterSubMod/docs/reports/in_progress/2026/07/20260726_methyl_allele_axis_backbone_coenrichment/methyl_class_x_linkage_annotation.tsv` | `523b2674...e099d`; 30,077 data rows / 30,077 unique `chrom:snv_pos` |
| phylo/methyl source | `InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/phylo_cpp_wg_full_records_v6.json` | `3838e820...1db6`; 34,736 records = TP 30,077 + FP 4,659 |
| linkage source | `InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/sm_linkage_genomewide.json` | `f71560a8...45eea`; union 35,332, pairs 53,094 |
| renderer | `InterSubMod/docs/reports/in_progress/2026/07/20260726_methyl_allele_axis_backbone_coenrichment/build_browser.py` | `a0d3a0a4...0e98`; reads `/tmp` and hard-codes 472 |

### 3.2 Ephemeral legacy intermediates

| 輸入 | audit 狀態 |
|---|---|
| `/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/efb6e3d8-c4af-43d8-ac97-9dffdbec60ed/scratchpad/methyl_backbone.json` | audit 時存在；SHA256 `92e9ccbd...aa33`；與 locked main JSON（排除 renderer 補入 `sources`）一致，但不可當長期 SoT |
| `/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/efb6e3d8-c4af-43d8-ac97-9dffdbec60ed/scratchpad/branch_methyl_cases.json` | audit 時存在；SHA256 `78b0ec76...21f`；14 cases 含 locus_detail/image path，但生成程式未落 repo |

### 3.3 Current artifact

| 角色 | 輸入 | audit evidence |
|---|---|---|
| dashboard | `/bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/index.html` | SHA256 `29a2f7da...9ba8`; bootData `sample=HCC1395`, `l1.n=19849` |
| receipt | `/bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/receipt.json` | SHA256 `697c742e...a0b2c`; `built_at=2026-08-08 01:13` |
| selfcheck | `/bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/SELFCHECK.md` | SHA256 `aac8fffb...357f`; 10 PASS / 1 FAIL / 1 SKIP |
| current source | `InterSubMod/scripts/build_drilldown_dashboard.py`, `InterSubMod/scripts/drilldown/sources/ism.py`, `InterSubMod/scripts/drilldown/emit/payload.py` | working-tree contract；artifact 不含新版 `observationScope`，不可反推 exact build revision |

## 4. Legacy cohort、axis 與 linkage 的實際定義

### 4.1 Cohort 與 join

- phylo v6 universe：34,736 = TP 30,077 + FP 4,659。
- linkage census universe：35,332 = TP 30,490 + FP 4,842。
- legacy TSV：phylo 的全部 30,077 TP，且都是 linkage TP；linkage 另有 413 TP 不在 phylo/TSV。
- coordinate transform：`snv_pos = phylo.pos + 5000`。`phylo.pos` 是 ±5 kb region 的 start，不是 mutation coordinate。
- TSV 30,077 rows 在 `chrom:snv_pos` 層級全 unique；15 欄中只有無 pair 的 9,542 rows 會合理地讓 `top_rel`/`partners` 空白。

### 4.2 Legacy methyl classes：文件定義不足，實際是 priority taxonomy

原始 `V_hp` / `V_allele` 是 coarse methyl cluster label 對 HP / ALT-support label 的 Cramér's V。來源程式先把兩值 round 到兩位小數，並以 `max(V_hp,V_allele)≥0.3` 設 `aligned`；legacy report 再於 locked records 上用 0.7 分類。

| 實際順序 | 類別 | n |
|---:|---|---:|
| 1 | `unstable → E_unstable` | 2,085 |
| 2 | `not aligned → D_split_unaligned` | 22,454 |
| 3 | stable+aligned 且 `V_hp≥0.7 → B_HP_ASM` | 4,025 |
| 4 | stable+aligned 且 `V_allele≥0.7 → A_ALLELE_clean` | 748 |
| 5 | 其餘 `C_other` | 765 |

獨立重建結果：30,077 / 30,077 完全一致，mismatch = 0。可是兩個門檻本身不是互斥：

| stable+aligned criterion | n | 最終類別 |
|---|---:|---|
| neither | 25,304 | C/D/E，依上游狀態 |
| allele-only | 748 | A |
| HP-only | 481 | B |
| allele+HP both | 3,544 | **全被 B-first 收進 B** |

此外有 98 loci 至少一軸剛好等於 round 後的 `0.70`，顯示 threshold 對 upstream rounding 有邊界敏感性。0.7 的選擇也未見 pre-specified rationale 或 sensitivity sweep。

### 4.3 Legacy linkage outcome：raw observed pair，不是 confirmed pair

`InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/scripts/sm_linkage_genomewide.py` 的 contract：

- 候選 pair：同染色體、距離 ≤50 kb。
- 只有 `coread≥2` 才寫入 `pairs`。
- `powered = coread≥6`。
- `both_somatic` 另由 matched-normal `normal_ref≥4` 且 normal VAF<0.05 判斷。

TSV `n_linkage_partners` 與「所有已寫入 raw pair」逐列 mismatch = 0；與 powered pair 有 7,799 rows mismatch。因此 headline estimand 是「至少一個 `coread≥2` raw pair」：

| Population | n | raw pair hit | powered hit | powered + both-somatic hit |
|---|---:|---:|---:|---:|
| all legacy TP | 30,077 | 20,535 (68.3%) | 17,749 (59.0%) | 17,141 (57.0%) |
| A | 748 | 678 (90.6%) | 642 (85.8%) | 619 (82.8%) |
| B | 4,025 | 2,750 (68.3%) | 2,346 (58.3%) | 2,247 (55.8%) |

這不是說 powered re-analysis 支持生物學關聯；它只證明原頁引用 `COREAD_MIN=6` 時，headline numerator 並未套用該門檻。

### 4.4 Region graph 與 472 的獨立重建

以 linkage union 35,332 nodes 與全部 53,094 raw pairs 建 undirected graph，再取「至少含一個 legacy annotation」的 connected components，可精確得到：

- regions total = 7,415；with A = 518；with ≥2 A = 104；A and B = 239。
- selection `A≥1 AND any(mutual_excl/nested) AND 2≤full-component n≤12` = **472**。
- 472 中，有任一 powered branch = 454；有任一 powered+both-somatic branch = 445。
- 7,415 中 1,399 含未 annotation nodes；1,341 含 FP、68 含 linkage-only TP。472 中相對為 71、68、3。
- 14 displayed cases 中 4 個包含 FP node：`chr1:115939392-116066164`、`chr2:21797623-21831456`、`chr2:18068480-18123326`、`chr7:51943525-52087777`。

這也解釋某些卡片 `n_ssnv > n_annotated`；SVG/table 只畫 legacy-annotated loci，但 header 的 region size 使用 full TP∪FP component。

## 5. 統計校正、missingness 與 selection bias

### 5.1 Legacy correction / inference ceiling

- legacy A/B taxonomy 使用 Cramér's V effect-size threshold，不是 p-value gate；此層沒有 multiple-testing correction。
- v6 records 雖另帶 `ls_hp_p` / `ls_al_p` 等 PERMANOVA 欄位，legacy class rule 沒使用它們；不可拿上游別的 FDR 流程替本頁補正。
- headline、CN×coverage strata、density strata 都沒有 CI、formal adjusted model、預先規格化 contrast 或 out-of-sample replicate。
- density strata 顯示的是 density 1、2、3–4、5–8、≥9；**density=0 的 A 57、B 1,095 被排除**。因此顯示 strata 的 denominator 只有 A 691/748、B 2,930/4,025。
- outcome 是「是否有任何 ≤50 kb pair」，而 density 又是「±50 kb 內其他 TP sSNV 數」。density 是 outcome opportunity 的直接構成量；高 density 下 any-hit 飽和至 100% 是機械性預期。這是應改 estimand 的證據，不足以單獨證明 causal confounder/root cause。

**評估：**將原本 co-enrichment claim 降為 descriptive 是正確且保守的；但 `REFUTED_BY_DENSITY` /「根因」仍比資料能支持的強。較準確用語是：`UNRESOLVED_AFTER_OPPORTUNITY_IMBALANCE; DO_NOT_INFER_ASSOCIATION`。

### 5.2 Legacy missingness

- linkage TP 30,490 中，30,077（98.65%）有 phylo/legacy annotation，413（1.35%）沒有；排除機制未在 browser 鎖定。
- component 定義又容許 FP/未 annotation TP 當 bridge；所以 site-level annotation complete 不等於 region-level composition complete。
- `top_rel` / `partners` 的 9,542 空值代表 raw-pair absence，而非 methyl missing；UI 應明示結構性空值。

### 5.3 Browser display selection

- 14 / 472 = 2.97% eligible regions 被展示；14 / 7,415 = 0.19% annotated components。
- full eligible set 中 A count 分布：A=1 378、A=2 76、A=3 11、A=4 5、A=5 1、A=6 1；因此有 **18 個 A≥3**，但頁面只含 14。
- audit 時 14 displayed 的全部 52 A-locus 均有兩張 heatmap；4 個未展示但同樣 A≥3 的 regions，其 12 個 A-locus 也都有 `distance_heatmap.png` 與 `cluster_heatmap.png`。這 4 個是：`chr10:132767774-132944853`、`chr19:51782649-51863021`、`chr1:111931189-111982877`、`chr1:57483998-57521711`。
- repo 內沒有 selection script、固定 top-N 或同分 tie-break contract，可解釋 18→14。因此 14 cases 應標 `curated display subset, selection rule incomplete`，不能作比例或代表性判讀。

## 6. Legacy 30,077 vs current v1 19,849 crosswalk

### 6.1 Coordinate overlap

| Metric | n / rate |
|---|---:|
| legacy unique `chrom:pos` | 30,077 |
| current topology unique `chrom:pos` | 19,849 |
| overlap | **16,566** |
| legacy-only | **13,511** |
| current-only | **3,283** |
| overlap / legacy | 55.079% |
| overlap / current | 83.460% |
| coordinate Jaccard | 49.658% |
| duplicates | legacy 0；current 0 |

注意：兩側 surface 都無 REF/ALT，以上只能稱 coordinate overlap。

### 6.2 Schema 明確不等價

| 維度 | Legacy browser | Current HCC1395_v1 | 可否直接 mapping |
|---|---|---|---|
| 母體 | Jan-2026 pileup TP 30,077 | exact-PS topology active-position union 19,849；11,590 regions | 否 |
| methyl window | ±5,000 bp（span 10,001） | ±1,000 bp（span 2,001） | 否；CpG/read composition 不同 |
| axis statistic | methyl cluster × HP/ALT label 的 Cramér's V | HP merged / ALT-REF / lineage 的 p/effect fields | 否 |
| category shape | exclusive 5-class、B-first | multi-label bitmask；同時 HP+ALT 合法 | 否 |
| significance | V≥0.7 heuristic；無 p/FDR | raw p≤0.05 exploratory；無 cohort-wide correction | 否 |
| linkage universe | TP∪FP、≤50 kb raw recorded pair component | exact-PS topology + strict endpoint edges | 否 |
| primary edge gate | record `coread≥2`；powered≥6 另存 | primary strict endpoint threshold（receipt 76,202/117,760 passed） | 否 |
| inference | HCC1395 single-sample descriptive | observation-only/no truth set | 都不可稱 validation |

### 6.3 Overlap 內 legacy class × current axis bitmask

current code：0=no current L4 summary；1=HP；2=ALT/REF；4=lineage；bitwise combination；8=至少一軸有測但都不顯著。此表是 cross-tab，**不是轉換表**。

| Legacy class | overlap n | 0 no L4 | 1 HP | 2 ALT | 3 HP+ALT | 4 lin | 5 HP+lin | 6 ALT+lin | 7 all | 8 none sig |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A_ALLELE_clean | 569 | 14 | 22 | 223 | 123 | 0 | 13 | 6 | 60 | 108 |
| B_HP_ASM | 2,197 | 48 | 268 | 443 | 630 | 7 | 144 | 20 | 447 | 190 |
| C_other | 453 | 23 | 73 | 108 | 61 | 9 | 25 | 0 | 23 | 131 |
| D_split_unaligned | 12,104 | 375 | 1,290 | 3,500 | 983 | 459 | 465 | 166 | 520 | 4,346 |
| E_unstable | 1,243 | 34 | 137 | 319 | 275 | 17 | 70 | 7 | 179 | 205 |

在 overlap 中，legacy A 共有 569 loci，其中 14 是 current axis code 0（無 L4 summary），所以 current L4-available 分母是 555；其 ALT raw-significant 為 412/555（74.2%）、HP 為 218/555（39.3%）。legacy B 以相同 L4-available 口徑有 2,149 loci，其中 current HP raw-significant 為 1,489/2,149（69.3%），但 current ALT 也有 1,540/2,149（71.7%）。這直接反證「legacy A=現在 ALT、legacy B=現在 HP」的一對一 mapping。Crosswalk JSON 的 set Jaccard 使用完整 shared-coordinate set，分母與這裡的 L4-available 百分比不同。

另外 overlap 的 current L4 availability 是 16,072/16,566（97.0%），遠高於 current 全母體 16,195/19,849（81.59%）；只分析 overlap 會產生 availability enrichment，不能代表 current universe。

## 7. Current artifact/source 的方法觀察

### 7.1 Current source 中可保留的 contract

`InterSubMod/scripts/drilldown/emit/payload.py` 的 audit-time working-tree source 已明定：

- global axes 只含 HP merged、ALT/REF、lineage；HP fine 只在 detail；cluster 軸因 double-dipping 不進結論。
- `raw p≤0.05` 只作 exploratory flag，未做 cohort-wide multiplicity correction。
- `AXIS_UNTESTED=16` 與 no-summary `0` 分開。
- `observation_scope()` 同時輸出 universe、summary join、testable、raw significant、detail available、missingness 與 legacy non-crosswalk warning。

這些是可重用的修正方向，但 `/bip7_disk/.../HCC1395_v1/index.html` 的 bootData 目前沒有 `observationScope`，axis legend 仍是舊版，故不能說 current artifact 已實現。

### 7.2 Current artifact 仍需修的證據鏈

- receipt：topology SHA256 `dc568237...28f6`、MLHP `7dc4084a...512`；ISM directory SHA256 為空，summary 29,572 rows 標「無 receipt 可比對」。
- current L1：19,849；summary join 16,195；detail dir 16,304；差異 109 也應在 funnel 說明。
- current selfcheck C10：cluster 20,903/20,904≈100%，已正確判 double-dipping，不可作證據。
- current receipt 同時保留 MLHP pattern 法的 `1,873/11,590 broken` 與 strict endpoint edges 的 `0/11,590 broken`。後者才是 component connectivity 的權威來源；前者應標 withdrawn/stale diagnostic，而非兩個都作有效 probe。

## 8. 可重用與不可重用方法

### 8.1 可重用（需照 caveat 重建）

| 方法 | Evidence path / key | 可重用範圍 | Caveat |
|---|---|---|---|
| universe→eligible→displayed funnel | legacy browser footer；current `InterSubMod/scripts/drilldown/emit/payload.py::observation_scope` | 所有 sample/dashboard | 每層要 machine-readable，不能 hard-code |
| source-locked JSON + required-key refuse | legacy `build_page.py` lines 18-29、`build_browser.py` lines 22-42 | renderer 防手打數字 | 只能保證 render input 完整，不保證分析正確或 input 可重生 |
| offset calibration | main legacy HTML `join offset +5000`；TSV `window_start`/`snv_pos` | 舊 cohort provenance | 僅限 w5000 legacy；current 不可套用 |
| confounder/opportunity sentinel | legacy density strata；current ISM×k sentinel | 提醒 selection/opportunity imbalance | 應換成 formal risk-set/offset analysis，不能直接下 causal verdict |
| per-locus heatmap + structure table | 14 browser case cards、base64 assets | candidate diagnosis / screenshot QA | 不可由 heatmap-selected subset估 population rate |
| strict endpoint edge receipt | current `InterSubMod/scripts/drilldown/sources/strict_edges.py` | exact-PS component/connectivity QA | 不得用 MLHP marginalized patterns替代 |
| observation-only claim ceiling | legacy L3 restrictions；current working-tree `claim_ceiling` | 多樣本一致 UI | claim ceiling 必須在 artifact 內，不只在 source |

### 8.2 不可直接重用

| 項目 | Evidence | 不可推論理由 |
|---|---|---|
| legacy A/B 類別作 current axis | 16,566 cross-tab | 母體、window、statistic、exclusive/multilabel contract都不同 |
| B 作 HP negative control | 3,544/4,025 B 同時 allele≥0.7 | control 被 exposure criterion 大量污染，且 priority 未揭露 |
| `90.6 vs 68.3` 作 confirmed linkage | raw-pair逐列 match；powered只17,749 | headline 沒套 `COREAD_MIN=6`/both-somatic |
| `REFUTED_BY_DENSITY` 作正式統計結論 | density=0 omitted；無 CI/model；outcome直接依賴 opportunity | 只能支持降級與重設 estimand，不能確認 root cause |
| 472/14 作 prevalence | hard-coded 472、14 cases來自 ephemeral JSON、18→14 tie-break缺失 | selection mechanism不完整且受圖片/人工選案影響 |
| legacy raw graph 作 exact-PS topology | TP∪FP、raw `coread≥2`、≤50kb | component semantics 與 current strict exact-PS 不同 |
| `/tmp` scratch / archive image path作 SoT | renderer lines 17-23；scratchpad path | 非持久、無 receipt、生成腳本缺失 |
| 單樣本 HCC1395 effect 外推多樣本 | legacy/readme明示 single sample | 無 between-sample variance、無 harmonized denominator |

## 9. Claim-to-source map

| Claim | Source / embedded key | 重算結果 | Verdict |
|---|---|---:|---|
| legacy TP n=30,077 | phylo v6 `set=TP`; locked `n_tp` | 30,077 | verified |
| class counts A/B/C/D/E | TSV `methyl_class` | 748/4,025/765/22,454/2,085 | verified |
| A/B class definitions | main HTML line 236 + phylo fields | displayed formulas overlap；actual B-first rule 0 mismatch | definition incomplete |
| all/A/B hit | TSV `n_linkage_partners>0` | 20,535/678/2,750 | arithmetic verified; raw-pair semantics |
| COREAD_MIN=6 | linkage `params.COREAD_MIN` | 6 | source verified; headline did not apply it |
| region metrics | locked `regions` | 7,415/518/104/239 | independently verified from union graph |
| 472 branch+methyl | browser JSON / build hard-code | graph rule reproduces 472 | value verified; provenance not reproducible |
| 14 display cases | browser JSON `n_cases` | 14 | verified; 18→14 rule unresolved |
| density strata | locked `density_strata` | shown bins sum A691/B2930 | arithmetic verified; zero-density omitted |
| density verdict | locked `verdict` | `REFUTED_BY_DENSITY` | over-strong inference |
| current n=19,849 | bootData `l1.n`; receipt topology counts | 19,849 | verified |
| current ISM join | receipt `ism_dirs.linkage` | 16,195/19,849=81.59% | verified |
| current window | receipt `window_size_bp` | ±1,000 bp | verified |
| cluster invalid | selfcheck C10 | 20,903/20,904≈100% | verified double-dipping sentinel |
| legacy/current overlap | TSV vs bootData L1 | 16,566 / 13,511 / 3,283 | verified coordinate-only |

## 10. 執行命令、I/O 與實際輸出片段

所有命令皆由 `/big7_disk/liaoyoyo2001/InterSubMod` 執行，皆為 read-only；未執行 legacy build scripts，避免覆寫產物或觸發其內含重跑/清理邏輯。

### 10.1 Hash 與 stable input identity

輸入：上述 legacy files、raw JSON、current `index.html/receipt.json/SELFCHECK.md`。  
命令：

```bash
sha256sum \
  /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/in_progress/2026/07/20260726_methyl_allele_axis_backbone_coenrichment/20260726_ssnv_branch_x_methyl_browser.standalone.html \
  /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/in_progress/2026/07/20260726_methyl_allele_axis_backbone_coenrichment/methyl_class_x_linkage_annotation.tsv \
  /big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/phylo_cpp_wg_full_records_v6.json \
  /big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/sm_linkage_genomewide.json \
  /bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/index.html \
  /bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/receipt.json
```

輸出：stdout only。實際片段：

```text
ecdc1890...c78ed  .../20260726_ssnv_branch_x_methyl_browser.standalone.html
523b2674...e099d  .../methyl_class_x_linkage_annotation.tsv
3838e820...1db6  .../phylo_cpp_wg_full_records_v6.json
f71560a8...45eea  .../sm_linkage_genomewide.json
29a2f7da...9ba8  /bip7_disk/.../HCC1395_v1/index.html
697c742e...a0b2c /bip7_disk/.../HCC1395_v1/receipt.json
```

### 10.2 Legacy class/outcome audit

輸入：legacy TSV、phylo v6、linkage JSON。  
命令核心：

```bash
python - <<'PY'
import csv, json
from collections import Counter
tsv='/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/in_progress/2026/07/20260726_methyl_allele_axis_backbone_coenrichment/methyl_class_x_linkage_annotation.tsv'
phy='/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/phylo_cpp_wg_full_records_v6.json'
lnk='/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/sm_linkage_genomewide.json'
with open(tsv,newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
P={(r['chrom'],int(r['pos'])):r for r in json.load(open(phy))}
def cls(q):
    if q['unstable']: return 'E_unstable'
    if not q['aligned']: return 'D_split_unaligned'
    if float(q['V_hp'])>=.7: return 'B_HP_ASM'
    if float(q['V_allele'])>=.7: return 'A_ALLELE_clean'
    return 'C_other'
print('rows/unique',len(rows),len({(r['chrom'],int(r['snv_pos'])) for r in rows}))
print('classes',Counter(r['methyl_class'] for r in rows))
print('rule mismatch',sum(cls(P[(r['chrom'],int(r['window_start']))])!=r['methyl_class'] for r in rows))
D=json.load(open(lnk)); E={k:Counter() for k in ('raw','powered','pbs')}
for e in D['pairs']:
    for p in (e['a'],e['b']):
        key=(e['chrom'],int(p)); E['raw'][key]+=1
        E['powered'][key]+=bool(e['powered'])
        E['pbs'][key]+=bool(e['powered'] and e['somatic_a'] and e['somatic_b'])
for name in E:
    print(name,sum(E[name][(r['chrom'],int(r['snv_pos']))]>0 for r in rows))
PY
```

輸出：stdout only。實際片段：

```text
rows/unique 30077 30077
classes Counter({'D_split_unaligned': 22454, 'B_HP_ASM': 4025, 'E_unstable': 2085, 'C_other': 765, 'A_ALLELE_clean': 748})
rule mismatch 0
raw 20535
powered 17749
pbs 17141
```

### 10.3 Coordinate crosswalk

輸入：legacy TSV 與 current `index.html#bootData`。  
命令核心：

```bash
python - <<'PY'
import csv,json,re
legacy='/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/in_progress/2026/07/20260726_methyl_allele_axis_backbone_coenrichment/methyl_class_x_linkage_annotation.tsv'
index='/bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/index.html'
s=open(index).read(); b=json.loads(re.search(r'<script type="application/json" id="bootData">(.*?)</script>',s,re.S).group(1))
l1=b['l1']; pos=[]; acc=0; prev=-1
for ci,d in zip(l1['chrom'],l1['dpos']):
    if ci!=prev: acc=d; prev=ci
    else: acc+=d
    pos.append((l1['chroms'][ci],acc))
with open(legacy,newline='') as f: L={(r['chrom'],int(r['snv_pos'])) for r in csv.DictReader(f,delimiter='\t')}
C=set(pos)
print(len(L),len(C),len(L&C),len(L-C),len(C-L),len(L|C))
print('Jaccard',100*len(L&C)/len(L|C))
PY
```

輸出：stdout only。實際片段：

```text
30077 19849 16566 13511 3283 33360
Jaccard 49.65827338129496
```

### 10.4 Report output

輸入：本稽核所有 stdout 與 source evidence。  
輸出：`InterSubMod/research/20260813_hcc1395_drilldown_validation/legacy_browser_method_audit.md`。  
驗證標準：檔案存在、HTML metadata 完整、包含 claim-to-source、P0/P1/P2、crosswalk 與 actual output snippets；不修改 `/bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/` 或 legacy report directory。

## 11. 最終判定

- **資料算術可信度：高。** 30,077 locus counts、headline raw-hit、region metrics 與 472 selection 均可由 stable raw JSON/TSV重算。
- **legacy 方法 contract 可信度：中低。** A/B priority、raw-vs-powered outcome、TP∪FP component 與 18→14 selection 未在 UI/locked JSON完整揭露。
- **legacy 科學推論可信度：低。** 單樣本、無 truth-set、control 污染、opportunity-dependent outcome、無 formal adjustment/CI/multiplicity。
- **current schema 可承接性：只可承接資訊架構，不可承接類別或數值。** 可以保留 funnel、claim ceiling、source badges、heatmap drilldown；不能把 legacy A/B/472/90.6%移植到 current v1。
- **多樣本 readiness：未達。** 先完成同 schema per-sample receipts、axis/missingness/opportunity-adjusted metrics，再談跨樣本一致性。

---

本文件為 HCC1395 全 22 常染色體的 read-only 方法稽核；未做外部 truth validation，亦未把 subset/case screenshot 當作 validation evidence。
