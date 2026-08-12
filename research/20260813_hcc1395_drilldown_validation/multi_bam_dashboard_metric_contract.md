<!--
建立時間: 2026-08-13
目標: 為一頁式多 BAM／多樣本分析總覽建立可實作、可稽核、fail-closed 的 metric 與 chart contract
處理範圍: 7 份 topology datasets（6 biological + 1 technical replicate）、HCC1395 v1/v3 bundle、HCC1395 methylation/lineage/LCA、legacy browser 方法附錄；不修改 generator/dashboard source，不生成新 bundle
關聯檔案:
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/results/metrics_audit.json
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/results/cohort_topology_metrics.csv
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/20260813_HCC1395_drilldown完整驗證與多樣本改進_01.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/legacy_browser_method_audit.md
任務類型: B — Comprehensive validation；全 7 topology datasets、全 HCC1395 v1/v3 保存 metrics，無 subset
服務目標: G3（read-level epigenetic 可解讀性）、G4（多樣本一致性與 reproducibility）、G5（可外部稽核的 metric provenance）
狀態: METRIC_CONTRACT_READY；SOURCE_IMPLEMENTATION_NOT_AUTHORIZED
-->

用 KPI Tree + Metric Contract + Evidence Ladder：**TL;DR：一頁總覽應明示「7 datasets = 6 biological samples + 1 HCC1395 technical replicate」，只把 topology 指標做七資料集比較；跨樣本 ISM/lineage 沒有等價資料時顯示 `N/A / PARTIAL`、絕不補值；6 個 biological samples 採 per-sample macro median + IQR，HCC1395 v1/v3 只放在帶 `BLOCKED` 權威標籤的診斷層（影響：高，信心：高）。**

# 多 BAM／多樣本 dashboard metric contract

> [!CAUTION]
> **Claim ceiling：descriptive observation + internal data-product QA。** 現存來源沒有同一 provenance chain 的 truth VCF／HighConf BED／TP、FP、FN，因此 precision、recall、F1 與 biological clone claims 都不得出現在 hero 或跨樣本結論。

> [!WARNING]
> **Multi-sample downstream status：PARTIAL。** 七列可比資料只涵蓋 topology；等價的 ISM/lineage drilldown 只在 HCC1395 有 v1/v3 legacy bundle。其他 dataset 必須顯示 `N/A — NO_EQUIVALENT_BUNDLE`，不可顯示 0、不可沿用 HCC1395、不可由 topology 欄位推估。

## 1. Dashboard brief 與核心決策

**主要讀者**：研究負責人、分析工程師與外部 reviewer 前的內部 gate owner。

**預設要回答的問題**：

1. 七份 dataset 中，哪些資料在相同 topology schema/model/AF basis 下可比較？
2. 6 個獨立 biological samples 的 topology coverage／identifiability 分布如何，哪一份需要優先 drill-down？
3. HCC1395 的 ISM、axis、LCA 與資產資訊可描述到哪裡，哪些 gate 阻止提升為 validation？
4. technical replicate、missing capability、不同 denominator 或不同 window 是否正在誤導比較？

**不回答**：caller accuracy、跨樣本 methylation/lineage effect、LCA causal effect、clone prevalence、跨版本 treatment effect。

### Step → Verify

1. 盤點所有保存資料表的 grain 與候選 key
   → 驗證：9 份 CSV 共 94 列；各宣告 key 重複列為 0。
2. 分離 `7-dataset topology` 與 `HCC1395-only downstream`
   → 驗證：cohort 表 7 列＝6 個 `technical_replicate=False` + 1 個 `True`；非 HCC1395 沒有等價 ISM/lineage row。
3. 先定義 numerator／denominator，再選 hero 與 chart
   → 驗證：85 個保存百分比由保存 counts 重算，與欄位值最大絕對差 `<5×10⁻⁷` percentage points。
4. 套用 aggregate／canonical／extreme／well-explained 四層閱讀
   → 驗證：每層都指向現存 sample row；technical replicate 不進 biological macro。

## 2. Source registry、grain 與 join policy

### 2.1 Canonical source precedence

1. `InterSubMod/research/20260813_hcc1395_drilldown_validation/results/metrics_audit.json`：snapshot claim ceiling、完整 row collections、`lca_ab_summary` 與既有 metric definitions。
2. 同目錄 CSV：dashboard chart/table 的 canonical rectangular inputs；不得從主報告 prose 抄數字。
3. 主報告與 `InterSubMod/research/20260813_hcc1395_drilldown_validation/legacy_browser_method_audit.md`：定義、scope、禁止推論與 evidence gate。
4. `artifact.json`：已發布 reader 的 presentation snapshot；可做 reconciliation，不作另一個 metric SoT。

### 2.2 Dataset contracts

| Dataset alias | Primary source | Grain / primary key | Rows | 跨 7 datasets？ | 用途與限制 |
|---|---|---|---:|---|---|
| `cohort_topology` | `results/cohort_topology_metrics.csv` | 1 row / dataset；PK=`sample` | 7 | **是** | 唯一可做完整七資料集比較的表；`sample` 是 dataset ID，不是已驗證的 BAM ID |
| `bundle` | `results/bundle_overview.csv` | 1 row / HCC1395 bundle；PK=`bundle` | 2 | 否 | v1/v3；`sample=HCC1395`，不可 join 後複製成兩個 biological samples |
| `coverage_k` | `results/methylation_coverage_by_k.csv` | 1 row / `bundle×k_bin` | 16 | 否 | HCC1395 ISM availability；8 個 ordered k bins / bundle |
| `axis` | `results/methylation_axis_metrics.csv` | 1 row / `bundle×scope×axis` | 20 | 否 | HCC1395 exploratory axis audit；五軸、兩 scope、兩 bundle |
| `lca_chrom` | `results/lca_ab_by_chrom.csv` | 1 row / HCC1395 autosome；PK=`chrom` | 22 | 否 | descriptive pre/post；不是 controlled A/B |
| `assets` | `results/asset_inventory.csv` | 1 row / `bundle×asset_type` | 12 | 否 | HCC1395 delivery/storage footprint |
| `integrity` | `results/input_integrity.csv` | 1 row / `bundle×path` | 6 | 否 | topology/MLHP hash 與 ISM missing/unverifiable status |
| `legacy_funnel` | `results/legacy_browser_funnel.csv` | 1 row / ordered `stage` | 4 | 否 | HCC1395 legacy selection ledger；不是 prevalence |
| `legacy_axis_crosswalk` | `results/legacy_current_axis_crosswalk.csv` | 1 row / `legacy_class` | 5 | 否 | coordinate-only descriptive cross-tab；不可作 class mapping |
| `visual_qa` | `results/legacy_browser_visual_audit.generated.json` | 1 object / surface×viewport | 4 views | 否 | HCC1395 runtime QA；不是 scientific evidence |

### 2.3 Join rules

- `cohort_topology` 不與 `bundle` 做無條件 `sample` join。HCC1395 在 `bundle` 有 v1、v3 兩列，直接 join 會把 cohort row 膨脹兩倍。
- HCC1395 downstream 區塊採獨立 data island：先選 `sample=HCC1395`，再以 `bundle` filter 查 `bundle`／`coverage_k`／`axis`／`assets`／`integrity`。
- `lca_chrom` 的 sample 是固定 HCC1395；不得靠缺少的 sample 欄位擴到其他 dataset。
- legacy 30,077 loci 與 current 19,849 loci 只有 coordinate overlap；缺 REF/ALT，且 window/statistic/taxonomy 不等價，禁止 join 後當同一 biological unit。
- `sample` 名稱不得被解析成 platform、basecaller 或 BAM identity。`HCC1395_DORADO` 的唯一可用契約是 `biological_id=HCC1395`、`technical_replicate=True`；現有表沒有 `bam_id`、input BAM SHA256 或 platform 欄位。

## 3. Identity 與 aggregation policy

| Dataset | `biological_id` | `technical_replicate` | Biological macro | Dataset chart |
|---|---|---:|---:|---:|
| COLO829 | COLO829 | False | Include | Include |
| H1437 | H1437 | False | Include | Include |
| H2009 | H2009 | False | Include | Include |
| HCC1395 | HCC1395 | False | Include | Include |
| HCC1395_DORADO | HCC1395 | **True** | **Exclude** | Include，灰色／空心標記 |
| HCC1937 | HCC1937 | False | Include | Include |
| HCC1954 | HCC1954 | False | Include | Include |

因此固定顯示：**7 datasets = 6 biological samples + 1 technical replicate**。

### Macro 規則

1. 先保留 `technical_replicate=False`，並要求每個 `biological_id` 恰有一個 canonical row；不符即 fail closed。
2. 先算每個 sample 的 numerator／denominator rate，再對 6 個 sample 取 `PERCENTILE_CONT(0.50)`；IQR 用相同線性插值口徑的 P25/P75。
3. 不做 pooled-loci micro rate；若未來要加，必須另列、另命名、同時顯示 sample weights，且不能取代 macro。
4. IQR 是 sample distribution，不是 confidence interval；`n=6` 必須和數值共置。

## 4. KPI tree

```text
可否安全比較與往下鑽？
├─ Hero：資料範圍與權威
│  ├─ 7 datasets / 6 biological / 1 technical
│  ├─ topology audit-ready datasets
│  ├─ biological macro tree coverage、unique-among-tree
│  └─ downstream status PARTIAL；canonical HCC1395 BLOCKED
├─ Diagnostic：差異從哪裡來？
│  ├─ opportunity：regions、distinct sSNV、region-position memberships
│  ├─ topology：tree coverage、unique rate 的兩種 denominator
│  └─ HCC-only：ISM availability by k、axis testability/effect、LCA、assets
├─ Guardrail：哪些數字不可升級？
│  ├─ technical replicate、source/hash、自檢、family/objective certificate
│  ├─ missing/untested/invalid、window、multiplicity、causal A/B gate
│  └─ truth-set、legacy non-equivalence、claim ceiling
└─ Detail：可回查 numerator、denominator、path、chromosome、axis 與 audit evidence
```

## 5. Hero metric contract

| ID | Metric / current snapshot | 公式與 grain | Numerator / denominator | Source fields | Scope | 建議 surface | 禁止推論 |
|---|---|---|---|---|---|---|---|
| H1 | **Dataset scope：7 / 6 / 1** | snapshot grain；`COUNT(sample)`、`COUNT(DISTINCT biological_id WHERE technical_replicate=False)`、`COUNTIF(technical_replicate)` | 7 total；6 biological；1 technical | `sample`, `biological_id`, `technical_replicate` | 7 datasets | 三值 KPI card；不可只顯示 7 | technical replicate 不是第七個獨立 biological sample |
| H2 | **Topology audit-ready：7/7** | dataset grain；row 同時滿足 `hash_match`, `sample_identity_all_match`, `receipt_all_pass`，且全 cohort 的 `schema_names/schema_versions/models/af_basis` 一致 | ready datasets / all datasets | 上述 7 欄 | 7 datasets | status card + denominator | 這是 engineering comparability，不是 biological validation、truth accuracy 或 unique clone proof |
| H3 | **Biological macro tree coverage：77.680% [P25 73.997, P75 78.600]，n=6** | sample rate=`100×tree_n/region_n`；再對 6 biological rows 取 median/IQR | 每 sample `tree_n / region_n`；macro 不合併 counts | `tree_n`, `region_n`, `tree_pct_all_regions`, `technical_replicate` | 6 biological | KPI card，IQR 作 chip；與 H4 分卡 | 不得改用 pooled `Σtree_n/Σregion_n` 取代；不代表 genome coverage 或 caller recall |
| H4 | **Biological macro unique-among-tree：62.386% [P25 46.662, P75 84.173]，n=6** | sample rate=`100×unique_tree_n/tree_n`；再取 median/IQR | 每 sample `unique_tree_n / tree_n` | `unique_tree_n`, `tree_n`, `unique_pct_among_tree`, `technical_replicate` | 6 biological | KPI card，IQR 作 chip | 不等於 unique clone prevalence、subclone purity 或 biological correctness |
| H5 | **Multi-sample ISM/lineage：PARTIAL / N/A** | categorical gate；沒有等價 per-dataset rows，不計 rate | **無合法 numerator/denominator** | 主報告 Multi-sample status；`cohort_topology` 無 ISM/lineage 欄位 | 7 datasets | warning status card | 禁止顯示 `0/7`、0%、空白即 0、或插補 HCC1395 值 |
| H6 | **Canonical downstream authority：v1 BLOCKED；v3 BLOCKED** | 1 row / HCC1395 bundle；狀態取保存值，不自行重分類 | 不適用 | `bundle`, `sample`, `validation_status_recomputed`, `selfcheck_fail`, `selfcheck_skip` | HCC1395 only | 常駐 authority ribbon；v1/v3 分列 | 不得把 source/runtime tests PASS 解讀成 bundle scientific PASS；不得外推到其他 6 datasets |

## 6. Diagnostic metric contract

| ID | Metric | 公式 / grain | Numerator / denominator | Canonical source fields | Scope | 建議圖表 | 禁止推論 |
|---|---|---|---|---|---|---|---|
| D1 | Region opportunity | dataset grain；保存 count | `region_n`；無 denominator | `sample`, `region_n`, `technical_replicate` | 7 datasets | 水平 bar；與 sSNV 分成 small multiples | region count 不是 variant count、coverage 或樣本品質排名 |
| D2 | Distinct topology sSNV | dataset grain；unique `(chrom, 1-based position)` in active positions | `distinct_ssnv_n`；無 denominator | `sample`, `distinct_ssnv_n`, `technical_replicate` | 7 datasets | 水平 bar；可 log scale但須明示 | 不等於 TP、mutation burden 或可直接跨 caller 比較的 truth count |
| D3 | Region-position memberships | dataset grain；active position 在 topology region 的 membership count | `region_position_pair_n`；無 denominator | `sample`, `region_position_pair_n` | 7 datasets | detail table / tooltip | 同一 sSNV 可屬多 region；不可當 distinct sSNV |
| D4 | Tree coverage per dataset | `100×tree_n/region_n`，dataset grain | `tree_n / region_n` | `sample`, `tree_n`, `region_n`, `tree_pct_all_regions`, `technical_replicate` | 7 datasets | 0–100% horizontal dot/bar | 不代表正確 tree 或 biological lineage coverage |
| D5 | Unique rate among tree-bearing | `100×unique_tree_n/tree_n`，dataset grain | `unique_tree_n / tree_n` | `sample`, `unique_tree_n`, `tree_n`, `unique_pct_among_tree` | 7 datasets | 與 D4 對齊的第二 panel | 不得用 `region_n` 作 denominator；不代表 unique biological clone |
| D6 | Unique rate among all regions | `100×unique_tree_n/region_n`，dataset grain | `unique_tree_n / region_n` | `unique_tree_n`, `region_n`, `unique_pct_all_regions` | 7 datasets | tooltip / detail column；不要和 D5 混名 | 不得與 D5 在未標 denominator 時比較 |
| D7 | Technical-pair descriptive delta | 同 `biological_id=HCC1395` 比較 technical row − canonical row；每 metric 各算一次 | rate 或 count 的 pairwise difference；不形成 pooled denominator | `sample`, `biological_id`, `technical_replicate`, D1–D6 fields | HCC1395 pair | 兩列 exact table或 paired dot；不要畫趨勢線 | 目前沒有 BAM hash/platform/covariate；不可稱 reproducibility、platform effect 或 DORADO effect |
| D8 | Overall ISM topology linkage | `100×methyl_topology_linked/distinct_ssnv`，bundle grain | linked topology sSNV / distinct topology sSNV | `bundle`, `methyl_topology_linked`, `distinct_ssnv`, `methyl_topology_linkage_pct` | HCC1395 only | v1/v3 exact table；不作勝負 headline | v1 ±1 kb、v3 ±5 kb；差異不是 treatment effect，summary availability 不是 significant effect |
| D9 | ISM linkage by active k | `100×numerator/denominator`，`bundle×k_bin` grain | linked distinct sSNV / topology distinct sSNV within k bin | `bundle`, `k_bin`, `numerator`, `denominator`, `percent`, `formula`, `unit` | HCC1395 only | 每 bundle 8-point ordered small-multiple line；共享 0–100% axis | 不得把 missing 當 0；v1/v3 不可作 controlled A/B；此圖是 selection-bias diagnostic |
| D10 | Axis testability / raw / BH rate | `100×count/tested_n`，`bundle×scope×axis` grain | `raw_p_le_0_05_n/tested_n`；`bh_fdr_q_le_0_05_n/tested_n` | `bundle`, `scope`, `axis`, `tested_n`, raw/BH count+percent, `multiplicity_family`, `interpretation_gate` | HCC1395 only | scope 固定後做 grouped dot/bar；cluster 灰化且帶 invalid label | raw/BH rate 不是 prevalence；不同 `tested_n` 不可只比百分比；BH 是 audit-time family |
| D11 | Axis effect distribution summary | 1 row / `bundle×scope×axis` 的五數摘要 | `effect_n` 為有 effect 值的 n；無單一通用 denominator | `source_effect_field`, `effect_unit`, `effect_n`, `effect_min/q1/median/q3/max`, `interpretation_gate` | HCC1395 only | detail table；只在相同 axis/unit 內比較 | `delta_beta`、`pseudo_F` 與未定義的 HP-fine F 不得放同軸或排名；cluster 不作 evidence |
| D12 | Bundle storage / delivery footprint | count 或 bytes，`bundle×asset_type` grain | `count`, `bytes`；無 performance denominator | `bundle`, `asset_type`, `count`, `bytes`; bundle `panel_unreported_bytes` | HCC1395 only | asset type 分面 bar；bytes 用相同單位 | bytes 增加不等於 scientific coverage、品質或 validation 增益 |
| D13 | LCA descriptive pre/post | chromosome grain；`net_new_lv_written=post_lv_written-pre_lv_written` | pre/post counts；無因果 denominator | `chrom`, pre/post/lca fields, `same_in_bam`, `same_threads`; JSON `lca_ab_summary` | HCC1395 only | gate table + 22-chrom paired plot放 detail | 7/22 input BAM 不同、22/22 threads 不同；4.969× 不得稱 LCA causal gain |

## 7. Guardrail metric contract

| ID | Guardrail / trigger | Current evidence | Source / fields | Dashboard behavior | 禁止推論 |
|---|---|---|---|---|---|
| G1 | Biological replicate policy | 7 datasets；6 biological；1 technical | `biological_id`, `technical_replicate` | macro 預設排除 technical；per-dataset chart仍保留並灰化 | 不把 HCC1395_DORADO 當獨立 biological n |
| G2 | Capability completeness | 七樣本 topology 有資料；七樣本等價 ISM/lineage 無資料 | scope matrix + 主報告 | 非 HCC sample 的 downstream card顯示 `N/A — NO_EQUIVALENT_BUNDLE` | 不補值、不沿用 HCC、不顯示 0 |
| G3 | Input provenance | topology/MLHP hash match；v1 ISM `MISSING`、v3 ISM `UNVERIFIABLE_DIRECTORY` | `integrity.status`, expected/actual SHA256/size | source status 紅/黃/綠；ISM card不可升級 | directory 存在不等於可驗 provenance |
| G4 | Bundle selfcheck authority | v1 1 FAIL + 1 SKIP；v3 1 FAIL + 0 SKIP；皆 BLOCKED | bundle selfcheck/status fields | 任一 FAIL 固定顯示 BLOCKED；不得被測試 PASS 覆蓋 | source regression PASS 不回溯改寫 immutable bundle status |
| G5 | Family / objective certificate | 7 rows `receipt_all_mutation_bearing_families_complete=False` 且 `receipt_all_units_objective_certified=False` | cohort boolean fields；count fields | 與 H2 同卡註記「technical PASS ≠ all families certified」 | `receipt_all_pass=True` 不證明所有 family/objective complete |
| G6 | Denominator visibility | tree、unique、ISM、axis rate 各有不同 denominator | numerator/denominator fields | tooltip/table 永遠同顯 numerator、denominator、unit、scope | 不顯示裸百分比；不混 `all regions` 與 `tree-bearing` |
| G7 | Missing / untested / invalid 三態 | no-summary、untested、cluster invalid 是不同狀態 | axis `interpretation_gate`; main/legacy definitions | enum：`MEASURED`、`UNTESTED`、`ABSENT`、`UNVERIFIABLE`、`INVALID`、`NOT_APPLICABLE` | 未測不得塗成陰性；invalid 不得併入非顯著 |
| G8 | Multiplicity family | BH 僅在單一 bundle×axis×declared scope audit-time 重算 | `multiplicity_family`, `scope`, `bundle` | q-rate 必須帶 family/scope；跨 sample 不 pool p/q | 不稱原 pipeline cohort-wide FDR |
| G9 | Cluster circularity | cluster axis `DOUBLE_DIPPING_INVALID_AS_EVIDENCE`，近 100% significant | `axis`, `interpretation_gate`, effect unit | 灰化、加斜線或移至 guardrail；不進 hero/aggregate | 不把顯著率當強訊號或驗證通過 |
| G10 | Version comparability | v1/v3 ISM window/source/build state不同 | `bundle`; bundle/report provenance | v1/v3 僅 descriptive side-by-side；禁止 delta badge | +109 linked loci、+0.549 pp 不代表效果或改善 |
| G11 | LCA causal gate | `same_in_bam_n=15/22`、`same_threads_n=0/22`、`causal_ab_gate=FAIL` | JSON `lca_ab_summary`; chrom gates | 所有 LCA圖加 `DESCRIPTIVE / CAUSAL GATE FAIL` | 不稱 intervention effect |
| G12 | Truth / accuracy availability | 無 TP/FP/FN/HighConf benchmark chain | 主報告 scope/omission table | precision/recall/F1 card不建立；顯示 `NOT MEASURED` | topology/ISM coverage 不能代理 caller F1 |
| G13 | Legacy non-equivalence | coordinate Jaccard 49.658%；window/statistic/taxonomy不同 | legacy crosswalk JSON + method appendix | legacy 只在 method detail；禁止全域 filter與 current聯動 | A≠ALT、B≠HP；14 cases 不是 prevalence |
| G14 | UI/runtime evidence boundary | generated QA=`DIRECT_GENERATED`、current desktop/mobile overflow false、browser errors 0；selfcheck仍 2 FAIL | generated visual audit JSON nested fields | footer顯示 runtime QA，和 science status分欄 | screenshot/無 error 不提升 scientific claim |
| G15 | BAM identity gap | source只有 dataset/sample identity與 topology output hash | cohort fields | dashboard標題稱 dataset，不宣稱每列有可驗 BAM identity | 不由 sample suffix推論 BAM、platform或basecaller effect |

## 8. Detail metric contract

| ID | Detail view | Exact fields / source | Grain | Visualization | 使用限制 |
|---|---|---|---|---|---|
| T1 | Seven-dataset topology lookup | `cohort_topology_metrics.csv` 全欄 | dataset | paginated/exact table；預設 7 列全顯示 | 不排序成 biological quality leaderboard；排序需標 metric |
| T2 | HCC1395 bundle lookup | `bundle_overview.csv` 全欄 | bundle | v1/v3 comparison table | 重複 topology 指標不可加總；window/source不同 |
| T3 | HCC1395 axis audit | `methylation_axis_metrics.csv` 全欄 | bundle×scope×axis | filterable exact table | `effect_unit`、`multiplicity_family`、`interpretation_gate` 不得隱藏 |
| T4 | HCC1395 chromosome LCA | `lca_ab_by_chrom.csv` + JSON `lca_ab_summary` | chromosome | gate matrix + exact counts | pre/post input/threads 不同，僅 descriptive |
| T5 | Asset inventory | `asset_inventory.csv` + bundle receipt byte fields | bundle×asset type | table + bytes bars | panel receipt漏記 bytes需顯示 |
| T6 | Input integrity | `input_integrity.csv` | bundle×path | status table；path、expected/actual hash/size可展開 | 空 hash不等於 match；missing/unverifiable分開 |
| T7 | Legacy selection funnel | `legacy_browser_funnel.csv`: `stage,n,pct_of_first_stage,pct_of_previous_stage` | ordered stage | stage bars + exact labels；不用裝飾性 funnel geometry | `14/472` 是 displayed selection，不是 population prevalence |
| T8 | Legacy/current crosswalk | `legacy_current_axis_crosswalk.csv` + `legacy_browser_crosswalk.json` | legacy class / snapshot | exact table；若做 heatmap需先以 `shared_loci` 正規化並保留 n | coordinate-only；不能產生一對一 class conversion |
| T9 | Browser visual QA | generated visual audit JSON：`evidence_status`, `viewport`, `widths.horizontalOverflow`, `browserErrors`, `selfcheck` | surface×viewport | engineering footer/table + screenshot link | QA viewport/runtime，不是數據或科學驗證 |

## 9. Chart contract 與欄位對照

| Chart ID / placement | Analytical question | Family / concrete form | 現存欄位（全部可直接對到） | Data sufficiency / fallback | Supported takeaway | 不支援的 claim |
|---|---|---|---|---|---|---|
| C1 / primary | 各 dataset 的 tree coverage 與 conditional uniqueness 差多少？ | 兩個對齊的 horizontal dot/bar panels；共用 0–100% scale | `sample`, `biological_id`, `technical_replicate`, `tree_n`, `region_n`, `tree_pct_all_regions`, `unique_tree_n`, `unique_pct_among_tree` | 7 categories 足夠；DORADO空心/灰色，HCC1395描邊；tooltip顯 n/d | topology heterogeneity、technical非獨立 | biological quality ranking、clone truth |
| C2 / diagnostic | dataset opportunity 是否差很多？ | `region_n` 與 `distinct_ssnv_n` 分開的 horizontal bar small multiples | `sample`, `technical_replicate`, `region_n`, `distinct_ssnv_n` | 7 rows；若長尾壓縮可明示 log scale；不可雙軸 | H2009 opportunity/count scale明顯較大 | mutation burden、recall、樣本品質 |
| C3 / four-layer strip | aggregate、canonical、extreme、well-explained 如何同時看？ | 四張 evidence cards，不畫折線 | C1/C2欄位 + `validation_status_recomputed` | 直接讀指定 rows；每卡附 selection rule | 見樹也見林，且不讓單例取代 cohort | 四張卡不是四個獨立 biological samples |
| C4 / HCC-only | ISM availability 是否隨 active k下降？ | 每 bundle 一個 8-point ordered line small multiple；相同 y-scale | `bundle`, `k_bin`, `numerator`, `denominator`, `percent`, `unit`, `formula` | 每 bundle正好8點；若缺 bin改 stage table，不連線插值 | 高 k availability較低，需 missingness guardrail | v1/v3 treatment comparison、methyl effect |
| C5 / HCC-only | 各 axis 有多少 tested、raw與audit-BH significant？ | grouped dot/bar；固定一個 `scope`，facet `bundle` | `bundle`, `scope`, `axis`, `tested_n`, `raw_p_le_0_05_n/pct`, `bh_fdr_q_le_0_05_n/pct`, `interpretation_gate` | 5 axes×2 bundles；cluster灰化；tooltip顯 n/d | denominator差異與 invalid axis | prevalence、causal/signature strength |
| C6 / technical control | HCC1395 canonical/technical row差多少？ | exact two-row table或 paired dot；不畫 trend | `sample`, `biological_id`, `technical_replicate`, D1–D6 fields | 只有2 rows，低於 scatter/trend門檻；table優先 | technical row非獨立且觀察值相近/不同可查 | reproducibility estimate、platform effect |
| C7 / HCC detail | v1/v3 delivery footprint由哪些資產組成？ | asset type faceted bars；count、bytes分圖 | `bundle`, `asset_type`, `count`, `bytes`; `panel_unreported_bytes` | 6 asset types×2 bundles | storage/delivery trade-off | science quality或validation gain |
| C8 / HCC detail | LCA pre/post gate在哪些 chromosome失敗？ | 22-row gate heatmap/table；可附 pre/post paired dots | `chrom`, `pre_lv_written`, `post_lv_written`, `net_new_lv_written`, `same_in_bam`, `same_threads`, `sample_identity_pass` | 22 chromosomes足夠；gate先於delta | input/parameter non-equivalence | causal LCA gain |
| C9 / method detail | legacy population如何縮到displayed cases？ | ordered stage bars + exact percentages | `stage`, `n`, `pct_of_first_stage`, `pct_of_previous_stage` | 4 ordered stages；bar比裝飾性漏斗更誠實 | selection attrition | case prevalence/representativeness |

### 明確不採用的圖表

- **7-point scatter**：低於建議的 12–20 observations，無法可靠顯示 relationship；C1/C2 使用 bars/dots。
- **sample 折線圖**：sample 無自然順序，不得製造「趨勢」。
- **v1/v3 單一差值 waterfall**：版本不受控，driver 也不具可加總性。
- **pie**：7 dataset 身分與多個 denominator 不適合 part-to-whole。
- **沒有 denominator 的 significant-rate leaderboard**：tested n 差異大，cluster 又 invalid。
- **同軸 effect plot**：`delta_beta`、`pseudo_F`、HP-fine stored F 不同單位。

### 視覺 encoding

- biological datasets：同一藍色 root；HCC1395 用外框/直接標籤，不靠另一個彩虹色。
- HCC1395_DORADO：中性灰 + 空心 marker + `technical` 文字標籤；不可只靠色彩。
- invalid／blocked：灰化、斜線或文字 badge；不預設用紅綠二分。
- 所有 percentage chart 顯示 `%`，來源 CSV 是 0–100 scale；不得再次乘 100。
- title 保持描述性；subtitle 必含 grain、denominator、`n=6 biological` 或 `HCC1395 only`。

## 10. 四層展示 contract

| Layer | Selection rule | Current row / exact snapshot | 為何放這層 | 限制 |
|---|---|---|---|---|
| **Aggregate** | `technical_replicate=False` 的 6 個 canonical biological rows；per-sample macro median/IQR | tree coverage 77.680% [73.997, 78.600]；unique/tree 62.386% [46.662, 84.173]；n=6 | 顯示 cohort 中心與spread，不讓大樣本支配 | IQR非CI；不作 pooled loci inference |
| **Canonical** | 任務 focal 且唯一有 downstream v1/v3 bundle 的 biological row | HCC1395：11,590 regions；9,130 trees；78.775% tree coverage；7,047 unique；77.185% unique/tree；v1/v3皆BLOCKED | 提供完整 drilldown 的閱讀入口 | canonical 是 evidence availability選擇，不代表 cohort中位或最佳樣本 |
| **Extreme observed** | 在 6 biological rows 中同時為 tree coverage與unique/tree最小者 | H2009：36,042 regions；23,128 trees；64.170% coverage；8,161 unique；35.286% unique/tree。相反端可註記 HCC1954 unique/tree=90.579% | 優先檢查 opportunity、resource guard與identifiability | n=6；稱 observed extreme，不宣稱正式統計 outlier或根因 |
| **Well-explained control** | 與 canonical 相同 `biological_id` 且 `technical_replicate=True` | HCC1395_DORADO：6,865 regions；77.320% coverage；78.466% unique/tree；相對 HCC1395為 −1.455 pp / +1.281 pp | 已知非獨立身分可解釋「為何不進 biological macro」 | 只解釋 replicate policy；沒有 BAM/platform欄位，不能解釋差異成因 |

## 11. One-page layout 與 interaction semantics

### 預設畫面

1. **Authority banner**：`PARTIAL multi-sample` + `HCC1395 BLOCKED` + claim ceiling。
2. **Hero strip**：H1–H5；H6 置於 authority banner，不讓使用者滾動後才看到。
3. **Cohort primary**：C1；下方 C2 或 compact opportunity table。
4. **Four-layer strip**：aggregate／canonical／extreme observed／well-explained control。
5. **HCC-only diagnostic**：C4 + C5；非 HCC sample 時不顯空圖，改顯 `N/A — NO_EQUIVALENT_BUNDLE`。
6. **Guardrail footer + detail drawer**：G1–G15 狀態；T1–T9 可展開。

### Filters

| Filter | Default | Applies to | Fail-closed behavior |
|---|---|---|---|
| Dataset | All 7 visible | C1/C2/T1 | macro仍只用6 biological；選 technical時加non-independent banner |
| Aggregate population | Biological only | H3/H4/aggregate layer | technical不得透過「All」偷偷進 macro；另提供dataset view |
| Sample | HCC1395 for downstream | HCC-only區塊 | 非HCC顯示N/A，不借用HCC source |
| Bundle | v1、v3 side-by-side | HCC bundle/ISM/axis/assets | 只有sample=HCC1395才啟用；禁止當全 cohort version filter |
| Axis scope | `topology_linked_distinct_ssnv` | C5/T3 | `all_summary_rows`另選；兩scope不得疊在同一rate |
| Axis | HP / hpfine / allele / cluster / lineage | C5/T3 | cluster永遠保留invalid gate；不能由filter移除警告 |

### Missing-value semantics

- `0`：只在已知 denominator > 0 且保存 numerator=0 時使用。
- `N/A — NO_EQUIVALENT_BUNDLE`：dataset 沒有同 contract 的 source row。
- `UNTESTED`：有 summary，但該 axis 無合法 test；不是 non-significant。
- `UNVERIFIABLE`：檔案/目錄存在但 hash/provenance 無法核對。
- `MISSING`：receipt 所指 source 已不存在。
- `INVALID`：數字可存在，但設計上不可作 evidence，例如 cluster double-dipping。
- `NOT_APPLICABLE`：metric 與當前 dataset/surface無關。

## 12. Data-quality findings 與 implementation gates

| Finding | Evidence | Severity / confidence | Implementation gate |
|---|---|---|---|
| 宣告 key 無重複 | 9 CSV、94 rows；各 composite key duplicate=0 | Low risk / High | 可作 snapshot，但仍需 CI key uniqueness test |
| 保存 rate 算術一致 | tree/unique/ISM/axis 共 85 checks；最大差 `<5×10⁻⁷` pp | Low risk / High | renderer可讀保存率，但建議以 counts重算並reconcile |
| HCC-only join會一對多 | `bundle_overview` 對 HCC1395 有 v1/v3兩列 | High / High | data islands；禁止裸 `sample` join |
| Bundle provenance缺欄 | `generator_commit`, `receipt_schema_name/version` 兩列皆空 | High / High | H6永遠BLOCKED；不得標 reproducible build |
| ISM input不封口 | v1=`MISSING`；v3=`UNVERIFIABLE_DIRECTORY` | High / High | HCC methyl只可legacy observation |
| 跨樣本 capability coverage不足 | cohort只有topology；無6個非HCC等價ISM/lineage rows | Critical for multi-omics / High | H5顯N/A/PARTIAL；不得算coverage rate |
| BAM identity不可驗 | cohort缺 `bam_id`, input BAM SHA256, platform | High for multi-BAM / High | UI稱dataset；下一版manifest必補 |
| technical replicate非獨立 | HCC1395_DORADO biological_id重複且technical=true | High / High | 只作control，不進n=6 macro |
| legacy/current schema不等價 | cohort/window/statistic/taxonomy不同；coordinate-only overlap | High / High | legacy只放method detail |

### 下一版資料契約必補欄位（目前不得假裝已有）

- `dataset_id`, `biological_id`, `replicate_id`, `replicate_type`, `bam_uri`, `bam_sha256`, `platform`, `basecaller`, `coverage_summary`。
- 每 sample 的 capability matrix：`topology/MLHP/ISM/lineage/LCA/panels/IGV` 各自 `AVAILABLE/PARTIAL/ABSENT/UNVERIFIABLE`、source URI、hash、schema/version。
- 固定 ISM `window_bp`、axis definition/effect unit、group n、validity、missing reason、multiplicity family/q-value provenance。
- truth benchmark contract：caller/callset、truth VCF、HighConf BED、command、TP/FP/FN、precision/recall/F1 與 receipt。
- 每次 dashboard snapshot 的 `generated_at`, source revision, extract hash；目前只有 bundle-level `built_at`，不足以代表整頁 freshness。

## 13. 完整 BAM 狀況仍需新增的資料欄

下列欄位在本次 `results/*.csv/*.json` dashboard snapshot **全部是 `NOT_COLLECTED`**。這裡定義的是未來採集契約，不代表數值為 0，也不代表 workspace 其他位置一定沒有原始檔。實作前必須由逐 BAM receipt／QC 輸出提供可驗 source、grain、denominator 與 hash。

> [!IMPORTANT]
> Dashboard 現在只能顯示 `NOT_COLLECTED`；不得產生 0、0%、空 bar、零長度 phase block 或「無 methylation」等視覺。只有 source row 實際存在且通過 identity/hash gate 後，才可把狀態改為 `MEASURED` 並啟用圖表。

| ID / 面向 | 未來必備欄位 | 定義、grain 與 denominator | Current status（全 7 datasets） | 未來可用 surface | 現在禁止的替代推論 |
|---|---|---|---|---|---|
| B1 / BAM identity + reference | `dataset_id`, `bam_uri`, `bam_sha256`, `bai_uri`, `bai_sha256`, `reference_build`, `reference_fasta_uri`, `reference_sha256`, `read_group_samples` | 1 row / dataset×BAM；BAM、index、reference各自 content hash；sample/RG identity另存，不由檔名解析 | **NOT_COLLECTED** | identity/provenance table；hash gate | topology output hash不可代替 input BAM/reference hash；`sample` suffix不可代理platform/basecaller |
| B2 / mapped + primary reads | `reads_total`, `reads_primary`, `reads_secondary`, `reads_supplementary`, `reads_duplicate`, `reads_mapped`, `reads_primary_mapped`, `flagstat_tool_version` | BAM grain；primary mapping rate=`reads_primary_mapped/reads_primary`；overall mapping rate=`reads_mapped/reads_total`；兩者分列且 N/D 同顯 | **NOT_COLLECTED** | exact KPI + 7-dataset horizontal bars，等資料完整後才啟用 | 不得用 region/sSNV count代理 read count；missing mapping rate不得顯示 0% |
| B3 / depth + IQR / breadth | `depth_scope`, `scope_bases`, `depth_mean`, `depth_p25`, `depth_median`, `depth_p75`, `depth_p95`, `covered_bases_ge_1/10/20/30`, `depth_tool_version`, `depth_command` | dataset×declared scope；IQR=`P25–P75`；breadth@x=`covered_bases_ge_x/scope_bases`；必須凍結 WGS/HC BED/callable scope與是否含duplicate | **NOT_COLLECTED** | median+IQR dot/interval；breadth bars | 不得從 BAM bytes、read count、topology regions推估 depth；不同 scope不可同圖 |
| B4 / read length + N50 | `read_length_population`, `eligible_read_n`, `eligible_bases`, `read_length_p25/median/p75`, `read_n50_bp`, `read_length_max`, `read_length_tool_version` | dataset×明示 read population；N50 是使累積 bases達50%的 read length，不是 median；須固定用query length或aligned length、primary/supplementary policy | **NOT_COLLECTED** | N50 KPI + length distribution summary；有原始 bins才畫histogram | 不得由 IGV asset、coverage 或 file size反推 read N50 |
| B5 / HP/PS tag coverage | `eligible_primary_mapped_reads`, `hp_tag_read_n`, `ps_tag_read_n`, `hp_and_ps_read_n`, `valid_hp_value_read_n`, `valid_ps_value_read_n`, `hp_tag_rate`, `ps_tag_rate`, `hp_ps_joint_rate`, `tag_parser_version` | dataset×BAM；各 rate 的 denominator固定為eligible primary mapped reads；HP/PS存在、值合法、共同存在分開 | **NOT_COLLECTED** | numerator/denominator KPI + paired bars | `ranked_n`、tree coverage或phased variant count不可代理read-level tag coverage；missing不得稱 unphased |
| B6 / phase blocks | `phased_variant_n`, `phase_block_n`, `phase_block_total_span_bp`, `phase_block_n50_bp`, `phase_block_median_bp`, `phase_block_max_bp`, `phase_block_scope`, `phasing_tool_version`, `phasing_command` | dataset×VCF/BAM phasing run；block N50以block span加權且定義必須固定；autosome/HC BED scope另列 | **NOT_COLLECTED** | block N50/median KPI + distribution；需有per-block rows才畫histogram | topology tree/region不是phase block；無truth時不得顯示switch error或phasing accuracy=0 |
| B7 / MM/ML + CpG | `eligible_primary_mapped_reads`, `mm_tag_read_n`, `ml_tag_read_n`, `valid_mm_ml_read_n`, `mm_ml_parse_error_n`, `cpg_call_n`, `unique_cpg_site_n`, `callable_cpg_site_n`, `methylated_cpg_call_n`, `beta_median/p25/p75`, `modkit_or_parser_version`, `mod_code`, `ml_threshold` | dataset×BAM×modification code；tag coverage以eligible reads為denominator；CpG call/site/read三種grain分開；beta需明示threshold與callable denominator | **NOT_COLLECTED** | tag coverage bars、CpG count table、beta interval；只在同mod/threshold下比較 | ISM summary rows或linked sSNV不可代理MM/ML/CpG完整度；missing tag不可當unmethylated |
| B8 / KDE provenance | `kde_status`, `kde_corrected`, `kde_method`, `kde_software_version`, `kde_parameters_json`, `kde_parameters_sha256`, `kde_input_uri/sha256`, `kde_output_uri/sha256`, `kde_completed_at` | dataset×methylation preprocessing run；`kde_corrected`須為明示 enum/boolean，不能由資料分布猜測 | **NOT_COLLECTED** | provenance status table；無合格數值圖 | 不得從 v1/v3、CpG數或effect distribution推斷KDE-corrected；unknown不等於False |
| B9 / VCF + truth benchmark | `callset_vcf_uri/sha256`, `caller`, `caller_version`, `caller_mode`, `truth_vcf_uri/sha256`, `highconf_bed_uri/sha256`, `reference_sha256`, `benchmark_command`, `benchmark_tool/version`, `tp`, `fp`, `fn`, `precision`, `recall`, `f1`, `benchmark_scope`, `benchmark_receipt_sha256` | dataset×callset×truth×scope；precision=`TP/(TP+FP)`、recall=`TP/(TP+FN)`、F1 harmonic mean；所有 counts與region policy同一receipt | **NOT_COLLECTED** | truth readiness status；來源完整後才建立accuracy cards/plots | topology distinct sSNV、ISM linkage或selfcheck PASS不可代理TP/FP/FN/F1；缺truth不得顯示0分 |
| B10 / runtime + resources | `run_id`, `command_argv`, `code_commit`, `working_tree_dirty`, `container_image_digest`, `host`, `started_at`, `completed_at`, `exit_code`, `wall_seconds`, `cpu_seconds`, `max_rss_bytes`, `read_bytes`, `write_bytes`, `output_bytes`, `thread_count` | 1 row / run；資源必須和確切command/input/output hashes綁定；failed/aborted run保留status，不與成功run取平均 | **NOT_COLLECTED** | run table、wall/RSS/output cost bars；只比較同scope run | bundle asset bytes不可代理runtime/RSS；缺run receipt不得稱效能改善 |

### BAM readiness gate

完整 BAM 狀況區塊只能在下列條件同時成立後啟用：

1. 每個要比較的 dataset 都有 B1 identity/reference row，BAM、BAI、reference hash非空且可驗。
2. 每一 metric family 都有明確的 eligible population、numerator、denominator、tool/version與scope。
3. 同一張跨樣本圖的 reference、region scope、read policy、tag policy與工具版本相容；否則分面或標 `NON_COMPARABLE`。
4. 任一 dataset 缺 row 時維持 `NOT_COLLECTED`；不能因其他 6 份有值而插補或以 cohort median代替。
5. VCF/truth、KDE、MM/ML、phasing各自有獨立 readiness gate；BAM identity通過不會自動讓下游 family變成 `MEASURED`。

## 14. 驗證命令、I/O 與實際結果

**輸入路徑**：

- `InterSubMod/research/20260813_hcc1395_drilldown_validation/results/*.csv`
- `InterSubMod/research/20260813_hcc1395_drilldown_validation/results/*.json`
- `InterSubMod/research/20260813_hcc1395_drilldown_validation/20260813_HCC1395_drilldown完整驗證與多樣本改進_01.md`
- `InterSubMod/research/20260813_hcc1395_drilldown_validation/legacy_browser_method_audit.md`

**執行命令（formula reconciliation）**：

```bash
python3 - <<'PY'
from pathlib import Path
import csv
b = Path('research/20260813_hcc1395_drilldown_validation/results')
checks = []
with (b/'cohort_topology_metrics.csv').open(newline='') as f:
    for r in csv.DictReader(f):
        checks += [
            (100*int(r['tree_n'])/int(r['region_n']), float(r['tree_pct_all_regions'])),
            (100*int(r['unique_tree_n'])/int(r['tree_n']), float(r['unique_pct_among_tree'])),
            (100*int(r['unique_tree_n'])/int(r['region_n']), float(r['unique_pct_all_regions'])),
        ]
with (b/'bundle_overview.csv').open(newline='') as f:
    for r in csv.DictReader(f):
        checks += [
            (100*int(r['regions_with_tree'])/int(r['regions']), float(r['tree_coverage_pct_all_regions'])),
            (100*int(r['unique_best_tree'])/int(r['regions_with_tree']), float(r['unique_pct_among_tree'])),
            (100*int(r['unique_best_tree'])/int(r['regions']), float(r['unique_pct_all_regions'])),
            (100*int(r['methyl_topology_linked'])/int(r['distinct_ssnv']), float(r['methyl_topology_linkage_pct'])),
        ]
with (b/'methylation_coverage_by_k.csv').open(newline='') as f:
    for r in csv.DictReader(f):
        checks.append((100*int(r['numerator'])/int(r['denominator']), float(r['percent'])))
with (b/'methylation_axis_metrics.csv').open(newline='') as f:
    for r in csv.DictReader(f):
        checks += [
            (100*int(r['raw_p_le_0_05_n'])/int(r['tested_n']), float(r['raw_p_le_0_05_pct'])),
            (100*int(r['bh_fdr_q_le_0_05_n'])/int(r['tested_n']), float(r['bh_fdr_q_le_0_05_pct'])),
        ]
deltas = [abs(a-b) for a,b in checks]
print(f'formula_checks={len(checks)} mismatch_gt_1e-5={sum(d>1e-5 for d in deltas)} max_abs_delta={max(deltas):.12g}')
PY
```

**執行命令（key uniqueness + biological identity）**：

```bash
python3 - <<'PY'
from pathlib import Path
import csv, collections
b = Path('research/20260813_hcc1395_drilldown_validation/results')
keys = {
    'asset_inventory.csv': ['bundle', 'asset_type'],
    'bundle_overview.csv': ['bundle'],
    'cohort_topology_metrics.csv': ['sample'],
    'input_integrity.csv': ['bundle', 'path'],
    'lca_ab_by_chrom.csv': ['chrom'],
    'legacy_browser_funnel.csv': ['stage'],
    'legacy_current_axis_crosswalk.csv': ['legacy_class'],
    'methylation_axis_metrics.csv': ['bundle', 'scope', 'axis'],
    'methylation_coverage_by_k.csv': ['bundle', 'k_bin'],
}
total = duplicates = 0
cohort = None
for name, key in keys.items():
    with (b/name).open(newline='') as f:
        rows = list(csv.DictReader(f))
    total += len(rows)
    counts = collections.Counter(tuple(row[k] for k in key) for row in rows)
    duplicates += sum(n - 1 for n in counts.values() if n > 1)
    if name == 'cohort_topology_metrics.csv':
        cohort = rows
biological = sum(row['technical_replicate'] == 'False' for row in cohort)
technical = sum(row['technical_replicate'] == 'True' for row in cohort)
print(f'CSV key audit: {len(keys)} files, {total} rows, duplicate key rows={duplicates}')
print(f'cohort identity: datasets={len(cohort)}, biological={biological}, technical={technical}')
PY
```

**實際輸出片段**：

```text
formula_checks=85 mismatch_gt_1e-5=0 max_abs_delta=4.97091839691e-07
CSV key audit: 9 files, 94 rows, duplicate key rows=0
cohort identity: datasets=7, biological=6, technical=1
```

**輸出路徑**：`InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_dashboard_metric_contract.md`

## 15. Acceptance checklist

- [x] 明示 7 datasets = 6 biological samples + 1 technical replicate。
- [x] macro 一律排除 HCC1395_DORADO；dataset chart仍保留該列。
- [x] 等價 ISM/lineage 缺失顯示 N/A/PARTIAL，不補值、不當 0。
- [x] Hero／diagnostic／guardrail／detail metrics 均有公式、grain、N/D、來源、scope 與禁止推論。
- [x] 每張建議 chart 均只引用現存欄位；稀疏 chart有 table/bar fallback。
- [x] aggregate／canonical／extreme observed／well-explained control 四層皆有 selection rule 與現存 row。
- [x] HCC1395-only metrics 不外推其餘 6 datasets；v1/v3不當 controlled A/B。
- [x] legacy/current 不作 class mapping；visual QA 不作 scientific evidence。
- [x] 完整 BAM identity/reference、alignment、depth/read N50、HP/PS、phase block、MM/ML/CpG、KDE、VCF/truth、runtime/resources 全部標 `NOT_COLLECTED`；不畫 0。
- [x] 未修改 generator、dashboard source、artifact 或其他報告。
