<!--
建立時間：2026-07-11 15:35–15:40 +08:00
目標：G3（read-level epigenetic 證據邊界）／G5（可外部驗證的業界級證據鏈）
任務類型：E Hotfix / Bugfix audit（唯讀；未修改 HTML、generator 或研究資料）
處理範圍：frozen HTML、companion JSON、historical raw SoT、current production status、generator provenance
關聯檔案：InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.standalone.html
-->

# 2026-07-11 分層重建全景報告：Data Quality / Provenance Audit

用 Evidence-to-Claim：先固定 grain 與來源，再逐項比較 headline、table、chart，最後把「顯示錯誤」和「來源本身尚不完整」分開判定。

## TL;DR

**Frozen HTML 的歷史數值層通過：8/8 headline cards（9 個顯示值）、14/14 charts、14/14 fallback tables、359 chart-table rows／3,472 cells 均與 companion data 和 raw JSON 重算一致；沒有發現總數、分母、百分比或類別守恆錯誤。科學發布仍為 NO-GO：0 P0、4 P1、3 P2（影響：高，信心：高）。**

最重要的問題不是 `47,377`、`39,885` 或 `84.19%` 算錯，而是：

1. HTML 的 `0/7 PASS` 對其 07:31 cutoff 與 `PASS_ONLY_ABORTED_v1` 是正確的，但目前 active normalized raw-all producer 已是 `3/7 PASS + H1437 START`；generic「Production probe」標籤會把兩個 evidence tracks 混在一起。
2. Historical outputs 雜湊 7/7 通過，但 `code.sha256` 實際是 **0/6 executable scripts match**；表面上的 `1/7` 是一個 input-manifest entry，不是任何 producer script。
3. HCC1395 有真實且已揭露的 `W_ret 7,928 → W_tree 7,927` 缺口：少 1 region／4 sSNV。
4. 6/7 historical HP/PS inputs 受 truth-BED conditioning；read-AF 0/7、L3 methylation 0/7、orthogonal truth 缺失。這些是 source completeness / eligibility blocker，不是 HTML 計算錯誤。

## 1. Verdict 與嚴重度

| 層級 | 數量 | 判定 |
|---|---:|---|
| P0 Critical | 0 | 沒有未受保護的錯誤 scientific promotion；頁首與 validation 均明示 `PUBLICATION NO-GO`。 |
| P1 High | 4 | Production track freshness、historical truth-BED eligibility、producer-code reproducibility、缺失的 read-AF/L3/orthogonal evidence。 |
| P2 Medium | 3 | HCC1395 1-region/4-site gap、per-sample provenance 欄位不足、audit 期間 current generator 與 frozen artifact 漂移。 |
| Numeric display | PASS | Frozen scope 內的 headline/table/chart 值與比例均可重算。 |
| Scientific release | **NO-GO** | 只能作 historical engineering observation；不可當 biological confirmation、clone/subclone confirmation 或 external handoff evidence。 |

嚴重度定義：P0＝會無保護地導向錯誤科學結論；P1＝會影響 release、status decision 或 exact reproducibility；P2＝有界且已揭露、但仍須修復的完整性／metadata 問題。

## 2. Scope、grain 與正確來源

### 2.1 Audited inputs

| Evidence role | Input path | Audit role |
|---|---|---|
| Frozen display | `InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.standalone.html` | 使用者實際看到的 headline、tables、charts、claim boundary。 |
| Companion snapshot | 同目錄 `20260711_分層重建全景數據觀察_01.data.json` | 報告重算後的 structured snapshot。 |
| Chart/source contract | 同目錄 `20260711_分層重建全景數據觀察_01.source_notes.json` | 14 charts 的 question、grain、denominator、source contract。 |
| Artifact validation | 同目錄 `20260711_分層重建全景數據觀察_01.validation.json` | Frozen artifact hash、artifact PASS 與 scientific NO_GO 的分流。 |
| Historical raw SoT | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/` | 7 sample dirs、35 MLHP parts、7 region views、7 layered reconstructions、verification、SHA manifests。 |
| Historical eligibility | 上述 root 的 `VALIDATION_SCOPE.md`；`InterSubMod/research/20260710_layered_reconstruction_v2/clairs_longphase_ssnv_contract_audit.json` | 證明 6/7 truth-BED conditioning，唯一例外為 COLO829。 |
| Frozen probe source | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_production_sidecars_PASS_ONLY_ABORTED_v1/` | HTML 07:31 cutoff 的 `0/7` 與 `E_METHOD_SCOPE_PASS_ONLY`。 |
| Current production source | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/` | Audit 時真正 active 的 normalized raw-all producer status。 |
| Builder provenance | `InterSubMod/research/20260710_layered_reconstruction_v2/scripts/build_layered_observation_report.py` | 現行 workspace generator；不是 immutable archived copy。 |

### 2.2 Canonical grains

| Layer | Canonical grain | Denominator / total | Audit rule |
|---|---|---:|---|
| Dataset | dataset row | 7 rows | HCC1395 與 HCC1395_DORADO 是同一 biological sample 的兩個 dataset rows。 |
| Biology | biological sample | 6 unique IDs | 任何 inferential statement 不可把 7 rows 當 7 個獨立樣本。 |
| Site funnel | caller-PASS biallelic sSNV | `U=568,080` | 五個 mutually exclusive exits 必須回加 U。 |
| Region/group | MLHP group / emitted region | `W_pre=48,963`、`W_ret=48,960`、`W_tree=48,959` | 不是 sSNV grain；不可把 `W/U` 當 site retention。 |
| Primary region | emitted region with HP1/HP2 primary | `W_primary=47,377` | Topology / hidden-region composition 的 parent。 |
| Reconstruction unit | primary HP lineage | `68,544` | L1 outcomes、solver-cap、unit-level determinacy 的 parent。 |
| Read evidence | retained pattern observation | per dataset | full + partial-X；不是 read-AF 或 CCF。 |

Grain audit 結果：7 dataset keys 全唯一；biological IDs 為 6 個，唯一重複是 `HCC1395 × 2`。HTML 多處明示「7 dataset rows / 6 biological samples」及 descriptive aggregate，因此這一項 **PASS**；若後續做跨樣本 inference，仍須以 6 biological samples 或明確 technical-replicate model 重算。

## 3. Headline item-by-item comparison

| Headline / KPI | HTML | Independent recomputation | Source | Result |
|---|---:|---:|---|---|
| Historical snapshot outputs | 7 | 7/7 sample `output.sha256` manifests PASS（49/49 listed files match） | historical sample roots | PASS；只證明 output bytes，未證明 current code。 |
| Biological samples | 6 | 7 dataset rows、6 unique biological IDs；`HCC1395` duplicated | `input_manifest.json` + raw verification | PASS。 |
| Production probe PASS | 0/7 | Frozen referenced `PASS_ONLY_ABORTED_v1` 為 0/7；audit 時 active raw-all root 為 3/7 | frozen probe vs current producer | **Frozen cutoff PASS；current/global label stale，見 DQ-01。** |
| Scientific release | NO-GO | 6/7 truth-BED、0/7 read-AF、0/7 L3、HCC gap、code drift、無 orthogonal truth | validation + scope SoT | PASS。 |
| W_primary | 47,377 | 26,210 single-primary + 21,167 double-primary | 7 region views | PASS。 |
| Complete candidate regions | 39,885 | `47,377 − 7,492 incomplete = 39,885` | region/topology recomputation | PASS。 |
| Complete % | 84.19% | `39,885 / 47,377 = 84.18642% → 84.19%` | same | PASS；rounding correct。 |
| C=1, Topo=1 | 10,832 | 7-sample raw region recomputation = 10,832 | region views + layered units | PASS。 |
| C>1, Topo>1 | 17,909 | 7-sample raw region recomputation = 17,909 | same | PASS。 |

## 4. Total、percentage 與 conservation checks

所有以下等式都使用相同 grain；無跨 site／region／unit 混加：

```text
Site U:
568,080 = 112,870 out-of-scope
        + 49,978 singleton
        + 222,821 cap-excluded
        + 11 read-unsupported
        + 182,400 retained

Autosomal branch:
455,210 = 49,978 singleton + 405,232 k>1 pre-cap sites
405,232 = 222,821 cap-excluded + 11 read-unsupported + 182,400 retained

Region bridge:
W_all_pre 98,941 = 49,978 k=1 components + 48,963 k>1 groups
W_pre 48,963 → W_ret 48,960 → W_tree 48,959 → W_primary 47,377

HP multiplicity:
W_tree 48,959 = 26,210 single-primary + 21,167 double-primary + 1,582 no-primary
W_primary 47,377 = 26,210 + 21,167

Region topology:
47,377 = 39,885 complete + 7,492 incomplete
47,377 = 10,832 (C=1,Topo=1) + 11,144 (C>1,Topo=1)
         + 17,909 (C>1,Topo>1) + 7,492 incomplete
47,377 = 7,325 hidden-zero + 32,560 hidden-positive + 7,492 incomplete

Primary units:
68,544 = 24,409 determined + 36,462 ambiguous + 7,542 solver-capped + 131 recurrence
```

Percentage checks：所有 chart shares 均在 `[0,100%]`；stacked parents 在未四捨五入前等於 100%；顯示誤差只來自 1 或 2 位小數四捨五入。沒有發現 numerator > denominator、負數、impossible `C=1 / Topo>1` 或 denominator substitution。

## 5. Chart item-by-item comparison

Audit 對每張圖做三向比較：raw/companion 重算 → exact-data table → static fallback table。14/14 的 fallback 與 exact table 完全相同，且 exact table 與重算值相同；合計 359 rows／3,472 cells，mismatch `0`。

| Chart id | Grain / plotted denominator | Compared rows | Key invariant | Result |
|---|---|---:|---|---|
| `site-funnel-composition` | sSNV / PASS U | 7 | 5 exits sum to U per dataset | PASS |
| `region-stage-counts` | explicit nested region/group stages | 7 | `W_pre ≥ W_ret ≥ W_tree ≥ W_primary` | PASS |
| `k-distribution` | MLHP groups / W_ret | 7 | k=2…8 counts sum to W_ret | PASS |
| `read-pattern-composition` | pattern observations | 7 | full + partial-X = retained pattern observations | PASS |
| `hp-multiplicity` | emitted regions / W_tree | 7 | single + double + none = W_tree | PASS |
| `hp-h3-strata` | emitted regions / W_tree | 48 | 6 displayed strata sum to W_tree；aggregate row也閉合 | PASS |
| `candidate-combinations-single` | one-primary regions | 64 | C=1…6, >6, incomplete sum to single-primary parent | PASS |
| `candidate-combinations-double` | two-primary regions | 64 | joint-C bins sum to double-primary parent | PASS |
| `region-topology-classes` | primary regions / W_primary | 96 | 3 feasible states + incomplete sum to parent；impossible state=0 | PASS |
| `hidden-node-classes` | primary regions / W_primary | 24 | zero + positive + incomplete sum to parent | PASS |
| `tree-outcomes` | primary reconstruction units | 7 | 4 L1 outcomes sum to primary units | PASS |
| `determinacy-levels` | plotted unit measures / non-capped complete primary units | 7 | exact-tree unique ≤ shape unique ≤ denominator | PASS；exact table另列 region/W_primary，分母有顯示。 |
| `solver-cap-rate` | primary reconstruction units | 7 | capped / primary | PASS |
| `sex-chromosome-share` | PASS sSNV universe | 7 | `(chrX+chrY+other)/U` | PASS；僅 census，不作 topology。 |

Chart contract 另有：14 unique chart IDs、14 Recharts hosts、14 fallbacks、14 source controls、0 duplicate HTML IDs。`source_notes.json` 14/14 IDs 與 HTML 一致。

## 6. Table item-by-item coverage

HTML 共 52 tables，52/52 有 caption、header scope；其中 14 chart fallbacks + 14 chart exact tables 已逐 cell 比較。其餘 table groups 如下：

| Table group | Scope | Audit comparison | Result |
|---|---|---|---|
| Aggregate S/U→W→primary | 14 rows | 與上述 6 組 conservation identities 逐項比較 | PASS |
| Formal-state interpretation | 5 rows | `C`、`Topo`、hidden/incomplete 的 semantic definitions 與 generator logic 比較 | PASS |
| Evidence-track status | 4 rows | Frozen 07:31 source 正確；current producer 已不同 | **PASS at cutoff；DQ-01** |
| Chart exact/fallback tables | 28 tables | 359 rows／3,472 exact cells；14/14 pair equality | PASS |
| Hidden × topology × HP breakdown | exact breakdown | 每個 parent subtotal 與 topology / hidden totals 比較 | PASS |
| Seven dataset dossiers | 7 × funnel table + 7 × interpretation table | companion sample objects、raw MLHP/region/layered JSON 比較 | PASS；HCC gap 保留而未掩蓋 |
| Complexity / L2 CN | 7 rows each | primary units、candidate totals、recurrence=131 及 per-sample CN split 比較 | PASS；CN unavailable 沒有被填成 0 evidence |
| HCC1395 drill-down | funnel、region bridge、13-layer metrics | raw HCC files與精確分母比較 | PASS within available region view；DQ-04 |
| Unit dictionary | 12 rows | site、region、primary unit、pattern observation 不混用 | PASS |

## 7. Five-dimension data-quality assessment

| Dimension | Evidence | Verdict |
|---|---|---|
| Completeness | 7/7 historical sample dirs、35/35 MLHP parts、7/7 region views、7/7 layered files存在；但 read-AF 0/7、L3 0/7、HCC少 1 region。 | **Partial / release NO-GO** |
| Uniqueness | 48,960/48,960 MLHP region keys unique；48,959/48,959 emitted-region keys unique；110,334/110,334 layered `(region,family)` keys unique；7 dataset keys unique。 | PASS |
| Validity | schema/gate 7/7 pass；counts nonnegative；shares in range；impossible topology state 0；method params一致。 | PASS within historical engineering scope |
| Consistency | Raw JSON→companion→headline/table/chart 無 numeric mismatch；site、region、unit funnels全部守恆。 | PASS；current production track是唯一 freshness mismatch |
| Integrity / provenance | Root hash 2/2、sample output manifests 7/7（49/49 files）PASS；historical producer executable scripts 0/6 match；6/7 truth-BED；current builder已與 frozen companion hash不同。 | **Partial / release NO-GO** |

## 8. Findings、正確來源與最小修正

### DQ-01 — Production status 的 generic label 已跨 evidence track（P1 High；confidence High）

- Frozen HTML evidence cutoff：`2026-07-11T07:31:04+08:00`。
- HTML/source notes 指向 `20260711_longphase_s_production_sidecars_PASS_ONLY_ABORTED_v1`：`0/7 PASS`、`E_METHOD_SCOPE_PASS_ONLY`；此值對其 source/cutoff 正確。
- Audit `2026-07-11T15:35:27+08:00` 的正確 current source 是 `20260711_longphase_s_raw_all_production_sidecars_v2`：HCC1395、HCC1395_DORADO、COLO829 已 PASS，H1437 START，故 latest sample state 為 `3 PASS + 1 START`；aggregate `_SUCCESS` 仍不存在。
- 影響：讀者可能把舊 run 的 `0/7` 當 current producer 狀態；但 NO-GO 結論仍正確，因 current run 尚未 7/7 + aggregate gate。
- 正確來源：current raw-all root 的 `run_status.tsv`、sample `_SUCCESS`、run-level `_SUCCESS`／`verification_summary.json`。
- 最小修正：保留兩張獨立 status cards，名稱直接帶 run ID：`PASS-only aborted v1: 0/7`、`normalized raw-all v2: N/7 as-of <timestamp>`；禁止用一個 generic `Production probe PASS` 覆蓋兩條 track。Active status 必須附 evidence timestamp，不可從 static companion 假裝 live。

### DQ-02 — Historical source 不具 comprehensive genome-wide eligibility（P1 High；confidence High）

- `VALIDATION_SCOPE.md` 與 contract audit 證明 6/7 historical LongPhase-S commands 使用 `--truth-bed`；唯一例外為 COLO829。Genome-wide paired ClairS PASS VCF 與 truth-BED-conditioned HP/PS tagging scope 不一致。
- HTML 已正確標示 historical engineering snapshot、6/7 truth-BED、publication NO-GO；所以這是 source eligibility blocker，不是隱藏的 HTML mismatch。
- 影響：歷史 rates/trees 不能 promotion 為 genome-wide read-tag evidence，不能作外部 handoff。
- 正確來源：`/big7_disk/.../20260710_232501_layered_reconstruction_v2/VALIDATION_SCOPE.md`、`InterSubMod/research/20260710_layered_reconstruction_v2/clairs_longphase_ssnv_contract_audit.json`。
- 最小修正：只在 normalized raw-all producer 7/7 sample PASS、run-level aggregate PASS、clean layered-v3 7/7 identity/conservation PASS 後，整批替換 historical quantitative panels；不可用目前 3/7 partial producer 混補舊報告。

### DQ-03 — Historical producer code 無 exact reproducibility（P1 High；confidence High）

- Output integrity：7/7 sample manifests PASS。
- `code.sha256`：7 entries 中只有 1 entry match；該 entry 是 `layered_v2_input_manifest.json`，不是 executable code。六支 scripts：`sm_linkage_genomewide.py`、`sm_multilocus_combinations.py`、`tree_enumeration_solver.py`、`layered_tree_reconstruction.py`、`build_region_view.py`、`verify_layered_v2.py` 全部 hash mismatch，即 **0/6 executable scripts match**。
- 影響：可以證明 frozen output bytes 未變，但不能用 current scripts exact reproduce historical outputs。
- 正確來源：historical run root 的 `code.sha256`；不要以 sample `output.sha256` 代替 code provenance。
- 最小修正：每次 run 將實際執行的 scripts 複製到 immutable run-scoped `code/`，manifest 雜湊指向該 copy，另記 repo commit / dirty patch。若 exact code 無法復原，historical run 永久保留「engineering snapshot only」，以 clean rerun 取代。

### DQ-04 — HCC1395 region-view 少 1 region／4 sSNV（P2 Medium；confidence High）

- `W_pre=7,931 → W_ret=7,928 → W_tree=7,927`；唯一 missing key：`chr6:104657214-104669173`，`k=4`、`n_full_cov_reads=3`、`n_populations=0`。
- HCC retained sites `25,639`，region-view site sum `25,635`，差 4。
- 相對大小：HCC region `1/7,928=0.01261%`；全 aggregate region `1/48,960=0.00204%`；HCC sites `4/25,639=0.01560%`；全 aggregate retained sites `4/182,400=0.00219%`。
- HTML 與 validation 已揭露此缺口，且 topology denominators 使用實際 emitted `W_tree/W_primary`，沒有把 missing region 偽造成 zero candidate；顯示邏輯正確。
- 正確來源：HCC1395 的 5 個 `mlhp_part_*.json` 與 `layered_region_view_HCC1395.json` key-set difference。
- 最小修正：讓 region-view producer 對每個 W_ret group 明確產出 `emitted` 或 `rejected(reason)` record；釐清 `n_populations=0` 是否應 fail/drop。修復後要求 `MLHP keys = emitted + explicitly rejected`，不可靜默遺失。

### DQ-05 — Read-AF、L3 methylation、orthogonal truth 缺失（P1 High；confidence High）

- candidate-tree read-AF ordering artifacts：0/7。
- L3 methylation：0/7 datasets evaluated、0 units。
- 無 single-cell／multi-region orthogonal truth。
- HTML 正確使用「not generated / not evaluated」，沒有把缺資料填成 biological zero。
- 影響：不能 rank/confirm candidate trees、不能換算 CCF、不能確認 clones/subclones，也不能主張 methylation 支持特定 topology。
- 正確來源：historical run root artifact inventory、7 個 `layered_reconstruction_<sample>.json` 的 `L3_methyl`、validation reasons。
- 最小修正：分三個獨立 evidence contracts 產出 canonical read-AF、bounded L3 methylation、orthogonal validation；每一層明定 grain、denominator、missingness。完成前維持 NO-GO，不得以 0 代替 missing。

### DQ-06 — Per-sample truth-conditioning provenance 不夠 machine-readable（P2 Medium；confidence High）

- Historical `input_manifest.json` 有 7 sample/biological IDs 與 input paths，但沒有每 sample 的 `truth_vcf_conditioned`、`truth_bed_conditioned`、tagging command hash。
- Companion sample objects也沒有這些 per-row provenance fields；`source_notes.json` 只提供 global `6/7`，必須跳到另外的 contract audit 才能知道 COLO829 是唯一例外。
- 影響：自動 QA 能阻止 global promotion，但不能只靠 report companion 做 per-dataset provenance slice。
- 正確來源：`clairs_longphase_ssnv_contract_audit.json` 的 per-sample `tagging_command`、`truth_bed`、`truth_bed_removal_executed`。
- 最小修正：在 future input manifest / companion schema 每 row 加 `tagging_source_run_id`、`tagging_command_sha256`、`truth_vcf_conditioned`、`truth_bed_conditioned`、`scope_role`；global 6/7 必須由 per-row fields 重算，而不是 hard-code。

### DQ-07 — Audit 期間 current generator 與 frozen artifact 產生 in-progress drift（P2 Medium；confidence High）

- Frozen companion 記錄 builder SHA-256：`94fa61a500d28de94874c35c3cc6e182b454290c88f520a81c317e0aae60bbda`。
- Audit 開始時 current builder 仍吻合；`15:38:57+08:00` shared workspace 並行修改後，`15:39:31` current hash 為 `161dced86ff51d1289d21b6d92435ad5949236c7fc0f9c403ddb80a6a648a10c`，但 frozen HTML/data 仍是舊 hashes。
- 這不證明 frozen artifact 當初用了錯誤 builder；它表示「current generator + frozen data」已不是同一 atomic release set。
- 正確來源：companion `report_generation_provenance.builder.sha256` 與 current file `sha256sum`。
- 最小修正：generator 穩定後，以單一 atomic command 重新生成 HTML、data、source-notes、validation，更新全部 hashes，再完整跑 numeric + browser QA。不可手改 frozen HTML 或只更新一個 companion。

## 9. Controls that passed

- Frozen HTML SHA-256：`2098e2493e621843b1e9692053c4396d17af0060a508abe4d52a75f18fa842af`，與 validation 記錄一致。
- Frozen data SHA-256：`78cec57adff20d6d3c21acd29e75499acee3be220fdc27012a0029bcbe9d3042`，與 validation 記錄一致。
- Root `verification.sha256`：2/2 entries PASS。
- Sample `output.sha256`：7/7 manifests、49/49 listed output files PASS。
- MLHP uniqueness：48,960 total / 48,960 unique keys；35/35 parts present。
- Region uniqueness：48,959 total / 48,959 unique emitted keys。
- Layered uniqueness：110,334 total / 110,334 unique `(region,family)` keys。
- Raw→companion field mismatch：0。
- Chart/table mismatch：0。
- HTML source controls：14/14 charts；52/52 table captions；duplicate IDs 0。

## 10. Reproducible commands與實際輸出

### 10.1 Raw / companion / HTML 重算

輸入：frozen HTML/companion JSON、historical run root 內 35 MLHP + 7 region-view + 7 layered JSON。

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python /tmp/audit_layered_report.py > /tmp/audit_layered_report_results.json
jq '{chart_rows_compared,chart_cells_compared,chart_mismatches,raw_sample_mismatches,aggregate,raw_uniqueness}' \
  /tmp/audit_layered_report_results.json
```

實際輸出片段：

```text
chart_rows_compared: 359
chart_cells_compared: 3472
chart_mismatches: []
raw_sample_mismatches: []
aggregate: U=568080, retained=182400, W_ret=48960,
           W_tree=48959, W_primary=47377, primary_units=68544
```

### 10.2 Frozen artifact與 current builder hashes

```bash
sha256sum \
  research/20260710_layered_reconstruction_v2/scripts/build_layered_observation_report.py \
  docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.standalone.html \
  docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.data.json
```

`2026-07-11T15:39:31+08:00` 實際輸出：

```text
161dced86ff51d1289d21b6d92435ad5949236c7fc0f9c403ddb80a6a648a10c  .../build_layered_observation_report.py
2098e2493e621843b1e9692053c4396d17af0060a508abe4d52a75f18fa842af  ...standalone.html
78cec57adff20d6d3c21acd29e75499acee3be220fdc27012a0029bcbe9d3042  ...data.json
```

### 10.3 Current production status

```bash
date --iso-8601=seconds
tail -n 20 /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/run_status.tsv
find /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples \
  -mindepth 2 -maxdepth 2 -name _SUCCESS -printf '%h\n' | sed 's#.*/##' | sort
test -f /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/_SUCCESS \
  && echo AGGREGATE_SUCCESS=1 || echo AGGREGATE_SUCCESS=0
```

`2026-07-11T15:35:27+08:00` 實際輸出片段：

```text
HCC1395         PASS
HCC1395_DORADO  PASS
COLO829         PASS
H1437           START
AGGREGATE_SUCCESS=0
```

### 10.4 Audit output

本次唯一 project write：

```text
InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/audit/20260711_data_quality_findings.md
```

未修改：target HTML、generator、companion data/source-notes/validation、historical outputs、current producer outputs。

## 11. Release gate recommendation

維持 `artifact_validation=PASS` 與 `scientific_release=NO_GO` 的雙狀態。下一版只有在以下條件全部成立時才可重新審核 promotion：

1. normalized raw-all producer 7/7 sample PASS 且 aggregate `_SUCCESS` / verification all-pass；
2. clean layered-v3 7/7 identities、group/region/site conservation 全通過；
3. run-scoped executable code snapshot 100% hash match；
4. HCC1395 missing region 已修復或有 explicit rejected-record accounting；
5. read-AF、L3 methylation、orthogonal truth 的 missingness 與 claim ceiling重新審核；
6. 重新生成全部 report artifacts，numeric audit 與 desktop/mobile browser QA 皆為 0 blocking errors。

在此之前，現有 HTML 適合作「方法與 historical engineering distribution 展示」，不適合作 current production dashboard 或 biological confirmation。

## 12. 2026-07-11 修正後 data-contract readback（16:24 Asia/Taipei）

### 12.1 根因與修正邊界

這次問題不是 frozen historical 數值算錯，而是三種不同資訊被壓在同一閱讀層：

1. 已中止的 PASS-only evidence track 與 active normalized raw-all track 沒有先後層級；
2. 每個數字旁的 `.json` source tooltip 讓 provenance 標記比數值本身更搶眼；
3. historical 7-row aggregate 沒有在首屏先提醒它只代表 6 個 biological samples，且不是 current canonical output。

修正後保留 historical 數值作可重算 baseline，但將 current status 放到首屏，並把所有 JSON 來源改成摺疊區或人類可讀「來源」連結。這是展示與語意修正，不提高科學證據層級。

### 12.2 可重算 regression test

輸入：

- `InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.data.json`
- `InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.validation.json`
- `InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/20260711_分層重建全景數據觀察_01.standalone.html`

命令：

```bash
python3 docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/audit/validate_panorama_data.py
```

輸出：

```text
InterSubMod/docs/reports/in_progress/2026/07/20260711_分層重建數據全景觀察_01/audit/after/data_quality.json
```

實際輸出片段：

```text
status=pass
checks_total=70
checks_passed=70
checks_failed=0
aggregate U=568080, retained=182400, W_tree=48959,
          W_primary=47377, primary_units=68544,
          complete=39885, incomplete=7492
production state=in_progress, pass=3/7, aggregate_complete=false
```

測試逐 dataset 驗證 site funnel、region bridge、HP partition、L1、topology、hidden、HP×H3、k；另驗證 7 rows／6 biological samples、HCC1395 technical replicate、aggregate complete/incomplete、normalized raw-all 角色與 gate、HTML/data hash、inline tooltip=`0`、JSON link label 不暴露 `.json`，以及所有 5 個 JSON href target 皆存在。

### 12.3 修正後 verdict

- **Artifact validation：PASS** — 70/70 data checks、37/37 Chromium checks，errors=`0`。
- **Production freshness：PARTIAL** — 2026-07-11 16:23 仍為 `3/7 PASS`，`H1437 START`，aggregate `_SUCCESS` 不存在。
- **Scientific release：NO-GO** — 6/7 truth-BED conditioning、historical executable hash mismatch、HCC1395 1 region／4 sSNV gap、read-AF/L3/orthogonal truth 缺失仍未被這次 UI hotfix 解決。

目前 evidence hashes：HTML=`069980a2de9c236aa075c6a17ab1150ce0df16a17856a3f39a411c89b5e9ffe6`；data=`2a22585776e23b66af1e80169540bb50df28e33382255e59132fd33d808cc4d8`；data-quality=`0a9bb8ad46b85608ec5928e30751880826f898ca2d4c68f0359e34f5391cb70b`。
