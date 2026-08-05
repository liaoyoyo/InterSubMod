<!--
建立時間: 2026-07-16 10:00 +08:00
目標: 獨立 red-team 教授版 PARTIAL HTML 與 fail-closed builder 的數據、分母、名詞及發布 gate
處理範圍: Task Type B；chr1-22；7 technical datasets / 6 biological samples；canonical v5 + M0 full receipt/TSV；不執行 BAM
服務目標: G3 / G4 / G5
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_read_linked_hypercube_exact_likelihood全量驗證報告_01.partial-preview.html
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_validated_html_report.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/m0_receipt.json
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_audit_full_v4/independent_verification.json
  - InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json
狀態: independent_red_team_complete; FINAL promotion needs revision
-->

# HTML 教授報告數據語意獨立 red-team 稽核

## 結論先行

**Verdict：目前 PARTIAL artifact 可安全作內部方法預覽，但 builder 尚不夠資格把未來輸出自動升級為 FINAL。** Canonical 與 M0 的既有數字經 TSV/JSON 獨立重算後一致；主要發布風險不是這兩層的算術，而是：完整漏斗未展示、M2 主 stratum 可因 dictionary insertion order 任選 HP1/HP2、T/V 等非母體比值被格式化成百分比、per-dataset M2 結果沒有被實質驗證，以及 FINAL gate 沒有要求 M0 independent verifier。

- 影響：**高**（若直接升 FINAL，教授可能把不同 grain 的數量當成同一漏斗或把單一 HP stratum 當全體）。
- 信心：**高**（數字由 canonical JSON、M0 TSV 與獨立 receipt 重算；gate 缺口以 adversarial fixture 實測）。
- 本稽核沒有修改共享 builder，也沒有執行 BAM。

## 1. 稽核輸入與內容身分

| 證據 | Scope | SHA-256 |
|---|---|---|
| `InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json` | current canonical；7 datasets / 6 biological samples / chr1-22 | `71da78b69f8afe5fb8e618179ab7b38a6940fdb17be6282f99a6ec4b720e5de7` |
| `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/m0_receipt.json` | M0 full 7-dataset receipt | `eba081a70f16c008f70cd97c85ec3bcbce41d3982eb6a55c5915e89149197699` |
| `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/m0_hp_lineage_likelihood_census.tsv.gz` | 64,973 eligible HP-unit rows | `9df74db30bc930fee4e0a6941f371bfe12b069808423bb53be5eaf3fc77c1a6c` |
| `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_audit_full_v4/independent_verification.json` | deep verifier；64,973 / 64,973 units | `912760104cb3b7dca4e18d56ac429f0cbc901c81800719425ebaf228c2ae735c` |
| `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/raw_all_receipt_closeout.json` | raw ClairS → LongPhase-S producer continuity | `d5360b7422ee1017d71f60ede51bc3283160e57c6abbb78daa3e74090362ab7c` |
| `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_validated_html_report.py` | 本次靜態稽核版本 | `e0715d906a548728c00bb462adc28048881743afd7cab0dafb0329157eade319` |
| `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_read_linked_hypercube_exact_likelihood全量驗證報告_01.partial-preview.html` | 本次 rendered PARTIAL | `7275573cd1a2e23baf69dfa9aa8c5cc143437a4e20f21954c0bae18546710e85` |

> Builder 與 HTML 仍在同一 research cycle 迭代；後續修改後必重新計算 hash 並重跑本稽核相關測試。

## 2. 已獨立確認的數字

### 2.1 上游與 legacy site / region 漏斗

| 層 | 數量 | 正確分母與比例 | 判讀 |
|---|---:|---:|---|
| normalized raw ClairS records | 638,259 | 起始母體 | 7 datasets 的 record 加總，不是 unique biological loci |
| LongPhase-S PASS | 582,820 | 91.31% / 638,259 | current tree-input records |
| chr1-22 biallelic sSNV | 469,849 | 80.62% / 582,820 | primary autosomal site scope |
| positional singleton sites | 50,432 | 10.73% / 469,849 | legacy 50 kb grouping 的 k=1 site；不是失敗 read |
| MAX_SNV=8 cap-excluded sites | 225,268 | 47.94% / 469,849 | legacy densest-8 計算 cap 排除；不是生物邊界 |
| read-unsupported sites | 0 | 0.00% / 469,849 | 本 canonical 口徑 |
| retained sites | 194,149 | 41.32% / 469,849 | `50,432 + 225,268 + 0 + 194,149 = 469,849` |
| legacy W_tree | 51,815 regions | 單位由 site 轉 region，不應除以 site | positional multi-site groups |
| W_primary | 50,215 regions | 96.91% / 51,815 | 至少一個 primary HP1/HP2 mutation lineage |
| no-primary-lineage | 1,600 regions | 3.09% / 51,815 | 與 W_primary 守恆 |
| candidate-complete | 42,240 regions | 84.12% / 50,215 | 全 primary units candidate-complete |
| candidate-incomplete | 7,975 regions | 15.88% / 50,215 | 不可直接等同 5,854 個 `MAX_SNV`-capped groups |

重要區分：`n_groups_capped_by_MAX_SNV=5,854` 是 site 選擇階段的 group 數；`candidate-incomplete=7,975` 是 region-level candidate completeness。兩者不是同一欄位、不是同一原因集合，也不得共用「capped」一詞而不加限定。

### 2.2 Current canonical topology

| Region-level joint class | 數量 | 占 W_primary 50,215 | 占 complete 42,240 |
|---|---:|---:|---:|
| joint T=1、Topo=1 | 11,582 | 23.06% | 27.42% |
| joint T>1、Topo=1 | 10,737 | 21.38% | 25.42% |
| joint T>1、Topo>1 | 19,921 | 39.67% | 47.16% |
| candidate-incomplete | 7,975 | 15.88% | 不適用 |
| 合計 | 50,215 | 100.00% | complete 三類合計 100.00% |

因此「Topo>1 很多」的正確數字是：`19,921 / 42,240 = 47.16%`（只看 complete），或 `19,921 / 50,215 = 39.67%`（全部 W_primary）。報告必須同時列出兩種分母；不能只顯示 stacked-chart 的 W_primary 比例。

### 2.3 M0 HP-lineage unit census

Unit 轉換必須先說清楚：`50,215 regions → 72,994 primary HP-lineage units`。M0 再分成 `64,973 eligible + 8,021 capped/excluded = 72,994 units`；64,973 units 分布於 46,385 regions，42,240 regions 是全部 primary units 都 eligible。

| 互斥分類 | 數量 | 占 64,973 eligible units | 相對分母 |
|---|---:|---:|---:|
| T=1、V=1 | 26,225 | 40.36% | 40.36% / all eligible |
| T>1、V=1（edge-only） | 101 | 0.16% | 0.26% / 38,748 T>1 |
| T>1、V>1：likelihood unique | 8,751 | 13.47% | 22.64% / 38,647 V>1 |
| T>1、V>1：likelihood tied | 27,270 | 41.97% | 70.56% / 38,647 V>1 |
| T>1、V>1：optimizer abstain | 2,626 | 4.04% | 6.79% / 38,647 V>1 |
| 合計 | 64,973 | 100.00% | 守恆 |

其他重算：raw parent-edge candidates `T=444,007`；distinct vertex sets `V=443,745`；64,973 rows 觸及 46,385 distinct dataset×region keys。以上與 M0 receipt 及 independent verifier 完全一致。

## 3. 發布 blocker（FINAL 前必修）

### P0-1 — FINAL gate 沒有要求 M0 independent verifier

`_validate_m0()` 只驗 receipt 自報 aggregate / by-dataset 守恆，沒有要求 TSV 存在、TSV SHA、64,973-row 深驗、frozen-solver重建或 independent verifier `verdict=PASS`。現有 final fixture 的 M0 甚至沒有 `outputs` 欄位，仍回傳 `([], [])` 並可生成 FINAL。

**影響**：HTML 可以在沒有 row-level M0 evidence 的情況下宣稱 Task-B gate pass。

**建議 patch contract**：新增必填 `--m0-independent-verification`，驗證 verifier schema/version、`verdict=PASS`、`n_errors=0`、receipt/TSV SHA、deep checked `64,973/64,973`、categorical conservation與 canonical input hashes，並列為新 evidence source。

### P0-2 — M2 headline stratum 是 order-dependent，且可能只展示 HP1

`_primary_ranking_cell()` 取第一個 key 名含 `hp` 的 basis；相同 payload 只改 dictionary insertion order，主圖會在 `PS_HP1` 與 `PS_HP2` 間切換。正式 M2 runner 預期輸出兩個 disjoint primary bases：`PS_HP1`、`PS_HP2`；目前報告不會自動做全 primary aggregate。

**影響**：標題與 scope 寫 7 datasets / 6 samples，但主要 status / topology / partial chart 可能只代表 HP1 或 HP2。

**建議 patch contract**：receipt 明示 `report_primary_basis` 與 `report_primary_threshold`；教授主表至少同時顯示 HP1、HP2，若合併則只對已證明 disjoint 的 unit counts 做守恆加總。不得靠 JSON key order 選結果。

### P0-3 — `T/V`、candidate rows/unit 等非母體比值被 `_metric_cell()` 當百分比

全部 strata 表把 raw T 的「相對分母」設為 V，會顯示 `100×T/V%`；candidate-table metadata 也把 candidate rows 除以 units 顯示百分比。這些不是 population share，可能大於 100%。即使文字括號寫「不是百分率母體」，畫面仍印出 `%`。

**建議呈現**：

- `T`、`V` 各列純 count；另列 `T/V`（trees per vertex set）與 `V/T`（edge-collapse ratio），明示 ratio、不是比例。
- candidate table 列 `rows`, `units`, `mean candidates/unit`，不顯示 `%`。
- W_tree、component total、parent choices 等沒有第二個合理母體時顯示 `relative denominator: N/A`，不要製造 100%。

### P0-4 — 使用者要求的完整漏斗目前沒有被報告

PARTIAL 只從 `582,820` 開始，並跳過 `469,849 → 50,432 singleton + 225,268 cap-excluded + 0 unsupported + 194,149 retained`、`51,815 → 50,215 + 1,600`、`72,994 → 64,973 + 8,021`。這些正是解釋 site、region、HP unit 為何不能直接相除的核心。

**建議 patch**：把 §2.1 與 §2.3 的 cross-grain bridge table 放進教授摘要；raw `638,259` 必須引用 production closeout receipt，不能手抄。

### P0-5 — per-dataset M2 結果既未實質驗證，也未完整展示

`_validate_m2_ranking()` 只要求每個 dataset 的 nested mapping 非空，沒有逐 cell 跑 `_validate_ranking_cell()`，也沒有重算 7 datasets 加總是否等於 aggregate。現有 final fixture 的 per-dataset cells 是空 `{}`，10/10 builder tests 仍 PASS。HTML 只展示各 dataset molecule funnel，沒有各 dataset 的 solver complete、V/T、unique/tie/abstain與 coarse topology。

**影響**：即使 FINAL，仍不能回答「每個樣本最後有哪些 topology、數量與比例」。

**建議 patch / test**：逐 dataset×basis×threshold 驗證完整 schema與守恆；加總逐欄等於 aggregate；教授頁顯示 7 datasets 的 primary threshold status/topology table。HCC1395 與 Dorado 必須並列但標成同一 biological ID，不可當兩個生物重複。

### P0-6 — M2 extraction 的 154-task gate 仍有可穿透缺口

`_validate_m2_extraction()` 只檢 `PASS + REUSED_PASS = 154` 且 `FAIL=0`，沒有要求所有 status 總和恰為154，也沒有要求 154-row unique `dataset×chr` task index。Adversarial fixture 加入 `SKIP=1` 後總 status=155，validator 仍回傳空 error list。

**建議 patch / test**：要求 status categories 為 allowlist、總和=154、task index長度=154、set精確等於7×22且每 row child schema/hash/parameters/output checks PASS。

## 4. 高優先語意修正

1. **`C` 名稱衝突**：current canonical JSON 的 `C_region` 是「跨 primary HP units 的 exact tree-candidate product」；新 M2 想把 `C` 用作 structural R/A/X pattern groups。兩者不可共用同一符號。建議新報告完全使用 `n_structural_pattern_groups`，並加 migration note：legacy `C_region` 已改名 `joint_T_region`，不與 read-group C 混用。
2. **`T` 必須帶 grain**：canonical T 是 joint region candidate count，M0 T 是 per-HP unit，M2 T 是 component×HP×PS unit。定義應寫成 `T(unit)`；raw totals 不能跨 grain 直接比較。
3. **candidate-incomplete ≠ MAX_SNV-capped**：7,975 與 5,854 必須分開；前者為 region completeness，後者為 legacy site-cap group count。
4. **Topo 需雙分母**：19,921 同時列 `39.67% / W_primary` 與 `47.16% / complete`；補一段「高歧義來自 partial evidence、structural symmetry與 snapshot edge non-identifiability，不代表存在很多真 clone」。
5. **M0 per-dataset relative denominator**：aggregate unique/tie/abstain 正確用 V>1=38,647；per-dataset table卻用 T>1（包含 edge-only）。最大偏差在 HCC1937：unique 顯示44.66%，正確 V>1-relative 是45.31%。應新增 `ds_v_gt1=unique+tied+abstain`。
6. **S8 pilot標示過強**：現有 HCC1954 chr22 receipt 是 schema `1.0.0`，不是 PS-aware 1.2。來源名稱與區塊需明示 `LEGACY DIAGNOSTIC / single chromosome / not PS-aware`；不要只稱「lossless M2 pilot」。HP表分母應叫「eligible alignments by raw HP value（含 missing `.`）」而不是 `HP-tagged eligible alignments`。
7. **symbolic 高效敘述需限縮**：group constraint 的「表示大小」隨 pattern groups線性增加，但求全域 minimum vertex set 仍是 combinatorial / NP-hard；不可讓「近似線性」看似 runtime claim。
8. **2ᵘ 必須先講 universe**：只有 full k-cube 是2ᵘ；observed-ALT structural universe 是 `2^popcount(free_mask & universe)`，且 fixed ALT 不在 universe 時可為0。標題應寫「full cube 概念上2ᵘ」。
9. **candidate table gate**：目前沒有重算 semantic SHA、沒有驗 candidate_id order/unique、unit_key與欄位一致性、`candidate_set_complete` boolean、winner/tie/relative-LL unit-level守恆。FINAL 前至少需由獨立 ranking verifier補足。
10. **canonical per-dataset gate**：目前只驗各 count 跨 dataset 加總；把一個 dataset complete +1、另一個 -1，造成個別 `complete+incomplete != W_primary` 仍可通過。應逐 dataset檢 complete partition與 topology partition。

## 5. 已通過且應保留的設計

- **PARTIAL fail-safe 有效**：缺 full M2 時，不建立 final-looking檔；CLI `generation_pass=true` 但 `all_pass=false`、`final_ready=false`，HTML固定顯示 `PARTIAL PREVIEW · NOT VALIDATION EVIDENCE`。
- **legacy M2 schema 不會升 FINAL**：extraction 1.1 / ranking 1.0 只能 diagnostic；final 要求1.2 / 2.0。
- **dataset / biological sample 有分開**：HCC1395與HCC1395_DORADO明示為兩個 technical datasets、同一 biological sample。
- **chr1-22 / chrX 理由合理**：性別、non-PAR hemizygosity、PAR、X-inactivation與腫瘤CN/LOH使chrX需要獨立copy-state/HP期望；不是因ONT不能讀chrX。報告也正確提醒autosome仍需CN sensitivity。
- **partial read核心語意正確**：每條pattern是一個group constraint，所有groups聯合求最小extra/hidden；不是逐completion「任一成功」。complete exact units保留全部全域同分vertex sets；incomplete/capped必須abstain。
- **VAF不重複計分合理**：joint molecule patterns已包含同reads的marginal allele-frequency資訊；π是固定vertex set下的latent expected proportions。再用同reads算VAF加第二個分數會重複使用證據。獨立caller/CN/甲基若另有預先定義模型，可作外部sensitivity，但不能默認成第二次read score。
- **claim ceiling 合理**：結果只稱 regional mutation-state candidates，不稱真cell clone數、唯一parent edge或已確認完整演化樹。
- Browser QA receipt（partial v5）為PASS：11 sections、5 SVG、9 details、desktop/mobile無overflow、無external/console/page error。

## 6. Partial-read 實際方法的最終判讀

### 結構層

對 `p∈{R,A,X}^k`，若有 `u` 個X，full cube中有 `2ᵘ`個概念completions；實作不逐一跑樹，而是保存：

```text
G(p) = {v : ((v XOR alt_mask) AND fixed_mask) = 0}
constraint: N ∩ G(p) ≠ ∅
```

所有full observations、root及partial groups一起進一次global objective；依extra vertex數由小到大找可行解。只有完整枚舉完最小層時，才能把全部同分minimum `N`輸出為complete candidate set。Frozen v5若遇`extra_cap/per_level_budget`會greedy fallback並標`capped`，不能冒充all-optimal。

因此使用者句子應修正為：

> 「partial read有多個相容state，但不為每個completion各建一棵樹；把每條partial read編成一個subcube group，聯合所有read求全域最少hidden/extra的candidate vertex sets。只在solver有completeness certificate時保留全部同分最小候選；無certificate則abstain，再由read-pattern likelihood對可比較的vertex sets排序。」

### likelihood層

未知X維度的emission factor為1（marginalized missing）；固定R/A依BQ conditional symmetric-substitution model計分，再對latent vertex state加總：

```text
P(pattern | V, π) = Σ_v π_v × P(fixed R/A cells | v, BQ)
```

這不是「只要任一completion像就算整條read成功」，也不是把每個completion當獨立read。相同vertex set的不同parent edges有相同snapshot likelihood，所以edges保留不可辨識。

### 更高效但仍精確的方向

1. bitmask symbolic group（已採用），避免顯式`2^Σu` joint completions；
2. active-bit compression、duplicate/dominance group reduction與mandatory-hit preprocessing；
3. sparse/lazy group rows，不建立PS×whole-chromosome dense矩陣；
4. k較大時使用有certificate的MILP/BDD/ZDD/branch-and-price或separator boundary-state DP；
5. 無global certificate則輸出local/abstain，不能把任意切8點稱global optimum。

## 7. Step → Verify 修正清單

1. 加入完整S/W/HP-unit crosswalk → 驗證：上述每一層守恆，且site→region轉換標為「不可相除」。
2. 移除order-dependent primary basis → 驗證：交換JSON中PS_HP1/PS_HP2順序後HTML headline及數字byte-identical。
3. 修T/V與candidate rows呈現 → 驗證：無任何非population ratio帶`%`；`T/V`明示ratio。
4. 加M0 independent receipt gate → 驗證：缺verifier、錯receipt SHA、少1 TSV row或deep count<64,973均exit 2且不產FINAL。
5. 強化extraction 154-task gate → 驗證：extra SKIP、duplicate task或缺chr22均exit 2。
6. 驗證及展示per-dataset M2 cells → 驗證：任一dataset空cell、七樣本加總與aggregate差1均exit 2；HTML有7-row status/topology table。
7. 修C/T/Topo名詞 → 驗證：全文不再把legacy `C_region`與`n_structural_pattern_groups`稱同一C；T每次都有unit grain。
8. 重跑HTML/browser/print QA → 驗證：FINAL ribbon、所有local sources存在、desktop/mobile無overflow、print無裁切、console/page error=0。

## 8. 實際執行命令與輸出片段

### 8.1 Builder tests

輸入：

- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_validated_html_report.py`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/fixtures/validated_html_report_bundle.json`

命令：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 -m unittest \
  research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/test_build_validated_html_report.py -v
```

輸出：terminal stdout；未建立持久分析輸出。

實際片段：

```text
Ran 10 tests in 0.574s
OK
```

解讀：現有測試對已寫的gate全部PASS；但fixture本身沒有M0 outputs/independent verifier、extraction task index及非空per-dataset ranking cells，故10/10不能證明本稽核指出的缺口已關閉。

### 8.2 Adversarial gate probe

輸入：同上fixture與builder。

命令：以Python只讀import builder，測`_validate_m0`、`_validate_m2_extraction`與`_primary_ranking_cell`。

實際片段：

```text
m0_has_outputs= False
m0_validate= ([], [])
extraction_status_sum= 155
extraction_errors_with_extra_SKIP= []
primary_basis_order_A= ('PS_HP2', '3')
primary_basis_order_B= ('PS_HP1', '3')
fixture_extraction_has_task_index= False
```

輸出：terminal stdout；未建立持久分析輸出。

### 8.3 Legacy site funnel重算

輸入：

`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/samples/*/layered_reconstruction_*.json`

命令：逐dataset以`jq`讀`input_funnel`，再由`awk`加總。

實際片段：

```text
TOTAL  469849  50432  225268  0  194149  469849  51815  5854  0  51815
```

欄位依序：scope、singleton、cap-excluded、read-unsupported sites、retained、accounted、precap groups、capped groups、unsupported groups、retained groups。輸出：terminal stdout；未建立持久分析輸出。

## 9. 分享狀態

- **目前 PARTIAL：Share internally with caveats。** Ribbon與footer已清楚阻止把空M2當0或當validation evidence。
- **未來 FINAL：Needs revision。** P0-1至P0-6至少關閉後，再做完整數據、browser與多agent終審。

