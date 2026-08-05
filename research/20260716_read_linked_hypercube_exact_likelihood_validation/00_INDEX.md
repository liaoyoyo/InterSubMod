<!--
建立時間: 2026-07-16 06:50 +08:00
目標: 導航 read-linked Hypercube Steiner 全量驗證的問題、方法、資料層、執行順序與證據
處理範圍: chr1-22、7 datasets / 6 biological samples、LongPhase-S recalibrated FILTER=PASS
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/pre-decision-audit.md
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/implementation-notes.md
-->

# Read-linked Hypercube Steiner exact＋likelihood validation

> **一句話目標**：不再用50 kb與densest-8當作真實read linkage；先由每條ONT molecule建立可稽核的R/A/O/D/S/L/X sparse pattern，再以symbolic group-Steiner求最少hidden state候選，最後只在可辨識的vertex sets間使用error-aware likelihood排序。

## 任務分類與研究目標

- **Task Type**：B Comprehensive validation；最終必須是全chr1–22 × 7 datasets，subset只可作pilot。
- **G3**：read-level mutation／HP evidence。
- **G4**：七dataset與HCC1395跨技術重現性。
- **G5**：保存input hash、solver gap、守恆與abstain原因，能被外部重算。

## 目前正式狀態（2026-07-17 17:56 +08:00）

> **Frozen v2 第一個正式 pilot＝NO-GO。** `HCC1395_DORADO × chr6` extraction receipt PASS；bootstrap=0 ranking於8小時上限被停止（exit 124），缺少完成receipt，獨立single-task verifier fail closed為`NO_GO`。四個ranking gzip為截斷效能證據，不可作候選樹、拓撲或生物比例。B20、H2009 chr2、154-task full、final numeric summary與FINAL topology HTML均未啟動；本次另已產生一份明確標示NO-GO的教授版pilot驗證HTML，不得誤稱FINAL。

- 終端驗證紀錄：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_M2_frozen_v2正式pilot_NO_GO驗證紀錄_01.md`
- 教授版NO-GO HTML：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_M2_frozen_v2_pilot_NO_GO驗證報告_01.html`（official validation／package／desktop＋mobile browser verification PASS）
- 方法教學HTML：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_聯合GroupConstraints到Likelihood排序方法教學_01.html`（Task Type F；解釋joint groups→全部minimum-extra candidates→likelihood；方法測試23/23與browser/print QA PASS；不是全量驗證證據）
- 本次可引用的數字範圍：只限frozen provenance、HCC1395_DORADO chr6 extraction漏斗與NO-GO gate；ranking rows只能標為`diagnostic_only`。
- 下一步：建立新release，先做exact-preserving cache／persistent solver／per-unit deadline／checkpoint，再重新pre-decision audit、freeze與雙pilot；不得resume或清空本失敗root。

## 三個資料層不可混稱

| 層 | 證據 | 可以回答 | 不可宣稱 |
|---|---|---|---|
| M0 | current-v5、MINREAD≥3的alignment-exposure patterns | 舊候選中多少是vertex ambiguity／edge-only ambiguity；工程runtime | lossless、molecule likelihood、true clone |
| M1 | 469,849份既有`reads.tsv` accepted R/A calls | 38,345,639 accepted site-read observations、read-linked component、R/A/X* patterns | X*原因、O/error model、全部raw reads |
| M2 | BAM target alignment＋CIGAR＋LongPhase-S sidecar exact join | R/A/O/D/S/L/X、BQ、molecule、>50 kb bridge、error-aware likelihood | 真實parent edge／cellular clone pairing（仍需獨立lineage evidence） |

## Partial read 的正式語意

partial pattern `p`不是多筆read，而是一個subcube group：

```text
G(p) = {v : ((v XOR alt_mask) AND fixed_mask) == 0}
```

若有`u`個X，完整立方體中有`2^u`個概念相容states。只有count達`structural_exact_pattern_minread`的distinct exact R/A/X patterns會形成structural full/group constraints；低於門檻的informative molecules仍進likelihood。每次MILP build／rebuild會對每個reduced active group列出active universe內的相容vertex indices作為該次sparse hit row係數；不是整個workflow只產生一次，但也不把它們當成`2^u`個獨立tree worlds，更不建立所有reads補值的笛卡兒積。候選state set `N`只需滿足：

```text
N ∩ G(p) != ∅
```

likelihood階段再計算：

```text
P(p | N) = sum_v pi_v P(p | v, error)，v in N
```

因此一筆partial read只計一次；EM／simplex optimizer會隱式分配它可能來自哪些states。

## 固定的求解順序

1. **資料抽取**：每條canonical molecule逐CIGAR呼叫目標sSNV，保留R/A/O/D/S/L/X與BQ。  
   → 驗證：eligible alignment＝exact sidecar joined；molecule/site與call mass守恆。
2. **read-linked boundary**：相鄰site cut的支持為「兩側各至少一個固定R/A call的unique molecules數」。  
   → 驗證：threshold 1／2／3／5全輸出；>50 kb支持另列；不把50 kb當ONT極限。
3. **候選生成**：root＋full states＋每個partial subcube coverage＋monotone unit-edge closure。  
   → 驗證：先最小hidden nodes；k≤8與frozen oracle一致；k9–12 MILP gap=0。
4. **候選等價化**：先依vertex set合併trees。  
   → 驗證：same vertices／different edges分數差=0，標`EDGE_NONIDENTIFIABLE`。
5. **資料排序**：以error-aware molecule-pattern likelihood比較不同vertex sets。  
   → 驗證：conditional binary-flip error grid、conditional-on-fixed-candidate-set ranking bootstrap、held-out與CN/LOH sensitivity；不把同reads VAF再加一次。
6. **報告**：分開數學optimality、資料stability與biological claim ceiling。  
   → 驗證：每個數字有unit、total denominator、relative denominator與abstain原因。

## k上限與更高效策略

- exact gate使用minread-specific `k_observed_alt_active`，不是raw component-site k；完整座標仍保留供provenance。
- `k_observed_alt_active≤8`：frozen exhaustive solver作oracle。
- `k_observed_alt_active=9–12`：symbolic MILP；保存incumbent、dual bound、gap、status。
- `k_observed_alt_active>12`：不再擷取densest-8冒充完整區域；先輸出read-linked完整component，再走overlapping local blocks／separator diagnostics。無法證明全域等價者標`LOCAL_ONLY`或`ABSTAIN_K_GT_EXACT_LIMIT`，不可稱global optimum。
- partial constraint的**列數**隨pattern groups線性增加，不建立跨read的joint-completion Cartesian product；但現行MILP仍會為每個reduced group列出至多`2^u_eff`個相容active-vertex係數，因此row nonzeros仍可能隨`2^u_eff`增加。

## 已完成的驗證

| 驗證 | 結果 |
|---|---:|
| symbolic vs explicit compatible states | 9,840 patterns；2,015,538 checks；0 mismatch |
| frozen solver vs MILP | 66可比較vertex-set checks；0 mismatch |
| k9–12 exact pilot | 12/12 PASS；最大4,096 variables；gap=0 |
| same-vertex/different-edge likelihood | 差`1.42e-14`，數值容差內為0 |
| method-contract migration audit-time focused tests | 149/149 PASS；active HTML builder修正完成後仍須由主流程重跑，不替代154-task真實M2 receipts |
| current research full regression | 162/162 PASS；新增RRA+AAA+RAX、AX+XA與deterministic exhaustive-reduction audit；builder targeted 23/23 PASS；不替代真實M2 receipts |
| current sSNV funnel independent verifier | 322/322 checks PASS；7/7 dataset逐層守恆、path/SHA/ratio均通過 |
| cross-source denominator conservation | 22/22 PASS；funnel＋canonical＋M0的grain、分母與互斥partition一致 |
| M0七資料集深度核對 | 64,973/64,973 eligible units；missing=0、extra=0；verdict PASS |
| M0 robust optimizer pilot | 10/10舊abstain units、690/690 candidate fits取得global KKT certificate；尚不可外推全部2,626 units |
| current-solver identity-bound pilot | solver SHA=`9dbaf8ec...`；receipt SHA=`5cab8c37...`並綁定runner／solver live SHA；9,840 patterns、2,015,538 checks、0 mismatch；k9–12最大runtime=`0.191109856 s` |
| presolve/no-good exhaustive audit | 61,340 cases、1,979,356 predicates、23,909 sparse與21,844 dense no-good pairs；全部0 mismatch；receipt SHA=`a60037e7...` |
| independent profile-likelihood verifier | verifier SHA=`4859598d...`；R/A/X＋BQ＋count＋states＋π獨立重算LL/KKT/relative weight/winner-tie；tamper負向測試PASS；真實154-child coverage待全量執行 |
| QA artifact provenance independent red-team | P0=0、P1=0、P2=3；sealing TOCTOU已修正；QA＋presentation freezer 24/24 PASS；正式流程仍須 same-frozen `verify-only` |
| latest HTML builder final-publication red-team | **GO**；P0=0、P1=0、P2=3；builder 37/37、numeric＋presentation freezer 25/25 PASS；兩張 ledger 各 8 overall＋56 dataset；0/49、7/35、21/49 denominator、S13/S16、no-replace、overwrite、PARTIAL 與 remote guard 全通過 |
| resource-gate科學provenance獨立red-team | **reader-test GO / formal actual-data WAIT**；P0=0、open P1=0、resolved P1=5組、P2=2；checkpoint deep rebuild與ranking child-derived consolidated-candidate lockstep verifier已修補；targeted 117/117、full 327/327、exact-11 11/11全部PASS，仍不代表154+154真實資料已完成 |
| HTML真實瀏覽器QA | v9 PARTIAL preview：19/19 checks PASS；11 sections、5 SVG、11 details；desktop/mobile/print全PASS |
| frozen v2 NO-GO professor HTML official QA | validation／package／verification全部PASS；source dialog與keyboard semantic click PASS；1440＋390 viewports；2 blocks＋1 native funnel chart |
| joint-group方法教學HTML QA | 23/23 targeted method tests PASS；1440／1024／390／320 px無document overflow；7/7 details互動與print展開PASS；外部resource、console/page errors皆0 |
| exact-preserving solver/cache remediation | **PROBE implementation PASS**；prepared base 25→1 build、HiGHS solve仍25；complete-only process cache；341/341 full tests PASS；v2保留歷史NO-GO，v3為current雙pilot runbook |

## PS-aware primary contract（schema 1.2／2.0）

- extraction primary unit＝`HP family × exact known PS × read-linked component × threshold`；missing PS與舊pooled/HP components僅作diagnostic/sensitivity。
- ranking structural minread是exact R/A/X pattern count，grid 1/2/3/5、primary=3；相容partial patterns不自動合併。
- partial `2^u`是full-cube conceptual count；求解不逐一建立completion-wise tree cases，但每次MILP build／每個reduced active group的sparse hit row仍會列出active universe內至多`2^u_eff`個相容vertex indices。
- source-bound method contract同步宣告`active_compatible_vertex_indices_materialized_for_sparse_rows=true`與`completion_wise_tree_worlds_materialized=false`；前者是單一constraint的sparse row係數，後者才是被禁止的逐completion樹世界。
- structural constraints只來自`count >= structural_exact_pattern_minread`的distinct exact R/A/X patterns；全部informative molecule projections（包含低於門檻者）仍進quality-aware likelihood，故minread不是read-discard filter。
- BQ採symmetric substitution並condition on observed∈{R,A}；同reads VAF不重複計分。
- full M2未獲准啟動：frozen v2第一個正式HCC1395_DORADO chr6 B0 ranking於8小時timeout，verifier=`NO_GO`。目前沒有PS-aware全量數字；任何舊pilot或截斷rows均不可當final。
- extractor已改為稀疏PS bridge，ranker已改為exact-PS route index，detail rows已串流；2026-07-18另完成prepared-base＋complete-only process cache。這不是persistent solver：SciPy/HiGHS仍逐候選重建model/presolve；實際wall與solver tail仍須v4雙pilot。
- per-unit solver timing已用`monotonic_ns`分開candidate generation／likelihood／unit total，並保存invoked-only exact nearest-rank p50/p95/p99；真實數值待pilot。
- 2026-07-16 19:30 live audit 仍偵測使用者 all-sSNV root／child 與 40 workers；child `read_bytes` 持續增加，故formal live preflight繼續等待其自然結束。沒有繞過gate、沒有中止使用者作業。

## 主要檔案

- 決策：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/pre-decision-audit.md`
- Living notes：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/implementation-notes.md`
- 正式凍結全量執行與驗證協定：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M2正式凍結全量執行與驗證協定_01.md`
- 歷史 Runbook v2（frozen-v2 NO-GO，不得resume）：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M2正式執行Runbook_02.md`
- 歷史 Runbook v3（prepared-base/cache source；被資源閘門 hotfix 取代）：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260718_M2_exact_preserving修補正式執行Runbook_03.md`
- Current PROBE Runbook v4（修正 all-sSNV cooccurrence audit 假陰性；雙pilot先行）：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260719_M2_exact_preserving資源閘門修補正式執行Runbook_04.md`
- Exact-preserving remediation audit／notes：`InterSubMod/research/20260718_hypercube_exact_preserving_solver_cache_remediation/pre-decision-audit.md`、`InterSubMod/research/20260718_hypercube_exact_preserving_solver_cache_remediation/implementation-notes.md`
- 最終數據彙整器實作與驗證：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_final_numeric_summary實作與驗證紀錄_01.md`
- Presentation code/evidence snapshot：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/freeze_presentation_snapshot.py`
- QA artifact provenance 獨立 red-team：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_QAartifact_provenance獨立redteam_01.md`
- 最新 HTML builder 最終發布與數據綁定獨立 red-team：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_HTMLbuilder最終發布與數據綁定獨立redteam_01.md`
- Resource-gate科學provenance一致性獨立 red-team：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_resource_gate科學provenance一致性獨立redteam_01.md`
- Partial-read 獨立稽核：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/partial-read-method-audit.md`
- Frozen v2正式pilot NO-GO終端驗證：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_M2_frozen_v2正式pilot_NO_GO驗證紀錄_01.md`
- Joint groups→minimum-extra→likelihood方法教學HTML：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260717_聯合GroupConstraints到Likelihood排序方法教學_01.html`
- 方法教學HTML browser/print QA receipt：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/20260717_joint_group_method_html_qa/qa_receipt.json`
- M2 extractor QA：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_M2_extractor程式與合成資料QA_01.md`
- Current raw→retained funnel receipt：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/current_funnel_v1/current_funnel_receipt.json`
- Current funnel independent verification：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/current_funnel_v1/independent_verification.json`
- 最新數字、分母與守恆獨立稽核：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_最新數字分母與守恆獨立稽核_01.md`
- Current-solver identity-bound pilot：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/pilot_v3_identity_bound/pilot_receipt.json`
- Exact solver reductions audit：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_Hypercube_exact_solver縮減與列舉效能稽核_01.md`
- Presolve/no-good exhaustive receipt：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/hypercube_reduction_exhaustive_v1/receipt.json`
- M2 full independent verifier：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/verify_full_m2_receipts.py`
- 原始文獻與主張邊界：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_方法學原始文獻與主張邊界_01.md`
- 教授版HTML（目前明示PARTIAL）：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/20260716_read_linked_hypercube_exact_likelihood全量驗證報告_02.partial-preview.html`
- 瀏覽器QA receipt：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/html_qa_partial_full_m0_v10/browser_qa_receipt.json`
- symbolic／MILP／likelihood：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py`
- M0 census：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/run_m0_likelihood_census.py`
- M2 extractor：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/extract_lossless_read_linkage.py`
- lossless contract：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/lossless_read_contract.py`
- tests：`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/tests/`

## Claim ceiling

最終可以稱：

> solver-certified regional mutation-state candidate sets，及在宣告的molecule/error/CN contract下可被read patterns區分的候選。

不可直接稱：

> 真實clone數、HP1與HP2的cell-level配對、唯一parent edge、完整癌症演化史。
