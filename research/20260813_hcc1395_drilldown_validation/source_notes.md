<!--
建立時間: 2026-08-13
目標: 保存 HCC1395 drilldown 技術報告的 section-to-source mapping、chart contract 與 omission reasons
處理範圍: HCC1395_v1/v3 全量工程稽核、七份 topology dataset 盤點、2026-07-26 legacy methyl browser 方法/crosswalk、direct-generated 視覺 QA、multi-BAM snapshot dashboard 與 Claude Code Round 1-5 execution record
關聯檔案:
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/20260813_HCC1395_drilldown完整驗證與多樣本改進_01.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/results/metrics_audit.json
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round1_review.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round2_review.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round3_review.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round4_review.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round5_review.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/legacy_browser_method_audit.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/legacy_browser_visual_audit.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_dashboard_artifact.json
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/20260813_multi_bam_analysis_overview.standalone.html
-->

# Source Notes：HCC1395 drilldown validation report

## 1. Reporting job

- **Question**: `/bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1/` 是否可提升，哪些指標、驗證、視覺與多樣本資料整理需調整？
- **Primary audience**: technical；主價值是 provenance、metric contract、統計有效性與跨樣本隔離驗證。
- **Decision supported**: 是否可引用 v1/v3 作研究或對外證據，以及下一輪應先修哪些 P0/P1/P2。
- **Scope**: HCC1395 v1/v3 bundle 全量盤點；22 條常染色體 LCA receipt 對照；7 份 topology datasets；legacy standalone 的 30,077 loci、7,415-region funnel、472 candidates、14 displayed cases；direct-generated source-complete light QA；7-dataset/6-biological-sample portable snapshot dashboard。等價多樣本 ISM/lineage bundles 與 full panels/IGV build 未生成，故科學/asset 結論仍為 **PARTIAL / NOT EVALUATED**。
- **Comparison basis**: v1（±1 kb ISM）對 v3（±5 kb ISM）只做描述性版本比較，不視為受控 A/B。
- **Claim ceiling**: descriptive observation and internal data-product QA; not truth-set validation。
- **Services**: G3（read-level epigenetic）、G4（多樣本一致性/reproducibility）、G5（外部可驗證工程品質）。

## 2. Technical-report required structure mapping

| Required role | Visible report section | Mapping note |
|---|---|---|
| Title | `HCC1395 drilldown 完整稽核與多樣本改進` | Direct title |
| Technical summary | `技術摘要：工程可用，但科學引用仍被擋住` | Answer-first；同時宣告 BLOCKED/PARTIAL |
| Key findings with visual evidence | `v1/v3 的可用邊界`、`甲基化缺失不是隨機`、`七份 topology 可比但異質`、`2026-07-26 standalone`、`視覺與權威狀態` | 每張圖前後有解讀、分母與限制 |
| Scope, data, metric definitions | `範圍、資料與指標口徑` | 定義 region、tree coverage、unique rate、linkage、BH FDR |
| Methodology | `方法與可重現命令` | 列完整輸入、命令、輸出與實際片段 |
| Limitations / robustness | `限制、反證與穩健性` | truth 缺口、LCA gate、cluster double-dipping、版本非受控 |
| Recommended next steps | `P0 / P1 / P2 改進清單` | 分離已完成 source hardening 與仍需重建項目 |
| Further questions | `下一輪需決定的問題` | truth scope、七樣本計算預算、canonical biological replicate policy |

## 3. Section → source map

| Report section / claim | Canonical source | Fields / rows used | Evidence class |
|---|---|---|---|
| v1/v3 BLOCKED、selfcheck counts | `InterSubMod/research/20260813_hcc1395_drilldown_validation/results/bundle_overview.csv` | `validation_status_recomputed`, `selfcheck_*` | L1 direct output |
| Region/tree/unique topology | same | `regions`, `regions_with_tree`, `tree_coverage_pct_all_regions`, `unique_best_tree`, `unique_pct_among_tree` | L1 direct output |
| ISM linkage / version delta | same | `methyl_summary_rows`, `methyl_topology_linked`, `methyl_topology_linkage_pct` | L1 direct output |
| Panel/IGV volume and receipt undercount | `results/asset_inventory.csv` + `results/bundle_overview.csv` | `png_*`, `igv_js`, `panel_unreported_bytes` | L1 direct output |
| Input hash/missing source | `results/input_integrity.csv` | topology/MLHP `HASH_MATCH`; v1 ISM `MISSING`; v3 ISM directory `UNVERIFIABLE_DIRECTORY` | L1 direct output |
| ISM availability by active k | `results/methylation_coverage_by_k.csv` | numerator, denominator, percent, `k_bin` | L1 derived by saved audit code |
| Axis raw/BH rates | `results/methylation_axis_metrics.csv` | topology-linked scope; tested/raw/BH counts and gates | L1 derived by saved audit code |
| LCA descriptive delta / causal gate | `results/lca_ab_by_chrom.csv` + `results/metrics_audit.json` | `same_in_bam`, `same_threads`, `lca_ab_summary` | L1 direct/aggregated output |
| Seven topology datasets | `results/cohort_topology_metrics.csv` | per-dataset tree/unique rates, schema/model/AF basis, technical replicate flag | L1 direct/derived output |
| Browser before/after | `figures/04_mobile_before_sticky_overlap.png`, `figures/01_mobile_fixed.png` | 390×844 smoke; patched-render fixture | L2 browser fixture |
| Original status contradiction | `figures/05_v1_selfcheck_before.png`, v1 receipt/selfcheck audit | static “no issue” vs FAIL/SKIP | L1/L2 |
| Independent adversarial review | `InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round1_review.md`、`InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round2_review.md` | provenance, C8/C10/C12, contamination challenge, multi-sample contract | L2 independent code/output review |
| Independent post-fix verification | `InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round3_review.md` | 37 tests; source fixes PASS; remaining technical debt | L2 independent post-fix review |
| Legacy method / selection audit | `InterSubMod/research/20260813_hcc1395_drilldown_validation/legacy_browser_method_audit.md` + `results/legacy_browser_{crosswalk,funnel}.*` | B purity、raw/powered linkage、FP inclusion、18→14 selection gap、coordinate overlap | L1 read-only reconstruction |
| Direct-generated browser QA | `results/legacy_browser_visual_audit.generated.json` | desktop/mobile dimensions、errors、tables、denominators、fake cells、authority overlap | L1 machine receipt |
| Independent final verification | `InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round4_review.md` | 43 tests + five JS syntax checks; method and runtime findings accepted; scientific gate remains BLOCKED | L2 independent final review |
| Multi-BAM metric / layout contract | `InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_dashboard_metric_contract.md` + `InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_dashboard_design_spec.md` | hero/diagnostic/guardrail/detail formulas、sample filter、null policy、missing BAM QC fields | L2 reviewed design contract |
| Multi-BAM artifact payload | `multi_bam_dashboard_artifact.json` | 18 snapshot datasets、10 cards、6 charts、10 tables、2 required access issues、19 sources、macro P25/P75、bounded BAM/tag rails、four-layer/opportunity anchors | L1 generated snapshot |
| Multi-BAM browser interaction | `results/multi_bam_dashboard_browser_qa.json` | All + 7 dataset states、non-borrowing、1440/1024/512/390/320、fold、keyboard/focus/live region、HTML/artifact/PNG hashes、0 browser/external-request errors | L1 machine receipt |
| Multi-BAM independent re-audit | `InterSubMod/research/20260813_hcc1395_drilldown_validation/independent_dashboard_reaudit.md` | fresh builder/payload/hash/8-state/non-borrowing/responsive/accessibility recheck；ACCEPT bounded snapshot only | L1 independent reviewer receipt |
| Claude Code dashboard attempt | `InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round5_review.md` | CLI attempts stopped before Read/Bash by session quota | execution receipt only; **not a review verdict** |
| Design decisions / deviation | `InterSubMod/research/20260813_hcc1395_drilldown_validation/implementation-notes.md` | immutable v1, observation-only ceiling, fail-closed, no seven-sample rebuild | Governance evidence |

## 4. Chart map

| Figure | Report segment | Analytical question | Family / encoding | Source | Supported claim | QA note |
|---|---|---|---|---|---|---|
| `02_cohort_topology_overview.png` | 七份 topology 可比但異質 | 相同 schema/model/AF 口徑下，各 dataset 的 tree coverage 與 unique rate 差多少？ | Two-panel horizontal bars；y=dataset；x=percent；gray=technical replicate | `results/cohort_topology_metrics.csv` | topology 可做逐樣本描述，但異質性大，不能直接 pool | HCC1395_DORADO 明示 † technical replicate；biological n=6 |
| `03_methylation_coverage_by_k.png` | 甲基化缺失不是隨機 | ISM linkage 是否隨 active k 降低？ | Two-series line；x=k bin；y=linked/topology sSNV %；color=v1/v3 | `results/methylation_coverage_by_k.csv` | k=8+ coverage 約 41%，需把 missingness 當選擇偏誤 | v1/v3 window 不同；不可讀成效果差 |
| `04_mobile_before_sticky_overlap.png` | 視覺與權威狀態 | 原始 390px UI 是否有 sticky overlap？ | Browser screenshot, before | Original v1 render | mobile layout defect exists | Screenshot is QA evidence, not scientific evidence |
| `01_mobile_fixed.png` | 歷史 fixture | 第一輪 source patch 後 mobile contract 是否改善？ | Browser screenshot, overlay fixture | Patched assets using v1 payload | 初步修正 overlap/body overflow | 歷史 PARTIAL；不是最終驗收 |
| `05_v1_selfcheck_before.png` | 視覺與權威狀態 | UI authority copy 是否反映 selfcheck？ | Browser screenshot | Original v1 selfcheck view | static success language contradicts 1 FAIL/1 SKIP | Original immutable evidence |
| `06_legacy_standalone_desktop.png` / `07_legacy_standalone_mobile.png` | Legacy 參考 | 舊頁的閱讀節奏與 responsive 現況為何？ | Browser screenshots | Original standalone HTML | progressive disclosure 可參考；mobile 有水平 overflow | 方法 claim 另由 audit 驗證，不由截圖決定 |
| `13_current_generated_desktop.png` / `14_current_generated_mobile.png` | 最終視覺驗收 | source-complete generated shell 是否呈現 claim ceiling、漏斗與 BLOCKED authority？ | Direct-generated browser screenshots | Full HCC v3 summary; `igv=skip, panels=0` | 1440/1440、390/390，無 overlap/errors | Canonical runtime QA；仍是 light asset policy |
| `15_current_generated_cooccur_desktop.png` | 分母 QA | co-occurrence tables 是否有真實分母且無假互動？ | Direct-generated screenshot + machine counters | Generated page | 4 tables、106 denominator labels、0 fake cells | UI contract，不是生物 validation |
| `16_current_generated_selfcheck_mobile.png` | 權威 QA | mobile selfcheck 是否如實顯示失敗？ | Direct-generated screenshot | Generated page | 2 FAIL / 0 SKIP → BLOCKED | 保留失敗是正確行為 |
| `17_current_generated_detail_desktop.png` | Progressive disclosure | 位點 detail 是否和 observation scope 分離？ | Direct-generated screenshot | Generated page | detail availability 與統計 evidence 分開呈現 | panels=0；不驗圖像像素內容 |
| `18_current_generated_methyl_detail_desktop.png` | 甲基 detail QA | 有 ISM summary 的位點能否同時呈現 read/CpG、raw p、effect、第三態與 circularity gate？ | Direct-generated screenshot + machine-selected axis-code locus | Generated page `chr1:1,320,793` | 105 reads、371 CpG；ALT/REF raw p=.001/effect=.029；HP-fine/cluster/lineage 不被塗成陰性 | 探索性單位點；不是 cohort 或因果證據；不含 raster panels |
| `19_multi_bam_dashboard_all_desktop.png` | Multi-BAM 總覽 | All scope 能否先呈現 blocker、scope、technical/biological gates 與 macro 指標？ | Portable dashboard screenshot | artifact `scope_summary` + `availability` | 7 datasets / 6 biological；7/7 technical、0/7 family/objective、77.7%/62.4%、1/7 downstream | snapshot `partial`；不是 truth validation |
| `20_multi_bam_dashboard_hcc1395_selected.png` | Sample filter | 選擇 HCC1395 後 hero 與下游 extension 是否同步？ | Browser interaction screenshot | sample_filter=HCC1395 | 1/1；78.8%/77.2%；downstream 1/1 | 非 HCC rows 另有 explicit empty-state assertion |
| `21_multi_bam_dashboard_mobile.png` | Responsive summary | 390px 是否保留 blocker、filter、總覽入口且無全頁 overflow？ | Portable browser screenshot | same artifact | client/scroll 390/390；console/page error=0 | tables 保留局部 scroll；不隱藏資料 |
| `22_multi_bam_dashboard_folded_details.png` | Progressive disclosure | 公式、All-scope 值與缺失欄位能否按需展開？ | Sandboxed no-script HTML block screenshot | scope/availability/axis contracts | 3 details 預設關閉；互動開啟 PASS | 靜態 All 摘要明示 selected sample 以上方卡為準 |
| `23_multi_bam_dashboard_topology_chart.png` | Cohort primary chart | 7 datasets 是否在同一 0–100% 軸保留兩個分母不同的 topology 指標？ | Native chart screenshot + semantic fallback | `topology_rates` | 14 rows＝7 datasets×2 metrics；DORADO 緊鄰 HCC1395 | 描述性比較；不是 clone accuracy |
| `24_multi_bam_dashboard_dorado_selected.png` | Technical control | DORADO 是否明示 replicate、排除 biological macro 且不借用 HCC 下游？ | Browser interaction screenshot | sample_filter=HCC1395_DORADO | datasets=1、replicate=1、macro n=0/n/a、downstream inventory=0/1 | 不宣稱 platform effect 或獨立 biological n |

Repeated screenshot use is intentional: the three screenshots answer different UI-contract questions. Quantitative cohort and missingness questions use distinct chart families from the screenshots.

## 5. Metric contracts and calculations

- `distinct_ssnv`: unique `(chrom, 1-based position)` in topology `active_positions`.
- `region_position_pair`: one active-position membership in one topology region.
- `tree coverage`: regions with representative tree / all regions.
- `unique rate among tree`: unique-best-tree regions / regions with representative tree. It is not unique clone prevalence.
- `ISM linkage`: topology distinct sSNV with a methylation summary / topology distinct sSNV.
- `BH FDR`: audit-time Benjamini–Hochberg recomputation within one axis × one declared scope; it does not prove the original pipeline correction contract.
- `observation scope ledger`: topology distinct sSNV → 有 ISM summary → 至少一個非循環 global axis 可檢定 → 任一 global raw p≤0.05。每層保存 numerator/denominator；cluster 不進此 decision layer。
- `detail availability`: 有可展開的 locus detail / topology distinct sSNV；與 observation ledger 平行，不等同顯著或驗證。
- `axis code provenance`: legacy immutable v1 crosswalk 只含 0–8；`AXIS_UNTESTED=16` 是後續 source contract，不能回填或推論舊 bundle。
- `legacy funnel`: 7,415 regions → 518 with ≥1 A → 472 branch+methyl candidates → 14 displayed；14 是 curated display subset，不是 prevalence。powered linkage 另以 `coread≥6` 與 both-somatic gate 分層。
- `macro/micro`: cohort CSV is per dataset. HCC1395_DORADO is a technical replicate, so seven datasets correspond to six biological samples; no pooled micro rate is used.
- `dashboard downstream inventory`: `downstream_available_count=0` means the selected scope has zero audited bundles after an explicit inventory check；sample-matched ISM/lineage biological fields remain `null` with `ABSENT_NO_EQUIVALENT_BUNDLE`，so inventory zero is never interpreted as zero biological signal.
- `dashboard macro`: unweighted median across rows with `technical_replicate=no`；All scope n=6，single biological dataset n=1；no loci are pooled.
- Version deltas in prose are arithmetic differences of saved rows. They are descriptive because v1 and v3 use different ISM windows.

## 6. Omission reasons

| Requested or plausible metric | Omitted from decision-ready conclusions because | Evidence needed to add it |
|---|---|---|
| Caller precision / recall / F1 | No SEQC2 truth VCF, HighConf BED, TP/FP/FN, som.py/hap.py output in the same provenance chain | Frozen callset + truth + region BED + benchmark command/receipt |
| Causal LCA gain | 7/22 shared chromosomes have different `in_bam`; 22/22 have different thread settings | Same BAM, same parameters, same code/input hashes; LCA flag as only treatment |
| Cluster-axis significance | Groups and tested distance are source-coupled; near-100% significance is double-dipping | Independent labels or held-out validation design |
| Biological clone / subclone claims | A minimum recurrence-allowed representative tree is not a unique biological clone tree | Orthogonal lineage validation and clone-identifiability assumptions |
| Multi-sample ISM/lineage effect | Only HCC1395 has drilldown v1/v3; seven-dataset inventory covers topology, not equivalent ISM/lineage dashboards | Sample-specific manifests and complete per-sample ISM/lineage builds |
| KDE-corrected status | Bundle receipts do not provide sufficient KDE provenance | Explicit preprocessing flag, parameters and hashed artifact |
| Original multiplicity provenance | Saved rows contain p-values but not a complete original multiple-testing family declaration | Persisted q-value, family ID, method and hypothesis count |
| Legacy A/B classification claim | B-first taxonomy makes 3,544/4,025 B loci also allele-high; legacy/current windows/statistics differ | Frozen common locus+allele universe, mutually interpretable class rules and held-out validation |
| Legacy 14-case representativeness | 18 eligible A≥3 regions exist and all four omitted cases have heatmaps; tie-break not encoded | Deterministic selection/ranking contract with inclusion/exclusion reasons |
| Full image delivery QA | Final generated QA deliberately used `igv=skip, panels=0` | Clean full build, file hashes/inventory, image parse and browser interaction checks |
| Broad browser compatibility | Only Chromium desktop/mobile was performed | Shared-reader CI or explicit Firefox/Safari/touch/print checks |

## 7. Draft-time source exclusions

- Claude Code Round 1–4 已納入；Round 5 因 CLI session quota 在任何 Read/Bash 前中止，只保存 execution receipt，不列為 reviewer verdict。既有 rounds 只改變 source/runtime-fix 驗收狀態，不改變 L1 metric、legacy bundle ceiling 或 truth availability。
- No unpublished or inferred truth-set result is included.
- Directory-level “1.7 GB / 4.6 GB” shorthand is omitted from core tables; asset categories from `asset_inventory.csv` are used instead.
- The packaging-time biological macro median/IQR calculation is not used as a headline because it is not stored as a dedicated result row; the report exposes all seven per-dataset values directly.
- **長篇 report** 的原 portable HTML 仍不列為 deliverable；其失敗產物保留於 `debug/renderer_packaging/`，Markdown 是 canonical technical report。另行建立的 **multi-BAM dashboard HTML** 已用 bounded content-width/gutter、mobile legend、filter focus/live-region 與 known-gap copy compatibility patch 修正，並通過 canonical 1440/390 verifier 與 40/40 interaction QA；兩者不可混為同一 artifact。

## 8. QA checklist for the Markdown draft

- Every embedded image path is relative and points to an existing file.
- All top-line bundle, methylation, axis, LCA and cohort values map to saved CSV/JSON rows.
- Original outputs remain immutable; only `InterSubMod/research/20260813_hcc1395_drilldown_validation/` report artifacts are written.
- The report never labels v1/v3 as truth-validated or externally publishable.
- Direct-generated receipt must report `current_mode=generated`, errors `[]`, desktop `1440/1440`, mobile `390/390`, 4 tables, 106 denominator labels, 0 fake cells, and an ISM-bearing machine-selected locus whose methyl section does not fall back to the missing-state message.
- Crosswalk JSON must explicitly state that immutable v1 axis codes predate `AXIS_UNTESTED=16`.
- All legacy headline/funnel counts are labelled raw, powered, both-somatic or displayed; none is called prevalence without a denominator.
- Multi-BAM artifact must retain 18 snapshot datasets、19 sources、10 cards、6 charts、10 tables、2 access issues and snapshot `partial`; portable verifier must report 35 rendered blocks、10 metrics、6 charts、10 tables、3 HTML blocks and pass 1440/390/source dialog；browser receipt must report 40/40 PASS、0 console/page/external-request errors and document widths equal at 1024/512/390/320.
- Browser QA must preserve All 7/6/7/0/0/77.7%/62.4%/1，HCC1395 1/1/78.8%/77.2%/1，and COLO829 downstream inventory 0 while HCC-only tables stay empty；HTML/PNG hashes and console/page error=0 must match the receipt.
- Folded detail block must contain no script，start with three closed `<details>`，and warn that selected-sample values come from dynamic cards.
- **PARTIAL** is attached to the multi-sample section; **BLOCKED** is attached to scientific/external use.
