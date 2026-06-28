<!--
建立時間: 2026-06-27
類型: provenance freeze record — clone/subclone 4 軸教學 HTML 的資料源凍結
build_branch: research/subclonal-reconstruction-202606
-->

# 凍結資料 Provenance（4 軸教學 HTML 用）

> **目的**：4 軸整合教學 HTML 的每個數字，其原始 JSON 都來自**未 merge 的來源 branch**
> `feat/summary-nreadsvalid @ 5308d9e0da8a6fc2d75fcdb77a7f55e8708fa97e`（2026-06-27 14:17）。
> 因該 branch pending-merge，數字在 trunk working tree 原本不可 grep。本目錄用 `git show <hash>:<path>`
> 把已驗證 JSON **凍結 materialize 到 trunk**，使數字落 trunk 可 grep + 可重現（§13-A 由構造防捏造）。
>
> **重現任一檔**：`git show feat/summary-nreadsvalid:<原始路徑> | diff - <本目錄檔>` 應 0 差異。

## 凍結清單（9 檔 + blob hash）

| 凍結檔（本目錄 data/）| 來源路徑（@5308d9e）| blob hash |
|---|---|---|
| sm_locus_master_summary.json | docs/methodology/20260627_clone_subclone_integrated_report/data/ | d65fe64 |
| sm_completeness_ledger.json | docs/methodology/_assets/20260618_subcluster_pilot/ | b8e0b41 |
| sm_summary.json | docs/methodology/_assets/20260618_subcluster_pilot/ | ffe6e50 |
| sm_phaseset_extension.json | docs/methodology/20260627_clone_subclone_integrated_report/data/ | 7b47b4e |
| sm_ccf_tiers.json | docs/methodology/20260627_clone_subclone_integrated_report/data/ | f1f7387 |
| sm_hp_contribution.json | docs/methodology/20260627_clone_subclone_integrated_report/data/ | 9f286cd |
| sm_methyl_corroboration.json | docs/methodology/20260627_clone_subclone_integrated_report/data/ | 54db18e |
| sm_methyl_reextract_ALL.json | docs/methodology/20260627_clone_subclone_integrated_report/data/ | 5310647 |
| sm_methyl_genetic_concordance.json | docs/methodology/20260627_clone_subclone_integrated_report/data/ | 705fb7e |

## 每數字 → 來源 key 對照（HTML 注入用）

| 數字 | 值 | 來源檔 : key | 加總核對 |
|---|---|---|:--:|
| sSNV 宇宙 | 35,332 = TP 30,490 + FP 4,842 | sm_locus_master_summary : src/seqc2 | sum_check_ok=true |
| 完整性帳本 | linked 21,554 / underpowered 5,458 / isolated 8,320 | sm_completeness_ledger : buckets | bucket_sum_check.ok=true |
| γ 類 FP-source linked | 3,204 | sm_completeness_ledger : linked_somatic_by_source.FP | — |
| 區域數 | 7,143 | sm_phaseset_extension : ps_reliability_per_region.total_regions | — |
| 樹形 full_tree | 677 (596+81) | sm_phaseset_extension : by_shape_reliable_vs_uncertain.full_tree | reliable+uncertain |
| 樹形 sibling/linear/co_linked | 1,235 / 1,908 / 858 | sm_phaseset_extension : by_shape_* | 同上 |
| structured / +single-lineage | 3,820 / 4,678 | 衍生：full+sibling+linear / +co_linked | 677+1235+1908=3820 |
| 乾淨 CN full_tree | 205 | region_trees doc §4 表（@5308d9e；無 JSON，doc §13-A 源）| — |
| HP 移除 allelic | 5,238/9,187 = 57.0% | sm_hp_contribution : mutual_excl_ablation | 3949+5238=9187 |
| HP enrichment | mutual_excl 0.86× / nested-indep 1.7–1.87× | sm_hp_contribution : per_relationship_same_hp | baseline 0.5 |
| CCF 梯度 | 支持 2,532(69.8%)/違反 205(5.7%)/tie 890(24.5%) | sm_ccf_tiers : ancestor_ge_descendant_VAF_gradient | 2532+205+890=3627 |
| CCF 決定性 | 92.5% (排除 tie) | sm_ccf_tiers : data_support_rate=0.925 | 2532/2737 |
| GMM | best_n=3, BIC[3990,-1567,-3018,-3001], ΔBIC(3v4)=17 | sm_ccf_tiers : gmm_bic_1to4 | -3018-(-3001) |
| CN-gain confound | somatic 52.8% (12,569/23,810); 乾淨 46.8% | sm_locus_master_summary : cn_somatic_pct | gain/total_somatic |
| PS 可信 | 6,623/7,143 = 92.7% | sm_phaseset_extension : reliable_rate | 含 71 no_ps |
| Tier-PS CCF 一致 | 53,939/127,183 = 42.4% | sm_phaseset_extension : tier_ps_extension | germline≠clonal |
| 甲基(既有ISM) | 8 測 4 = 50%; cis 0/4 | sm_methyl_corroboration : n_tested/n_corroborated | 小樣本 cherry-pick |
| 甲基(重抽ALL) | 740 測 49 = 6.6%; cis 0/740 | sm_methyl_reextract_ALL : n_tested/n_corroborated | 49/740 |
| 甲基獨立 recover | PERMANOVA 可算 1/754, recover 0 | sm_methyl_genetic_concordance : n_testable_permanova | double-dip |

## 已知 own-data 漂移（誠實標註，merge 前待修）

- region_trees doc §0/§7 仍寫「~69% 在 CN-gain」= stale；master TSV 重算為 **52.8%**（06_integrated_narrative §4 round-4 已查並改正，但 region_trees doc 殘留未同步 = 同根因第四次復發）。本 HTML 一律採 **52.8%**。
- structure 4,678 vs 3,820 = 定義差（含/不含 858 單 lineage），非矛盾；HTML 兩者並陳。

## 2026-06-28 補充：output-completeness 修補 + 模型驗證 + somatic 原則

**新凍結/衍生（消除 prose-only + 硬編 §13-A 違規）**：
| 檔 | 來源 | 內容 |
|---|---|---|
| `data/regions.tsv` | git show 5308d9e（branch tracked）| 7,143 區逐區明細（region/chrom/start/end/span/n_sSNV/cn/tree_shape/n_nodes/n_nested/n_sibling/max_depth/has_cycle/n_populations）|
| `data/region_shape_distribution.json` | `scripts/derive_region_buckets.py` derive 自 regions.tsv | 6-bucket（含**前缺的 no_confirmed_structure 2443 / inconsistent 22**）+ cn 分布 + n_populations 分布 |
| `data/clean_subset.json` | 同上 derive | 乾淨 CN(LOH+neutral) 各 shape：**full_tree 205**（取代 build script 硬編）/ linear 763 / sibling 408 / co_linked 357 |
| `data/sm_configuration_census.json` | git show 5308d9e | 2×2 細胞模式 census（mutual_excl/nested/co_linked... same/diff HP）|
| `data/chr17_subclone_data.json` `chr17_tree_data.json` | git show 5308d9e | chr17 per-read 機器真值（**3 sSNV、9/20/19**）|
| `scripts/sm_*.py` `derive_region_buckets.py` | git show + 本地 derive | 計算 sSNV 共現的實際 pipeline 原始碼 |

derive 錨點驗證全通過：6-bucket 加總=7143、structured 3820、+colinked 4678、clean full_tree=**205**（符 region_trees §4）。

**somatic 定義原則（2026-06-28 用戶定案）**：build 用 ClairS→longphase-S 的 **TP∪FP union**（`load_union` 讀 `filtered_snv_{tp,fp}_*.vcf.gz`），somatic 與否由 **normal 比對**（`is_somatic` normal-VAF）決定；**SEQC2 的 TP/FP 標籤只用於觀察/評估，絕不進前處理或定義**。經查程式碼**已遵守**（src 只記錄不當 filter）。唯一待修＝genome-wide `<5%` vs chr17 builder `==0` 不一致（兩者皆 normal-based、皆非 SEQC2）→ 統一後 chr17 = **3 個 somatic sSNV**。

**模型驗證裁決全文** → `InterSubMod/docs/methodology/20260628_reconstruction_model_verification_01.md`（workflow wf_f2b070ea-64c，15 agent）。核心修正：CNV/LOH = condition-on confound + 輔助刻畫（非刪除、非獨立證據 rung）；VAF = consistency-check 非獨立；HP = 鑑別器只在互斥；主張3（甲基判時序）INVALID；classify() 噪聲容忍不對稱 + somatic 雙定義 = 待修 code bug。

**門檻定案（2026-06-28，ε=2%）** → `InterSubMod/docs/methodology/20260628_sSNV_linkage_threshold_decision_eps2_01.md`。cell 為真 ⟺ `count > coread×0.02`；三路驗證（FP 裁判 +17.9 / 結構穩定 full_tree 677→616 / 塌陷集中 coread 74-vs-27 CN-gain 74%）。新產物：`data/`{`pairs_eps2_annotated.tsv`(逐對可手驗)、`eps2_final_band.json`、`threshold_comparison.json`、`region_threshold_impact.json`、`weak_pair_observation.json`、`fp_and_vaf_threshold.json`、`sSNV_combination_enumeration.json`、`region_shape_distribution.json`、`clean_subset.json`、`regions.tsv`、`lists/*.tsv`(10)}；`scripts/`{`apply_eps2_canonical.py`、`compare_thresholds.py`、`region_threshold_impact.py`、`observe_weak_pairs.py`、`fp_and_vaf_threshold.py`、`enumerate_sSNV_combinations.py`、`derive_region_buckets.py`、`diagnose_problems.py`}。band: co_linked 11,750→16,048、nested 13,113→11,352（重分類 15.7%）。

## 對抗稽核歷史（來源 branch，5 輪）

06_integrated_narrative §4 記錄 4 輪 fresh-context workflow NEEDS_WORK → 全數修正（baseline 0.443→0.5、HP 框架=鑑別器、CCF 含 tie 誠實口徑、GMM 取代目視、PS 分母含 no_ps、甲基 n=8→重抽 740、64%→52.8% own-data 矛盾）；+ trunk 統一敘述 workflow wck0bu3iq（8 agent，0 真矛盾）。三輪共通根因 = 修正只落 prose 未同步 JSON/HTML（§13-A 違規）→ 本 HTML 由 generator 從凍結 JSON 注入，杜絕該根因。
