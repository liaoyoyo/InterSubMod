<!--
建立時間: 2026-07-23 17:00 Asia/Taipei
目標: 記錄 HCC1395 biological-ID specificity 與多位點融合樹 HTML 報告的來源、圖表契約、命令與 claim ceiling
處理範圍: 7 technical datasets × chr1-22；21 組配對；HCC1395 × HCC1395_DORADO 為唯一 same-ID positive pair
關聯檔案:
  - InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/biological_id_specificity_report/build_artifact.mjs
  - InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/biological_id_specificity_report/deliver_portable_report.mjs
  - InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/biological_id_specificity_report/artifact.json
  - InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/biological_id_specificity_report/data/report.sqlite
  - InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/biological_id_specificity_report/20260723_HCC1395_biological_ID_specificity與融合樹驗證_01.html
-->

# HCC1395 biological-ID specificity 與融合樹驗證：來源筆記

## 任務分類與研究目標

- Task Type：B（comprehensive validation）
- Scope：7 technical datasets × chr1–22；完整 21 組 pairwise comparisons
- 服務目標：G4 多樣本一致性與 reproducibility；G5 可被外部驗證
- 主要問題：同 biological ID 的局部結構是否顯著比不同 cell lines 相近，以及目前證據是否足以升格為 global clone/subclone fusion tree。

## 假設與判定

1. HCC1395 × HCC1395_DORADO 是唯一預註冊 same-ID positive pair；其餘 20 組是 negative controls。
2. Cross-dataset comparison 排除 PS/HP label orientation，只比 exact allele、無向 read-linkage 與 component co-membership。
3. 高 conditional concordance 必須同報 parent-universe coverage。
4. SEQC2 Figure 4 只作 reference architecture，不作本地 exact directed-tree truth。
5. 同 cell line 不保證同 aliquot、passage、library 或 molecules。

## 主要輸入

- `InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/strict_pair_validation/results/strict_pairwise_metrics_01.json`
- `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage全資料集報告_01/artifact.json`
- `InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/report/artifact.json`
- `InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/clone_region_bridge_v1/summary.json`
- KB：`kb-02-samples-hcc1395`（runtime-fact, verified, last_verified 2026-07-11）
- KB：`kb-04-databases-seqc2-truth-set`（reference-summary, verified, last_verified 2026-07-12）
- 外部主來源：Fang et al. 2021, PMCID PMC8532138；官方 SEQC2 v1.2.1 README。

四個定量輸入已由 builder 以 SHA-256 fail-closed 綁定：

| Input | SHA-256 |
|---|---|
| strict pairwise metrics | `aaadfa8751fbe6b78efb4b5070aeef9f6bbdfc14e67755bad52b6b9d443cfa41` |
| exact-PS all-dataset artifact | `6497253851cd3fb3a6f043d90aee3685ae33ecd4f0cc084cd268325df04154a6` |
| raw-VAF / PyClone artifact | `1094b39be9d08a5314224a7ee9a56619a36fb22a6e4aacbf5def95b4088109ad` |
| clone-region bridge | `8a8b489ea5b77bd19e6d1160374b77fdddac0eac3cdb760347603bb1f83ba3c8` |

Pairwise provenance：HCC1395 與 HCC1395_DORADO receipts 同用 pre-hotfix builder `7260a763…`、graph core `df3d6f…`、extractor `2ca7ccb…`，可降低 pair 內 code-version confounding。All-7 cohort 並非完全同 builder；HCC1954 使用 post-hotfix `912721f9…`。Data-specific no-trigger equivalence 支持本次數值未被該修補改變，但不證明演算法正確。

## Chart map

| Section | Analytical question | Family / type | Fields | Supported takeaway | Palette |
|---|---|---|---|---|---|
| biological-ID specificity | same-ID 是否遠高於 unrelated？ | comparison / grouped bar | layer × pair_class → Jaccard | HCC pair 七層 rank 1/21 | blue + gold |
| shared structure | 低 Jaccard 是錯邊還是解析不對稱？ | comparison / grouped bar | layer × measure → fraction | 與 DORADO 為較稀疏局部 subset 相容 | blue + olive |
| phase resolution | 差異是否符合 PS fragmentation？ | comparison / bar | metric → DOR/HCC ratio | 更多 PS boundaries、較短 linkage reach | single blue |
| denominator audit | 高一致是否有足夠 coverage？ | uncertainty / grouped bar | evidence × measure → fraction | conditional fidelity 高、global coverage 低 | blue + gold |

## 關鍵數值

- Candidate alleles：76,721 shared；Jaccard 92.76%；rank 1/21。
- Shared-candidate projection：DORADO→HCC direct-edge containment 98.28%；co-membership containment 90.05%。
- Joint-W coverage：10,154 / 76,721 = 13.23%。
- Raw VAF：Pearson 0.9343；CCC 0.9339；MAE 0.0624。
- Independent PyClone-VI：ARI 0.5389；κ 0.5442；minor-set Jaccard 0.3810。
- Direction bridge：596 / 598 same direction，但 jointly determinate coverage = 598 / 8,096 = 7.39%。
- Informative multi-cluster regions：21 / 34 partition exact = 61.76%；34 / 5,720 = 0.59% parent-universe coverage。

## Claim ceiling

可主張：

- HCC1395 same-ID pair 的 exact allele、raw VAF 與局部無向 read-linkage skeleton 具有很強、且不同 biological ID 不具備的重現性。
- 現有結果與 gross sample mix-up、all-pairs homogenization 或 DORADO 大量選到與 HCC 衝突的局部無向 edges 不相容，但不能排除較細微 processing bias。
- Shared-opportunity 上的結果與 DORADO 是較稀疏局部無向 subset / contraction 的解釋相容；這不是已證明的 directed topology。

不可主張：

- 完整 global clone tree 已相同或已由 SEQC2 Figure 4 證實。
- 98.28% 是 parent→child edge accuracy 或全域 coverage。
- PyClone-VI cluster 是演化樹 truth。
- 同 cell line 的殘差全是演算法造成。

## Step → Verify

1. 讀取四個既有 artifacts → 驗證：四個 SHA-256 與凍結值一致、strict `all_validations_pass=true`、HCC pair 與 DORADO rows 存在。
2. 產生 canonical `artifact.json` → 驗證：第一個 markdown `#` 與 `manifest.title` 一致；snapshot 為 plain row arrays。
3. 呼叫 Data Analytics `validate_artifact` → 驗證：validator 回傳 valid。
4. 執行 packaged shared renderer、static-chart extractor 與 browser verifier → 驗證：receipt 的 `stages.verification` 為 `passed`；HTML 無外部 runtime dependency。

## 執行命令

```bash
node InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/biological_id_specificity_report/build_artifact.mjs

node InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/biological_id_specificity_report/deliver_portable_report.mjs \
  --input InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/biological_id_specificity_report/artifact.json \
  --output InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/biological_id_specificity_report/20260723_HCC1395_biological_ID_specificity與融合樹驗證_01.html
```

## 實際輸出與驗證 receipt

- Canonical artifact：`InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/biological_id_specificity_report/artifact.json`
  - 13 datasets；88 reviewed rows；20 sources；29 blocks；4 charts；9 tables。
  - Data Analytics `validate_artifact`：`ok=true`、`snapshot_status=ready`。
  - SHA-256：`87a3098f311bf059389432861d6b672126921cce41f4e82c4d2aa4d98a771897`
- Executed-query database：`InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/biological_id_specificity_report/data/report.sqlite`
  - 13 個實際 SQLite 查詢；每個 chart/table source 均保存 SQL。
  - SHA-256：`190e043ae8e05ac086a1aeb86c7d8e033b0355a7ddfdbfb96e8c50d3645265c0`
- Portable HTML：`InterSubMod/research/20260723_hcc1395_crosssource_topology_resolution/biological_id_specificity_report/20260723_HCC1395_biological_ID_specificity與融合樹驗證_01.html`
  - Shared renderer validation/package/static charts/browser verification：全部 `passed`。
  - Rendered counts：29 blocks、4 charts、9 tables、2 custom SVG。
  - Viewports：1440 px 與 390 px；source dialog 與 keyboard semantic click：`passed`。
  - SHA-256：`56ce740085bb1f3773ad840e3f1e0d363703fe6b4aed8f3f44ba46ecb94c4ca1`

## Renderer workaround

原始 `deliver_portable_artifact.mjs` 的 enhanced-reader QA 找到 shared runtime top bar 在有垂直 scrollbar 時因 `width: 100vw` 多出 8 px，錯誤碼為 `horizontal_overflow`。`deliver_portable_report.mjs` 沒有另寫 renderer；它仍呼叫同一 packaged `buildPortableArtifact`、`extractPortableChartSvgs` 與 `verifyPortableArtifact`，只對 outer document 與 packaged runtime 注入 `html,body{overflow-x:clip!important}`。修正後 desktop、mobile、source interaction 與 external-request checks 全部通過。
