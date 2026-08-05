<!--
建立時間: 2026-07-18 23:45 +08:00
目標: 全量驗證 positional-singleton sSNV 中 focal-ALT read 的甲基子結構，並建立 HCC1395 深入 HTML 報告
處理範圍: 7 datasets / chr1-22 / 50,432 singleton dataset-sites；HCC1395 8,279 sites 深入分析
關聯檔案:
  - InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/pre-decision-audit.md
  - InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/implementation-notes.md
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/implementation-notes.md
-->

# Positional-singleton focal-ALT 甲基子結構驗證

> **Task Type B — Comprehensive validation。** 本 topic 使用已完成且 source-attested 的
> 50,432-site singleton audit 做完整母體重算，HCC1395 的 8,279 個 singleton sites 不取 subset。
> 報告可確認 read-level residual epigenetic partition；正式 genetic cooccurrence、matched-normal、
> CN/purity/CCF 與 cellular lineage 尚未完成，因此不可把 M2 PASS 直接命名為 clone/subclone。

## Pre-registration（Confirmatory）

| 預測 H | 否證條件 | Decision threshold |
|---|---|---|
| HCC1395 singleton 母體可精確重現 8,279 sites，且 M1/M2 計數與 source-attested summary 一致 | 獨立重算任一核心計數或 digest 不一致 | 任一不一致即停止發布 HTML 數值結論 |
| M2 PASS 位點可在原始 focal-ALT read 距離／甲基矩陣中重建兩群視覺結構 | assignment read id 無法一對一 join 原始矩陣，或群組數／read 數不一致 | 任一目標位點 join 不完整即降為資料缺口，不畫 heatmap |
| M2 measured-axis guardrail 排除 HP/strand/read-geometry/CpG-count 的已測對齊 confound | 任一已測軸為 aligned，或 M2 axis 不 determinate | 該位點不得列為 clear residual partition |
| 單 sSNV 無足夠遺傳資訊證實 cellular clone/subclone | 未完成第二獨立 genetic marker、matched-normal/REF、CN/purity/CCF | confirmed clone/subclone 固定為 0（未驗證，不是真陰性） |

## 研究目標對應

- **G3**：確認 read-level epigenetic partition 的可觀察性、限制與區域診斷。
- **G4**：同一完整 contract 下報告 7 datasets，並對 HCC1395 做可重現深挖。
- **G5**：提供 source digest、可重跑 script、machine-readable results 與 standalone HTML。

## Scope 與資料契約

- 全 sSNV screen：469,849 dataset-sites，LongPhase-S recalibrated `FILTER=PASS`。
- positional singleton：同 dataset/chrom、相鄰 gap `>50 kb` 才切開的傳遞 component，component size=1。
- 統計單位：`(dataset, chrom, pos, ref, alt)`。
- HCC1395 singleton 母體：8,279 sites。
- M1：`coarse_ng >= 2`、非 unstable、`modal_assignment_ari_min >= 0.8`。
- M2：8 個 measured axes（HP exact/family、strand、read start/end/length、MAPQ、CpG-called）均需可判定且無群組對齊 confound。
- Claim ceiling：`M2_read_level_residual_epigenetic_partition`。

## 可重現性

- Git HEAD：`0ee2fa1b31fcf6af670efd301251b5b3a24c1a99`
- 上游 audit summary SHA-256：`6d4128be0535f2e16b1382cf038c558054ac42d0d7bde75ab4854b7a5be7aedc`
- 上游 site audit SHA-256：`2d5d24790918fca34a32c313b9965f1de8c186c031e31d1a643febba721d90ce`
- 上游 M2 PASS cases SHA-256：`2fe42d55924d4d73f8a0d7b436cc8c517e9d06834da9fd7a1e882ce24c4f1dd0`
- 隨機性：不使用隨機抽樣；heatmap display reads 使用 deterministic medoid-nearest selection。
- Next spaced check：2026-08-17。

## Step → Verify

1. 重算 50,432-site 與 HCC1395 8,279-site 母體 → 驗證：dataset-level 計數逐欄等於 frozen summary。
2. exact join HCC1395 兩個 M2 PASS loci → 驗證：108/108 與 109/109 core reads 都存在於 read、甲基、距離矩陣。
3. 產生統計、read-distance 與 methylation heatmap datasets → 驗證：無非有限值；矩陣對稱／對角為 0；顯示抽樣規則可重現。
4. 建立 artifact + standalone HTML → 驗證：artifact validator PASS；portable builder browser verification PASS。
5. 獨立重算核心比例與 claim ceiling → 驗證：receipt 列出 exact numerator/denominator 與 NOT_RUN gates。

## 交付物

- Standalone HTML：
  `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/20260718_HCC1395_singleton_ALT甲基子結構驗證報告_01.html`
- 完整 HCC1395 8,279-site audit：
  `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/hcc1395_singleton_site_audit.tsv.gz`
- 七資料集摘要：
  `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/all_dataset_singleton_summary.tsv`
- 兩個 HCC1395 M2 PASS 位點摘要：
  `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/hcc1395_m2_pass_locus_summary.tsv`
- 分析驗證收據：
  `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/validation_receipt.json`
- HTML 瀏覽器 QA 收據：
  `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/html_qa_receipt.json`
- Heatmap rendering QA 與六張視覺截圖：
  `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/heatmap_rendering_qa.json`
  與
  `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/heatmap_qa/`

## 結論層級

- **可確認**：HCC1395 有 734 個 singleton sites 通過 M1 穩定多群候選；其中 2 個通過
  已量測 HP、strand、read geometry、MAPQ 與 CpG-called guardrails，屬清楚的 ALT-read
  residual epigenetic partition candidates。兩點各有 distance、shared-CpG 與 methylation
  彩色 heatmap；每張顯示 8 個 deterministic representative ALT reads，但 effect summary
  使用全部 108/109 core reads。
- **不可確認**：這 2 個位點是否是 cellular clone/subclone、其數量、父子順序或唯一真實演化樹。
  原因是 singleton 不含 local second sSNV，且 G1/G2/formal R1、matched-normal、CN/purity/CCF
  與 cellular lineage 均尚未執行。

## Independent derivative sidecar v4

本節僅新增另一條獨立、no-clobber evidence path，不取代上方 parallel implementation：

- Technical HTML：
  `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/report.html`
- Canonical artifact：
  `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/artifact.json`
- Machine summary / denominators：
  `InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/derivative_validation_v4/machine_summary.json`
  與 `census_breakdown.tsv`、`status_denominators.tsv`。
- Target join / case / display-selection audits：
  `target_join_audit.tsv`、`case_statistics.tsv`、`representative_selection_audit.tsv`。
- Source inventory / chart map：
  `source_inventory.json`、`source_inventory.tsv`、`chart_map.json`、`chart_map.tsv`。
- Build / validation / browser receipts：
  `build_receipt.json`、`canonical_validation_receipt.json`、
  `portable_builder_compat_receipt.json`、`browser_qa_v3/browser_qa_receipt.json`、
  `focused_test_receipt.json`。
- Final anchors：artifact SHA `32183456381e8e0ec6fae5dfe19581ec7730c1909d565e92b3365b15b9649f7d`；
  HTML SHA `31981e8f3ff374e26b27ea7e497f2c2b98f637a8a7d1187dcb3ab2e53adcc371`；
  tests `22 passed`；desktop/mobile QA PASS。
