<!--
建立時間: 2026-07-27 09:37
更新時間: 2026-07-27（Task B 全範圍驗證與 read 聚集圖完成）
目標: 預先凍結多 sSNV R/A/X pattern 與 read-level 甲基關聯的全面驗證設計
處理範圍: 7 technical datasets / 6 biological samples / chr1-22 / exact PS × exact raw HP
關聯檔案:
  - InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/pre-decision-audit.md
  - InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/implementation-notes.md
  - InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/analysis-contract.md
  - InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/20260727_多sSNV_pattern與甲基關聯全面驗證_01.md
  - InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/
-->

# 多 sSNV pattern × 甲基關聯全面驗證

**Created**: 2026-07-27
**Status**: validated
**Task type**: B — Comprehensive validation
**Scope**: 7 technical datasets、6 biological samples、chr1–22
**Primary grain**: `(dataset, chrom, exact phase_set, exact hp_raw, region/unit)`

> **Claim ceiling**：本專案只驗證 LongPhase-S callset-anchored、
> pattern-conditioned regional methylation association。任何結果均不能單獨證明
> cellular subclone、mutation ancestry、因果方向或唯一 lineage tree。

## 1. Pre-registration

| 假說 | 否證條件 | Decision threshold |
|---|---|---|
| H1：至少存在一個 exact raw-HP full-four-state unit，其 R/A state 可解釋 regional Bernoulli methylation structure | report-wide frozen family 中沒有任何單元同時通過 effect、BY、PERMDISP 與 confound gates | `q_BY≤0.05`, `R²≥0.10`, standardized effect `≥0.5`, distance contrast `≥0.10` |
| H2：至少存在一個 state（如 RA）呈現 within-state methyl cohesion，且與一個以上其他 state 分離 | 所有 eligible state 的 within-vs-between effect 在等 N、共同 CpG rarefaction後消失或反向 | primary `n/state≥5`, `N≥40`；`n≥8/10` sensitivity 同方向 |
| H3：較長的 complete R/A signature 可保留 pairwise 以外的多點關聯訊息 | k≥3 frozen family 無任何 claim-eligible unit，或全部只由 geometry/strand/RG/dispersion 解釋 | unit-level restricted tests；X 只作 subcube 描述 |
| H4：sample-level finding 可跨 biological samples 重現 | claim-eligible finding 只出現在一個 biological sample | 至少 3/6 biological samples 同方向才稱 cross-sample；HCC1395 technical pair只算 n=1 |

不在上表的 locus ranking、狀態組合與圖例均標 `exploratory`，不得事後升格為
confirmatory family。

## 2. Frozen candidate families

1. **Primary full4**：同一 exact PS × exact `hp_raw` 中，RR/RA/AR/AA
   每態 `n≥5` 且總 `N≥40`；`n≥8` 與 `n≥10` 為敏感度。
2. **Secondary pair contrast**：至少兩個 complete states各 `n≥5` 且總
   complete `N≥40`；`n≥8` 為 robustness gate、`n≥10` 另報。此 family
   只作 comprehensive exploratory census，不升格 robust claim。
3. **Long complete signature**：k≥3、full-cover `N≥40`、至少兩個 complete
   signatures 各 `n≥5`。
4. **Partial R/A/X**：保留為 subcube/terminal-interval evidence，不填補 X，
   不映射成完整 topology node。

候選只由 sSNV、read coverage、PS/HP 與 callability決定，不得用甲基 effect
先挑位點。

## 3. Analysis contract

- marker 按 genomic coordinate 排序；R=REF、A=ALT、X=非固定/未觀察。
- methyl window 為每個 marker ±5 kb 的 union；marker 間無 gap cutoff。
- 所有待比 state 使用同一 unit-level callable CpG basis。
- Bernoulli read-pair 最少 3 個共同 CpG；unit 至少 10 個 distinct callable CpG；
  每個 state cell finite coverage ≥80%。
- invalid distance 保留 NA；禁止轉成 1.0。
- restricted permutation strata = exact `strand × read_group`；不可交換則
  `NOT_EVALUABLE`，禁止 fallback 到 unrestricted。
- primary 另做等 N downsampling、共同 CpG rarefaction、marker-created CpG /
  marker ±100 bp masking 與 distal sensitivity。
- confirmatory與secondary frozen families分開使用 BY，並各自以 Holm作
  family-wise sensitivity；state-pair只回報 effect與topology relation，
  不另作未預註冊的 pairwise顯著性宣稱。
- PERMANOVA 報 pseudo-F、R²、requested/realized permutations；PERMDISP
  顯著者不得列 robust association。

## 4. Topology projection

- complete R/A signature 可映射既有 mutation-state node。
- Hamming-1 且已由 global-best candidate family unanimous 支持的 transition，
  才可加 methyl association halo；halo 不改 edge incidence、AF score、support
  或 selected tree。
- Hamming >1 只畫 pair band/table。
- topology 與 methyl authority分開綁定、分開出 receipt。

## 5. 研究啟動五問

1. **Thread D**：是，read-level epigenetic 主線。
2. **Thread B 撤回範圍**：否；不做 TO filter，但沿用其 HP/self-phasing警告。
3. **KDE-corrected**：本題不使用 Coverage_Multiple/KDE 作主要統計；
   若加入 CN/LOH annotation，僅使用 corrected authority。
4. **VCF caller AF**：topology overlay沿用 frozen caller-AF authority；
   methyl association本身不以 AF 選樣。
5. **長計算/C++/搬移/NO-GO**：是；新輸出只寫 observation workspace，
   不修改 archive/canonical payload。

## 6. Step → Verify

1. 建 input/artifact catalog → 驗證：154/154 sparse files、read digest零碰撞、
   每個 marker tile有 attested source或明確 `MISSING`。
2. 建 exact raw-HP pattern census → 驗證：所有 pattern由同一 molecule row重算，
   family-collapse sensitivity與 primary分開。
3. join methyl long-form → 驗證：qname SHA-256一對一、同一 read metadata一致、
   tile union不拼接既有不同 read-universe distance matrix。
4. 跑統計 → 驗證：deterministic seed、restricted permutation realized count、
   BY/Holm與 confound gate可重算。
5. 產生 report/workstation sidecar → 驗證：數字由 TSV/JSON注入、association-only
   標籤、桌機/手機/print QA。

## 7. 研究目標

- **G3**：把 ISM read-level methylation 與多 sSNV molecular pattern正式連結。
- **G4**：全 7 technical datasets × chr1–22與可重跑 receipts。
- **G5**：公開 schema、claim ceiling、confound gates與可審核 HTML。

## 8. Validated result

Task B 全範圍執行完成：

- 21,601,606 sparse rows → 8,893,098 candidate read projections →
  154,132 exact raw-HP strata。
- Formal `n≥5/n≥8/n≥10` units = `1,045/957/910`；full-four = 0；
  `k≥3` = 162。
- 811 evaluable；3 `ROBUST_ASSOCIATION`、627
  `EVALUABLE_NO_ROBUST_ASSOCIATION`、181 `CONFOUNDED`、234
  `NOT_EVALUABLE`。
- 三個 robust 單元：H1437 chr22 `AARR/RRAR/RRRA`、H2009 chr3
  `AAA/RAA`、HCC1937 chr10 `AAR/RRA`。
- Exact 二位點 RA 對 AA/AR/RR 的可評估 memberships 分別
  `266/186/3`，robust 均為 0。
- Artifact catalog 2,313/2,313 PASS；Bernoulli parity 2,313/2,313 PASS；
  Python 88/88、C++ 258/258、report 18/18、browser QA PASS。
- 557/1,045 units 有原始 read × read matrix；已加入 label-independent
  UPGMA 聚集圖。其餘 488 units 明示 matrix 不可用原因，不以 state mean
  補造 read matrix。

交付：

- [完整結果與解釋](20260727_多sSNV_pattern與甲基關聯全面驗證_01.md)
- [互動 HTML](../../output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/20260727_multisite_pattern_methyl_association_06.html)
- [Browser QA](../../output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/qa_v06/browser_qa.json)
- [Sidecar receipt](../../output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/sidecar/pattern_methyl_sidecar.v1.receipt.json)
