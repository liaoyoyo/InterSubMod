<!--
建立時間: 2026-07-12T04:54:04+08:00
目標: 定義 HCC1395 逐位點拓撲 containment 技術報告的敘事、圖表與 artifact 欄位
處理範圍: HTML technical report；chr1-22；5,720 exact-coordinate complete-both regions
關聯檔案:
  - InterSubMod/research/20260712_hcc1395_pair_site_topology_containment_validation/pre-decision-audit.md
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_matches.tsv
-->

# Report plan：HCC1395 逐位點 topology cross-source validation

## Reporting job

- Audience：**technical**。
- Delivery mode：**portable self-contained HTML**；canonical `artifact.json` → packaged Data Analytics report builder。
- Question：兩個 HCC1395 cross-source technical datasets 是否在相同 genomic sSNV 上重現 exact 或 induced-substructure topology？
- Decision-useful answer：區分 read full candidate space、VAF top set、coverage、compatibility、conflict、conditional null 與 claim ceiling。
- Fixed denominator：5,720 exact-coordinate complete-both region pairs。
- Success：site-aware robust/nested signal高於 complexity-matched null，且不只由 k=2、大 candidate set 或 HP mismatch 驅動。

## Visible structure（technical-report mapping）

1. Title + visible PARTIAL ribbon。
2. Technical summary：直答有多少 exact／contained／overlap／conflict／not-evaluable，read 與 VAF 分開。
3. Key findings：fixed-denominator outcome、candidate uncertainty bounds、complexity strata、null。
4. Scope / definitions：same-cell-line cross-source、shared-site projection、HP mapping、candidate-set relation、VAF heuristic。
5. Methodology：solver re-enumeration、recurrence/hidden handling、conditional shuffle、block bootstrap。
6. Limitations / robustness：small-k、candidate multiplicity、HP split/merge、VAF comparison-count、historical-v2、no biological truth。
7. Next steps：clean-v3 closeout、depth matching、split-half、single-cell/synthetic truth。
8. Further questions：跨 cell line algorithmic-prior control、CN/purity-corrected ordering。

## Chart map

| Section | Question | Family / form | Fields | Supported claim | Palette |
|---|---|---|---|---|---|
| Fixed-denominator outcome | 5,720 區各落在哪一類？ | Composition / 100% horizontal stacked bar | layer, mapping, outcome, n, share | coverage-adjusted yield；read vs VAF | multi-category：deep blue / cyan / gold / orange / grey，另有文字標籤 |
| Candidate uncertainty | any-match 與 reciprocal/full-set 差多少？ | Uncertainty / dot + interval | layer, lower, upper, n | candidate multiplicity 形成 optimistic upper bound | hard two-root：blue + neutral |
| Complexity strata | signal 是否只來自簡單區？ | Comparison / faceted dot | k_bin, HP pair, candidate bin, rate, n, CI | k=2／大型候選不應掩蓋複雜 strata | single-root blue shades + open markers |
| Conditional null | observed 是否高於結構先驗？ | Uncertainty / observed dot + null interval | metric, observed, null_mean, q025, q975, p | same-region pairing超過 matched random | blue observed + gold null，shape区分 |
| Site-set inventory | 子結構規則實際適用多少區？ | Composition / horizontal bar | relation, n, share | 96.75% equal；relaxation只直接處理186區 | single-root blue shades |
| Representative regions | exact／contained／conflict各長什麼樣？ | Audit table + small topology diagrams | region, sites, HP, read/VAF outcome, reason | 解釋分類，不作抽樣估計 | neutral + status glyph |

所有圖：零基準（比例圖 0–100%）、明示 n/denominator、相鄰解釋段、canonical source metadata、mobile stacked layout；不以綠紅作唯一區分。

## Required machine-readable artifacts

- `data/hcc1395_site_topology_pair_outcomes.tsv`
- `data/hcc1395_site_relation_evidence.tsv.gz`
- `data/hcc1395_topology_compatibility_metrics.tsv`
- `data/hcc1395_topology_complexity_strata.tsv`
- `data/hcc1395_topology_null.tsv`
- `data/hcc1395_topology_examples.tsv`
- `data/hcc1395_topology_containment_summary.json`
- `data/hcc1395_topology_containment_checks.tsv`
- `data/hcc1395_topology_signature_profiles.jsonl.gz`

## Claim language

- Allowed：same-cell-line cross-source technical reproducibility；read-pattern-compatible candidate space；VAF-supported heuristic top set；shared-site induced substructure；高於 matched null。
- Forbidden：true clone tree、biological accuracy、independent biological replicate、VAF confirmed ancestry、癌症基因／藥物註記證明 topology。

## Provenance Footer

- Commit hash：`6067568637088838a9f518955e41d222f057e4f1`
- Build time：2026-07-12T04:54:04+08:00
- Skills：`data-analytics:build-report`、`data-analytics:visualize-data`、`frontend-design`

