---
id: ism-kb-11-full-method-process-external-validation-20260712
name: "InterSubMod full method/process external validation — 2026-07-12"
description: "濃縮現行 raw-all→LongPhase-S→partial-read regional tree 全流程、外部一致/差異/衝突、證據天花板、缺口與應用路線，供後續研究與論文敘述直接查詢。"
status: active
last_verified: 2026-07-12
content_nature: reference-summary
doc_type: reference
verified_scope: "7/7 raw-all producer/receipt closeout + fresh layered-v3 2026-07-12 18:16 partial status + external primary literature through 2026-07-12; biological clone validation remains absent"
decision_refs: [ISM-A01, ISM-M02, ISM-M03, ISM-M05, ISM-M06, ISM-M07, ISM-R01, ISM-R03, ISM-R04, ISM-R05, ISM-R06, ISM-D01, ISM-D02]
related_ids: [ism-kb-11-external-literature-index]
tags: [external-literature, method-validation, raw-all, longphase-s, partial-reads, regional-tree, cn, ccf, methylation, claim-boundary]
canonical_paths: [11_external_literature/14_20260712_full_method_process_external_validation.md]
alias_paths: []
---

# InterSubMod 全方法／流程外部驗證 — 快速知識卡

## 一句結論

現行方法在**窄版 estimand**上與外部研究一致：從long-read sSNV linkage在parental HP family內列舉recurrence-free regional mutation-state candidate trees並明示拒答；目前只完成7/7 producer contract，fresh downstream與biological clone truth尚未完成，因此正式Results與clone-level claim仍NO-GO。

完整報告：[/big7 external landscape 12](/big7_disk/liaoyoyo2001/external_validation/_landscape/12_intersubmod_full_method_process_external_validation_20260712.md)

## 現行 canonical chain

```text
normalized paired ClairS raw-all biallelic sSNVs
  → LongPhase-S bidirectional recalibration + alignment HP/PS sidecar
  → same-run _sc.vcf FILTER=PASS (canonical tree input)
  → positional region × HP1/HP2 family × partial/full R/A/X read states
  → all minimum compatible recurrence-free regional mutation-state trees
  → determined / ambiguous / capped / incompatible
  → CN/read-AF/methylation bounded post-hoc interpretation only

ClairS PASS → sensitivity/audit subset, not the complete recalibration universe
```

## 2026-07-12 operational truth

- raw-all producer/receipt：7/7 PASS。
- 638,259 input records；582,820 `_sc PASS`；LowQual→PASS 32,184；PASS→LowQual 17,444。
- 164,253,537 sidecar rows；unknown HP=0；identity conflict=0。
- fresh layered-v3在18:16 snapshot為4/7 sample完成、2 active、1未啟動，root無`_SUCCESS`；四個完成sample的71,783 non-capped eligible units通過V1–V7，但未經run-level final verifier。
- sensitivity、comparison、final independent verifier、post-run immutable receipt尚未完成。

## 最重要方法更正

1. **不是每棵樹都由full-span read直接觀測**：region先按positional gap成component，solver可由互相重疊partial reads完成約束。寫作必用`partial-read-constrained, model-determined regional mutation-state tree`。
2. **HP family不是clone**；PS只作phase-block QC。
3. **唯一相容解不是唯一生物史**；finite-sites violation、CN/loss、calling/mapping/phasing error與sampling仍可造成替代解釋。
4. **CN仍是post hoc**：只取region midpoint解讀recurrence；read AF未加purity、allele-specific CN、multiplicity時不是CCF。
5. **methylation目前未進fresh run**，也只能bounded annotation。

## 外部一致、張力與真正衝突

| 類型 | 知識 |
|---|---|
| 一致 | ClairS/LongPhase-S支持paired、phase-aware upstream；SomaticHaplotype/Vo/Sakamoto支持physical linkage；PhyClone/Pairtree支持顯式uncertainty；CNAqc/DeCiFer支持CCF校正；Rosenski/Blood支持methylation confound controls。 |
| 差異 | PhyClone/Pairtree估global cellular clones/frequencies；InterSubMod估regional per-HP mutation states。Bayesian posterior與exact compatible-set不是同一物件。 |
| 衝突1 | TreeClone 2019已用same-read mutation pairs推subclone phylogeny，所以不能宣稱首次使用co-occurrence。 |
| 衝突2 | MethylTree、EPI-Clone、PCDH、Nomura 2026證明特定single-cell/barcode regime可用methylation trace lineage，所以不能blanket否定。 |
| 衝突3 | CNAqc數學上直接否定`read/VAF=CCF`；missing/mis-fit CN也不能當neutral。 |
| 衝突4 | 目前partial-read實作否定「所有determined trees都是full-molecule direct observation」。 |

## 本輪核心新來源

- TreeClone 2019：physical mutation-pair prior art與novelty邊界。
- PhyClone 2025：最新bulk CN/purity-aware Bayesian clone-tree comparator。
- CNAqc 2024：VAF/CCF/CN/multiplicity數學與QC錨點。
- ASMS 2026 VOR：no-phasing read clustering + 5,000-permutation significance prior art。
- Nomura 2026：longitudinal single-nucleus methylation/RNA phylogeny與orthogonal design。
- Lee 2026：single-cell-derived tumoroid WGS phylogeny，適合作biological validation target。

## 下一步 gate

### P0

- 完成fresh 7/7 canonical、sensitivity、backbone comparison、independent verifier、immutable receipt。
- 凍結clean commit/tag與完整caller/patch provenance。

### P1

- DeepSomatic second-caller + alternate phaser。
- MAPQ/BaseQ/MINREAD/MAX_SNV、region-gap/read-link與full-span connectivity敏感度。
- bounded-loss/error-aware solver；strict model保留baseline。
- CNAqc-style purity/CN/multiplicity QC；未核定`smcnbed未列=neutral`前保持unavailable。

### P2

- single-cell-derived colony/tumoroid、multi-region或time-point truth，逐mutation-state/node/edge驗證。
- patient cohort與low-purity/FFPE generalisation。
- methylation matched-normal、germline/imprinting、cell-type、CN/LOH、5mC/5hmC negative-control battery。

## 記憶與搜尋規則

1. 先查`/big7_disk/liaoyoyo2001/external_validation/_schema/paper_claims.tsv`，pending claim不能由舊narrative解鎖。
2. 再查semantic catalog/mapping，最後回讀CONTEXT與一手paper/repo。
3. live run百分比不進長期記憶；只記run root、state、gate與`_SUCCESS`。
4. 外部paper可`calibrate/complement/tension`內部結果，不能替本地producer、solver或biological truth作`supports`偷渡。

## 可寫／不可寫

可寫：`candidate regional mutation-state trees within parental haplotype families`、`partial-read-constrained`、`all minimum compatible solutions under a recurrence-free model`、`explicit refusal`。

不可寫：`biologically validated clone tree`、`first mutation-co-occurrence phylogeny`、`HP=clone`、`read AF=CCF`、`missing CN=neutral`、`methylation cannot trace lineage`、`methylation validates clones/variants`。
