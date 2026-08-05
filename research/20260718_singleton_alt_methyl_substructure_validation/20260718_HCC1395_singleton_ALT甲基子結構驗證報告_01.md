<!--
建立時間: 2026-07-19T00:49:47+08:00
目標: 完整驗證 HCC1395 positional-singleton focal-ALT read 甲基子結構
處理範圍: 7 datasets / 50,432 singleton sites；HCC1395 8,279 sites 全量；兩個 M2 PASS loci read-level exact join
關聯檔案:
  - InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/artifact.json
  - InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/validation_receipt.json
-->

# HCC1395 positional-singleton ALT 內甲基子結構完整驗證

> **結論：** HCC1395 8,279 個 singleton sites 中，M1 stable multigroup 為
> **734（8.87%）**；通過八個 measured-axis guardrails 的 M2 PASS 為
> **2（0.0242%）**。兩點可稱 read-level residual epigenetic partition candidate；
> G1/G2/formal R1、matched-normal 與 CN/CCF 未執行，因此 **不能稱已證實 clone/subclone**。

## 分母

- HCC1395 chr1-22 autosomal biallelic LongPhase-S PASS sSNV：79,687。
- 50 kb positional components：16,501；其中 singleton components/sites：8,279，多位點 components：8,222。
- singleton 真值層：TP 7,242、FP 153、UNASSESSED 884。

## 核心數字

- M1 evaluable：8,074/8,279（97.52%）。
- M1 stable multigroup：734/8,279（8.87%；95% Wilson CI
  8.27%–9.50%）。
- M2 PASS：2/8,279（0.0242%；95% Wilson CI
  0.0066%–0.0880%）。
- M2 FAIL：0；NOT_EVALUABLE：732；NOT_RUN：7,545。
- Established cellular clone/subclone：0，語意是必要驗證未跑，不是真陰性。

## 兩個 M2 PASS 例子

| Locus | Core ALT reads | A/B | Between median | Pooled within | Ratio | Methyl Δ B−A | Interpretation |
|---|---:|---:|---:|---:|---:|---:|---|
| chr14:86,272,476 A>T | 108 | 86/22 | 0.445 | 0.235 | 1.89× | -0.345 | global high/low methylation partition |
| chr22:47,466,517 A>G | 109 | 88/21 | 0.366 | 0.234 | 1.57× | -0.028 | CpG-pattern partition with similar global methylation means |


兩點所有 core reads 的 latest HP/PS 都是同一值，因此 HP 沒有把兩群切開，但也不能提供獨立 clone 佐證。
caller AF（chr14 0.827；chr22 1.000）只作 focal allele burden context，未用於甲基分群。

## Claim ceiling

- 可用：`M1 stable focal-ALT methyl multigroup`。
- 兩個清楚例子可用：`M2 read-level residual epigenetic partition candidate`。
- 不可用：confirmed clone/subclone、clone number、parent-child、linear/branching ancestry、unique true tree。

## 可重現輸出

- 完整 HCC1395 8,279-row audit：`InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/hcc1395_singleton_site_audit.tsv.gz`
- 驗證 receipt：`InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/results/validation_receipt.json`
- Canonical report artifact：`InterSubMod/research/20260718_singleton_alt_methyl_substructure_validation/artifact.json`
